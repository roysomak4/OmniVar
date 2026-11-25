#pragma once

// Containers with JDK iteration order, for every place where upstream iterates a map
// and the order can change results (or the order of output lines):
//
//   HashOrderMap<K, V>   HashMap / HashSet: buckets by (h ^ h >>> 16) & (capacity - 1),
//                       insertion order within a bucket, capacity doubling at load factor 0.75,
//                       lo/hi split on resize, and the small-table resize that treeifyBin does.
//                       Tree bins (>= 9 colliding keys with capacity >= 64) are not modelled:
//                       such a bin keeps insertion order and tree_bins() reports it.
//   LinkedMap<K, V>     LinkedHashMap (insertion order, no access order).
//
// Values have stable addresses, like object references, and survive erase() (a pointer
// taken before the removal keeps working, as an object reference would). Keys: int, char,
// std::string (hashCode as on the JVM). Verified against the JVM by tests/unit/test_stats_containers.cpp.

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "util/hash_order.hpp"

namespace omnivar::util {

inline std::int32_t hash_code(int v) { return v; }
inline std::int32_t hash_code(char c) { return static_cast<unsigned char>(c); }  // Character.hashCode
inline std::int32_t hash_code(const std::string& s) { return string_hash_code(s); }

template <class K, class V>
class HashOrderMap {
    struct Node {
        K key;
        V value;
        std::uint32_t hash;
        std::int32_t next;
        bool alive;
    };

public:
    HashOrderMap() = default;

    [[nodiscard]] std::size_t size() const { return size_; }
    [[nodiscard]] bool empty() const { return size_ == 0; }
    [[nodiscard]] bool tree_bins() const { return tree_bins_; }

    [[nodiscard]] V* get(const K& key) {
        const std::int32_t n = find(key);
        return n < 0 ? nullptr : &nodes_[static_cast<std::size_t>(n)].value;
    }
    [[nodiscard]] const V* get(const K& key) const {
        const std::int32_t n = find(key);
        return n < 0 ? nullptr : &nodes_[static_cast<std::size_t>(n)].value;
    }
    [[nodiscard]] bool contains(const K& key) const { return find(key) >= 0; }

    // HashMap.put: replaces the value of an existing key (keeping its position), else inserts.
    V& put(const K& key, V value) {
        const std::int32_t n = find(key);
        if (n >= 0) {
            nodes_[static_cast<std::size_t>(n)].value = std::move(value);
            return nodes_[static_cast<std::size_t>(n)].value;
        }
        return insert_new(key, std::move(value));
    }

    // get, inserting a default-constructed value if absent (the getOrElse/computeIfAbsent idiom).
    V& get_or_create(const K& key) {
        const std::int32_t n = find(key);
        if (n >= 0) {
            return nodes_[static_cast<std::size_t>(n)].value;
        }
        return insert_new(key, V{});
    }

    // HashMap.remove
    bool erase(const K& key) {
        if (table_.empty()) {
            return false;
        }
        const std::uint32_t h = spread(key);
        const std::size_t b = h & (table_.size() - 1);
        std::int32_t prev = -1;
        for (std::int32_t n = table_[b]; n >= 0; n = nodes_[static_cast<std::size_t>(n)].next) {
            Node& node = nodes_[static_cast<std::size_t>(n)];
            if (node.hash == h && node.key == key) {
                if (prev < 0) {
                    table_[b] = node.next;
                } else {
                    nodes_[static_cast<std::size_t>(prev)].next = node.next;
                }
                if (tails_[b] == n) {
                    tails_[b] = prev;
                }
                // Like an object reference, a value obtained before remove() stays valid
                // (the realigner keeps using variations it has removed from their map).
                node.alive = false;
                --size_;
                return true;
            }
            prev = n;
        }
        return false;
    }

    // Visits entries in JDK iteration order: f(const K&, V&). The map must not be structurally
    // modified during the visit (upstream would throw ConcurrentModificationException).
    template <class F>
    void for_each(F&& f) {
        for (std::int32_t head : table_) {
            for (std::int32_t n = head; n >= 0;) {
                Node& node = nodes_[static_cast<std::size_t>(n)];
                const std::int32_t next = node.next;
                f(static_cast<const K&>(node.key), node.value);
                n = next;
            }
        }
    }
    template <class F>
    void for_each(F&& f) const {
        for (std::int32_t head : table_) {
            for (std::int32_t n = head; n >= 0; n = nodes_[static_cast<std::size_t>(n)].next) {
                const Node& node = nodes_[static_cast<std::size_t>(n)];
                f(node.key, node.value);
            }
        }
    }

    // Keys in iteration order (new ArrayList<>(map.keySet())).
    [[nodiscard]] std::vector<K> keys() const {
        std::vector<K> out;
        out.reserve(size_);
        for_each([&](const K& k, const V&) { out.push_back(k); });
        return out;
    }

private:
    static constexpr std::size_t kDefaultCapacity = 16;
    static constexpr std::size_t kTreeifyThreshold = 8;
    static constexpr std::size_t kMinTreeifyCapacity = 64;

    static std::uint32_t spread(const K& key) {
        const auto u = static_cast<std::uint32_t>(hash_code(key));
        return u ^ (u >> 16);
    }

    std::int32_t find(const K& key) const {
        if (table_.empty()) {
            return -1;
        }
        const std::uint32_t h = spread(key);
        for (std::int32_t n = table_[h & (table_.size() - 1)]; n >= 0; n = nodes_[static_cast<std::size_t>(n)].next) {
            const Node& node = nodes_[static_cast<std::size_t>(n)];
            if (node.hash == h && node.key == key) {
                return n;
            }
        }
        return -1;
    }

    V& insert_new(const K& key, V value) {
        if (table_.empty()) {
            resize();
        }
        const std::uint32_t h = spread(key);
        const auto idx = static_cast<std::int32_t>(nodes_.size());
        nodes_.push_back(Node{key, std::move(value), h, -1, true});
        const std::size_t b = h & (table_.size() - 1);
        std::size_t chain = 1;
        if (table_[b] < 0) {
            table_[b] = idx;
        } else {
            for (std::int32_t n = table_[b]; n >= 0; n = nodes_[static_cast<std::size_t>(n)].next) {
                ++chain;
            }
            nodes_[static_cast<std::size_t>(tails_[b])].next = idx;
        }
        tails_[b] = idx;
        // HashMap.putVal: a bin that reaches TREEIFY_THRESHOLD + 1 nodes calls treeifyBin, which
        // resizes tables smaller than MIN_TREEIFY_CAPACITY instead of making a tree bin.
        if (chain >= kTreeifyThreshold + 1) {
            if (table_.size() < kMinTreeifyCapacity) {
                resize();
            } else {
                tree_bins_ = true;
            }
        }
        if (++size_ > threshold_) {
            resize();
        }
        return nodes_[static_cast<std::size_t>(idx)].value;
    }

    // HashMap.resize: doubles the table; each bin splits into lo/hi lists in order, which is
    // the same as re-appending nodes in old iteration order.
    void resize() {
        const std::size_t new_cap = table_.empty() ? kDefaultCapacity : table_.size() * 2;
        std::vector<std::int32_t> heads(new_cap, -1);
        std::vector<std::int32_t> tails(new_cap, -1);
        for (std::int32_t head : table_) {
            for (std::int32_t n = head; n >= 0;) {
                Node& node = nodes_[static_cast<std::size_t>(n)];
                const std::int32_t next = node.next;
                node.next = -1;
                const std::size_t b = node.hash & (new_cap - 1);
                if (heads[b] < 0) {
                    heads[b] = n;
                } else {
                    nodes_[static_cast<std::size_t>(tails[b])].next = n;
                }
                tails[b] = n;
                n = next;
            }
        }
        table_ = std::move(heads);
        tails_ = std::move(tails);
        threshold_ = new_cap * 3 / 4;
    }

    std::deque<Node> nodes_;
    std::vector<std::int32_t> table_;
    std::vector<std::int32_t> tails_;
    std::size_t size_ = 0;
    std::size_t threshold_ = 0;
    bool tree_bins_ = false;
};

// LinkedHashMap with insertion order. Values live on the heap so that rename() can
// move one to another key without changing its address (upstream: put(newKey, map.get(oldKey));
// remove(oldKey), where references to the object stay valid). Small maps (the usual 1-4 alleles
// at a position) are searched linearly; a hash index is built once a map exceeds kIndexAt live
// entries. Either way the iteration order is the insertion order.
template <class K, class V>
class LinkedMap {
    struct Node {
        K key;
        std::unique_ptr<V> value;
        bool alive;
    };
    static constexpr std::size_t kIndexAt = 8;

public:
    LinkedMap() = default;
    LinkedMap(LinkedMap&&) noexcept = default;
    LinkedMap& operator=(LinkedMap&&) noexcept = default;
    // A copy holds new objects (upstream: new LinkedHashMap<>(other) copies references, but no
    // caller relies on sharing the values of a copied map).
    LinkedMap(const LinkedMap& other) { *this = other; }
    LinkedMap& operator=(const LinkedMap& other) {
        if (this != &other) {
            nodes_.clear();
            index_.clear();
            size_ = 0;
            other.for_each([&](const K& k, const V& v) { put(k, v); });
        }
        return *this;
    }
    ~LinkedMap() = default;

    [[nodiscard]] std::size_t size() const { return size_; }
    [[nodiscard]] bool empty() const { return size_ == 0; }

    [[nodiscard]] V* get(const K& key) {
        const std::size_t n = find(key);
        return n == kNone ? nullptr : nodes_[n].value.get();
    }
    [[nodiscard]] const V* get(const K& key) const {
        const std::size_t n = find(key);
        return n == kNone ? nullptr : nodes_[n].value.get();
    }
    [[nodiscard]] bool contains(const K& key) const { return find(key) != kNone; }

    V& put(const K& key, V value) {
        if (V* v = get(key)) {
            *v = std::move(value);
            return *v;
        }
        return append(key, std::make_unique<V>(std::move(value)));
    }
    V& get_or_create(const K& key) {
        if (V* v = get(key)) {
            return *v;
        }
        return append(key, std::make_unique<V>());
    }
    // remove(key); the value stays valid for pointers taken before (see above).
    bool erase(const K& key) {
        const std::size_t n = find(key);
        if (n == kNone) {
            return false;
        }
        nodes_[n].alive = false;
        if (!index_.empty()) {
            index_.erase(key);
        }
        --size_;
        return true;
    }
    // put(to, get(from)); remove(from): the same object, now under `to`. When `to` exists its
    // entry keeps its position and takes the object; otherwise `to` is appended. False (and no
    // change) when `from` is absent.
    bool rename(const K& from, const K& to) {
        const std::size_t src = find(from);
        if (src == kNone) {
            return false;
        }
        if (from == to) {
            return true;
        }
        auto object = std::make_unique<V>();
        std::swap(object, nodes_[src].value);  // the source keeps a fresh (dead) value
        if (const std::size_t dst = find(to); dst != kNone) {
            nodes_[dst].value = std::move(object);
        } else {
            append(to, std::move(object));
        }
        erase(from);
        return true;
    }

    template <class F>
    void for_each(F&& f) {
        for (Node& n : nodes_) {
            if (n.alive) {
                f(static_cast<const K&>(n.key), *n.value);
            }
        }
    }
    template <class F>
    void for_each(F&& f) const {
        for (const Node& n : nodes_) {
            if (n.alive) {
                f(n.key, static_cast<const V&>(*n.value));
            }
        }
    }
    [[nodiscard]] std::vector<K> keys() const {
        std::vector<K> out;
        out.reserve(size());
        for_each([&](const K& k, const V&) { out.push_back(k); });
        return out;
    }

private:
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);
    struct Hash {
        std::size_t operator()(const K& k) const { return std::hash<K>{}(k); }
    };

    std::size_t find(const K& key) const {
        if (!index_.empty()) {
            auto it = index_.find(key);
            return it == index_.end() ? kNone : it->second;
        }
        for (std::size_t i = 0; i < nodes_.size(); ++i) {
            if (nodes_[i].alive && nodes_[i].key == key) {
                return i;
            }
        }
        return kNone;
    }

    V& append(const K& key, std::unique_ptr<V> value) {
        nodes_.push_back(Node{key, std::move(value), true});
        ++size_;
        if (!index_.empty()) {
            index_.emplace(key, nodes_.size() - 1);
        } else if (size_ > kIndexAt) {
            for (std::size_t i = 0; i < nodes_.size(); ++i) {
                if (nodes_[i].alive) {
                    index_.emplace(nodes_[i].key, i);
                }
            }
        }
        return *nodes_.back().value;
    }

    std::deque<Node> nodes_;
    std::unordered_map<K, std::size_t, Hash> index_;  // empty until size() > kIndexAt
    std::size_t size_ = 0;
};

}  // namespace omnivar::util
