#include "util/hash_order.hpp"

namespace omnivar::util {

// Model of the JDK 8+ HashMap putVal/resize for distinct keys without removals.
std::optional<std::vector<std::size_t>> hashmap_order(const std::vector<std::int32_t>& hashes) {
    constexpr std::size_t kDefaultCapacity = 16;
    constexpr std::size_t kTreeifyThreshold = 8;
    constexpr std::size_t kMinTreeifyCapacity = 64;

    // HashMap.hash(): h ^ (h >>> 16)
    auto spread = [](std::int32_t h) {
        const auto u = static_cast<std::uint32_t>(h);
        return u ^ (u >> 16);
    };

    std::vector<std::vector<std::size_t>> table;  // empty until first put (lazy allocation)
    std::size_t threshold = 0;
    std::size_t size = 0;

    auto resize = [&]() {
        const std::size_t old_cap = table.size();
        const std::size_t new_cap = old_cap == 0 ? kDefaultCapacity : old_cap * 2;
        std::vector<std::vector<std::size_t>> next(new_cap);
        for (std::size_t j = 0; j < old_cap; ++j) {
            for (std::size_t idx : table[j]) {
                // preserves relative order: lo list stays at j, hi list moves to j + old_cap
                const std::size_t b = (spread(hashes[idx]) & (new_cap - 1));
                next[b].push_back(idx);
            }
        }
        table = std::move(next);
        threshold = new_cap * 3 / 4;
    };

    bool treeified = false;
    for (std::size_t idx = 0; idx < hashes.size(); ++idx) {
        if (table.empty()) {
            resize();
        }
        auto& bin = table[spread(hashes[idx]) & (table.size() - 1)];
        bin.push_back(idx);
        // putVal: binCount >= TREEIFY_THRESHOLD - 1 when the bin now holds >= 9 nodes
        if (bin.size() >= kTreeifyThreshold + 1) {
            if (table.size() < kMinTreeifyCapacity) {
                resize();  // treeifyBin resizes small tables instead of treeifying
            } else {
                treeified = true;
            }
        }
        if (++size > threshold) {
            resize();
        }
    }
    if (treeified) {
        return std::nullopt;
    }
    std::vector<std::size_t> order;
    order.reserve(hashes.size());
    for (const auto& bin : table) {
        order.insert(order.end(), bin.begin(), bin.end());
    }
    return order;
}

}  // namespace omnivar::util
