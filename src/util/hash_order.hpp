#pragma once

// String/Integer hash codes and HashMap iteration order, for loops whose order is observable.
//
// hashmap_order(hashes) returns the iteration order of a HashMap (default capacity 16,
// load factor 0.75) after inserting distinct keys with the given hashCode()s in the given order,
// with no removals. It models HashMap.putVal / resize (bins keep insertion order, split into
// lo/hi halves on resize, and a bin reaching 9 entries forces a resize while the table is below
// 64). It does not model tree bins: if a bin would be treeified (>= 9 entries with table >= 64)
// the order is not reproducible and std::nullopt is returned.
//
// Verified against the JVM: tests/unit/data/hashmap_order_vectors.txt.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace senri::util {

// String.hashCode() of an ASCII/Latin-1 string (bytes taken as UTF-16 code units).
[[nodiscard]] constexpr std::int32_t string_hash_code(std::string_view s) {
    std::uint32_t h = 0;
    for (char c : s) {
        h = 31 * h + static_cast<unsigned char>(c);
    }
    return static_cast<std::int32_t>(h);
}

// Integer.hashCode()
[[nodiscard]] constexpr std::int32_t integer_hash_code(std::int32_t v) { return v; }

// Indices into `hashes` in HashMap iteration order (see above).
[[nodiscard]] std::optional<std::vector<std::size_t>> hashmap_order(const std::vector<std::int32_t>& hashes);

}  // namespace senri::util
