#pragma once

// Ports of VariationUtils
// and small static helpers of VariationRealigner used by several stages.

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

#include "cli/config.hpp"
#include "core/scope.hpp"
#include "core/variation.hpp"

namespace omnivar::core {

using RefMap = std::unordered_map<int, char>;

// Map<Integer, Character>.get: null if absent.
[[nodiscard]] inline std::optional<char> ref_get(const RefMap& ref, int pos) {
    auto it = ref.find(pos);
    if (it == ref.end()) {
        return std::nullopt;
    }
    return it->second;
}

// ref.get(pos).toString() / unboxing: NullPointerException if absent.
[[nodiscard]] char ref_at(const RefMap& ref, int pos);

// VariationUtils.isEquals / isNotEquals (null-safe Character comparison)
[[nodiscard]] inline bool is_equals(std::optional<char> a, std::optional<char> b) { return a == b; }
[[nodiscard]] inline bool is_not_equals(std::optional<char> a, std::optional<char> b) { return a != b; }

// isHasAndEquals(char, ref, index): ref has a base at index and it equals ch.
[[nodiscard]] inline bool is_has_and_equals(char ch, const RefMap& ref, int index) {
    auto r = ref_get(ref, index);
    return r && *r == ch;
}
// isHasAndNotEquals(Character, ref, index)
[[nodiscard]] inline bool is_has_and_not_equals(char ch, const RefMap& ref, int index) {
    auto r = ref_get(ref, index);
    return r && *r != ch;
}
// isHasAndEquals(ref, index1, str, index2): str.charAt only evaluated when the ref base exists.
[[nodiscard]] bool is_has_and_equals(const RefMap& ref, int index1, std::string_view str, long index2);
[[nodiscard]] bool is_has_and_not_equals(const RefMap& ref, int index1, std::string_view str, long index2);

[[nodiscard]] int strand_bias(int forward_count, int reverse_count, const Config& conf);

// VariationUtils.findconseq; dir 3/5 checks adaptor seeds, 0 does not.
[[nodiscard]] std::string findconseq(Sclip& soft_clip, int dir, const RunContext& run, const Config& conf);

// VariationUtils.joinRef(ref, from, to): bases from..to inclusive, skipping missing positions.
[[nodiscard]] std::string join_ref(const RefMap& ref, int from, int to);

// VariationUtils.adjCnt(varToAdd, variant[, referenceVar]) and correctCnt.
void adj_cnt(Variation& to_add, const Variation& variant, Variation* reference_var = nullptr);
void correct_cnt(Variation& v);

[[nodiscard]] bool islowcomplexseq(std::string_view seq);

[[nodiscard]] bool ismatchref(std::string_view sequence, const RefMap& ref, int position, int dir, int mm_limit = 3);

}  // namespace omnivar::core
