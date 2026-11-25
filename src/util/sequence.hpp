#pragma once

// Base-sequence helpers matching htsjdk SequenceUtil and the upstream Utils.

#include <algorithm>
#include <string>
#include <string_view>

namespace omnivar::util {

// htsjdk SequenceUtil.complement(byte): A<->T, C<->G in both cases; anything else unchanged.
[[nodiscard]] constexpr char complement(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return b;
    }
}

// Utils.complement(String)
[[nodiscard]] inline std::string complement(std::string_view s) {
    std::string out(s);
    std::transform(out.begin(), out.end(), out.begin(), [](char c) { return complement(c); });
    return out;
}

// Utils.reverse(String)
[[nodiscard]] inline std::string reverse(std::string_view s) { return std::string(s.rbegin(), s.rend()); }

}  // namespace omnivar::util
