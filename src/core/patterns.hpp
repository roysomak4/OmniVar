#pragma once

// Hand-written equivalents of the upstream regex patterns (Patterns) used by the counting
// core. Each function documents the regex it replaces; "find" means Matcher.find() (leftmost
// match). Strings here never contain line terminators, so '.' is any character.

#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace omnivar::core::patterns {

inline bool is_digit(char c) { return c >= '0' && c <= '9'; }
inline bool is_atgc(char c) { return c == 'A' || c == 'T' || c == 'G' || c == 'C'; }
inline bool is_atgnc(char c) { return is_atgc(c) || c == 'N'; }

// Length of the run of digits starting at i.
inline std::size_t digit_run(std::string_view s, std::size_t i) {
    std::size_t j = i;
    while (j < s.size() && is_digit(s[j])) {
        ++j;
    }
    return j - i;
}

// ANY_SV "<(...)>" find: group 1 (three characters).
inline std::optional<std::string> any_sv(std::string_view s) {
    for (std::size_t i = 0; i + 4 < s.size(); ++i) {
        if (s[i] == '<' && s[i + 4] == '>') {
            return std::string(s.substr(i + 1, 3));
        }
    }
    return std::nullopt;
}

// "<tag(\d+)" find, e.g. DUP_NUM "<dup(\d+)", INV_NUM "<inv(\d+)": group 1 (digits).
inline std::optional<std::string> tag_num(std::string_view s, std::string_view tag) {
    for (std::size_t p = s.find(tag); p != std::string_view::npos; p = s.find(tag, p + 1)) {
        const std::size_t d = digit_run(s, p + tag.size());
        if (d > 0) {
            return std::string(s.substr(p + tag.size(), d));
        }
    }
    return std::nullopt;
}
inline std::optional<std::string> dup_num(std::string_view s) { return tag_num(s, "<dup"); }
inline std::optional<std::string> inv_num(std::string_view s) { return tag_num(s, "<inv"); }

// SOME_SV_NUMBERS "<(...)\d+>" find
inline bool some_sv_numbers(std::string_view s) {
    for (std::size_t i = 0; i + 5 < s.size(); ++i) {
        if (s[i] != '<') {
            continue;
        }
        const std::size_t d = digit_run(s, i + 4);
        if (d > 0 && i + 4 + d < s.size() && s[i + 4 + d] == '>') {
            return true;
        }
    }
    return false;
}

// BEGIN_MINUS_NUMBER "^-(\d+)": group 1
inline std::optional<std::string> begin_minus_number(std::string_view s) {
    if (s.empty() || s[0] != '-') {
        return std::nullopt;
    }
    const std::size_t d = digit_run(s, 1);
    if (d == 0) {
        return std::nullopt;
    }
    return std::string(s.substr(1, d));
}

// BEGIN_MINUS_NUMBER_CARET "^-\d+\^"
inline bool begin_minus_number_caret(std::string_view s) {
    if (s.empty() || s[0] != '-') {
        return false;
    }
    const std::size_t d = digit_run(s, 1);
    return d > 0 && 1 + d < s.size() && s[1 + d] == '^';
}

// BEGIN_DIGITS "^(\d+)": group 1
inline std::optional<std::string> begin_digits(std::string_view s) {
    const std::size_t d = digit_run(s, 0);
    if (d == 0) {
        return std::nullopt;
    }
    return std::string(s.substr(0, d));
}

// "&([ATGC]+)" (AMP_ATGC) find: group 1
inline std::optional<std::string> amp_atgc(std::string_view s) {
    for (std::size_t p = s.find('&'); p != std::string_view::npos; p = s.find('&', p + 1)) {
        std::size_t j = p + 1;
        while (j < s.size() && is_atgc(s[j])) {
            ++j;
        }
        if (j > p + 1) {
            return std::string(s.substr(p + 1, j - p - 1));
        }
    }
    return std::nullopt;
}

// CARET_ATGNC "\^([ATGNC]+)" find
inline bool caret_atgnc(std::string_view s) {
    for (std::size_t p = s.find('^'); p != std::string_view::npos; p = s.find('^', p + 1)) {
        if (p + 1 < s.size() && is_atgnc(s[p + 1])) {
            return true;
        }
    }
    return false;
}

// HASH_GROUP_CARET_GROUP "#(.+)\^(.+)" find: the leftmost '#' that has a '^' at least two
// characters later with at least one character after it; greedy group 1 extends to the last
// such '^'.
inline std::optional<std::pair<std::string, std::string>> hash_group_caret_group(std::string_view s) {
    for (std::size_t h = s.find('#'); h != std::string_view::npos; h = s.find('#', h + 1)) {
        if (s.size() < h + 4) {
            return std::nullopt;
        }
        for (std::size_t p = s.size() - 2; p >= h + 2; --p) {
            if (s[p] == '^') {
                return std::pair{std::string(s.substr(h + 1, p - h - 1)), std::string(s.substr(p + 1))};
            }
        }
    }
    return std::nullopt;
}

// BEGIN_ATGC_END "^[ATGC]+$"
inline bool begin_atgc_end(std::string_view s) {
    if (s.empty()) {
        return false;
    }
    for (char c : s) {
        if (!is_atgc(c)) {
            return false;
        }
    }
    return true;
}

// B_A7 "^.AAAAAAA" / B_T7 "^.TTTTTTT"
inline bool begin_any_7(std::string_view s, char base) {
    if (s.size() < 8) {
        return false;
    }
    for (std::size_t i = 1; i < 8; ++i) {
        if (s[i] != base) {
            return false;
        }
    }
    return true;
}

// String.replaceFirst("^-\\d+", "")
inline std::string remove_begin_minus_number(std::string_view s) {
    if (!s.empty() && s[0] == '-') {
        const std::size_t d = digit_run(s, 1);
        if (d > 0) {
            return std::string(s.substr(1 + d));
        }
    }
    return std::string(s);
}

// String.replaceFirst("\\^(\\d+)?", ""): first '^' and the digits right after it.
inline std::string remove_first_caret_digits(std::string_view s) {
    const auto p = s.find('^');
    if (p == std::string_view::npos) {
        return std::string(s);
    }
    const std::size_t d = digit_run(s, p + 1);
    return std::string(s.substr(0, p)) + std::string(s.substr(p + 1 + d));
}

}  // namespace omnivar::core::patterns
