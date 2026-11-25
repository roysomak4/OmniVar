#pragma once

// A small backtracking regular-expression engine with JDK regex semantics (the dialect of
// the upstream sample and delimiter options) for the subset needed: sample-name patterns (Patterns.SAMPLE_PATTERN,
// SAMPLE_PATTERN2, user -n), and String.split delimiters (user -d). std::regex is not used (slow, differs across
// libraries).
//
// Supported: literals; '.' (any char except \n and \r); ^ and $ (no MULTILINE: $ also matches
// before a final line terminator); character classes [...] with negation, ranges, and the
// escapes below; \d \D \w \W \s \S \t \n \r \f \e \a and backslash + non-alphanumeric
// (literal); groups (...) and (?:...); alternation |; quantifiers * + ? {n} {n,} {n,m}, greedy
// or lazy (?). Anything else (backreferences, lookaround, possessive quantifiers, flags,
// \b, \p{..}, class intersections) is rejected at compile time.
//
// Matching is on bytes. The upstream dialect matches UTF-16 chars, so '.' differs only next to non-ASCII text.

#include <bitset>
#include <cstddef>
#include <expected>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace omnivar::util {

class Regex {
public:
    struct Match {
        std::size_t begin = 0;
        std::size_t end = 0;
        // groups[0] is the whole match; unmatched groups are nullopt.
        std::vector<std::optional<std::pair<std::size_t, std::size_t>>> groups;

        [[nodiscard]] std::optional<std::string_view> group(std::string_view text, std::size_t i) const;
    };

    [[nodiscard]] static std::expected<Regex, std::string> compile(std::string_view pattern);

    // Matcher.find() starting at `from`: leftmost match, upstream backtracking preference order.
    [[nodiscard]] std::optional<Match> find(std::string_view text, std::size_t from = 0) const;

    [[nodiscard]] std::size_t group_count() const { return groups_; }
    [[nodiscard]] const std::string& pattern() const { return pattern_; }

    struct Node;  // implementation detail

    Regex(Regex&&) noexcept;
    Regex& operator=(Regex&&) noexcept;
    Regex(const Regex&) = delete;
    Regex& operator=(const Regex&) = delete;
    ~Regex();

private:
    Regex() = default;
    std::string pattern_;
    std::size_t groups_ = 0;  // number of capturing groups
    std::unique_ptr<Node> root_;
    // Precomputed at compile(): every match starts at 0; characters every match contains.
    bool anchored_ = false;
    std::bitset<256> required_;
};

// String.split(regex) with limit 0: leading empty strings kept (except for a zero-width match at
// index 0), trailing empty strings removed; no match -> the whole input as one element.
// The single-character fast path of String.split is the same as compiling the regex.
[[nodiscard]] std::vector<std::string_view> regex_split(std::string_view text, const Regex& re);
[[nodiscard]] std::expected<std::vector<std::string_view>, std::string> regex_split(std::string_view text,
                                                                                    std::string_view regex);

}  // namespace omnivar::util
