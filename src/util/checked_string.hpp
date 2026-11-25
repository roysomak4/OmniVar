#pragma once

// String operations with bounds checks. The caller catches exceptions per read and per position
// (print_exception_and_continue) and keeps whatever that read had already counted, so it has to
// fail at the same points as upstream: these helpers throw RecoverableError on an out-of-range
// index or a missing value.

#include <stdexcept>
#include <string>
#include <string_view>

namespace omnivar::util {

struct RecoverableError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

[[noreturn]] inline void throw_index(long index, std::size_t length) {
    throw RecoverableError("StringIndexOutOfBoundsException: index " + std::to_string(index) + ", length " +
                           std::to_string(length));
}
[[noreturn]] inline void throw_npe(std::string_view what) {
    throw RecoverableError("NullPointerException: " + std::string(what));
}

// String.charAt
[[nodiscard]] inline char char_at(std::string_view s, long index) {
    if (index < 0 || static_cast<std::size_t>(index) >= s.size()) {
        throw_index(index, s.size());
    }
    return s[static_cast<std::size_t>(index)];
}

// Base quality at index (charAt(i) - 33), chars treated as unsigned (16-bit char semantics).
[[nodiscard]] inline int qual_at(std::string_view q, long index) {
    return static_cast<unsigned char>(char_at(q, index)) - 33;
}

// String.substring(begin, end)
[[nodiscard]] inline std::string substring(std::string_view s, long begin, long end) {
    if (begin < 0 || end > static_cast<long>(s.size()) || begin > end) {
        throw RecoverableError("StringIndexOutOfBoundsException: begin " + std::to_string(begin) + ", end " +
                               std::to_string(end) + ", length " + std::to_string(s.size()));
    }
    return std::string(s.substr(static_cast<std::size_t>(begin), static_cast<std::size_t>(end - begin)));
}
// String.substring(begin)
[[nodiscard]] inline std::string substring(std::string_view s, long begin) {
    return substring(s, begin, static_cast<long>(s.size()));
}

[[nodiscard]] inline std::string substr(std::string_view s, long idx) {
    const long len = static_cast<long>(s.size());
    if (idx >= 0) {
        return substring(s, std::min(len, idx));
    }
    return substring(s, std::max(0L, len + idx));
}

[[nodiscard]] inline std::string substr(std::string_view s, long begin, long len) {
    const long n = static_cast<long>(s.size());
    if (begin < 0) {
        begin = n + begin;
    }
    if (len > 0) {
        return substring(s, begin, std::min(begin + len, n));
    }
    if (len == 0) {
        return "";
    }
    const long end = n + len;
    if (end < begin) {
        return "";
    }
    return substring(s, begin, end);
}

// Utils.charAt(String, int): negative index counts from the end; out of range below 0 -> (char)-1.
[[nodiscard]] inline int utils_char_at(std::string_view s, long index) {
    if (index < 0) {
        const long i = static_cast<long>(s.size()) + index;
        if (i < 0) {
            return 0xFFFF;  // (char) -1
        }
        return static_cast<unsigned char>(s[static_cast<std::size_t>(i)]);
    }
    return static_cast<unsigned char>(char_at(s, index));
}

// String.replaceFirst with a pattern that is a literal string (e.g. "&", "#", "\\^").
[[nodiscard]] inline std::string replace_first(std::string_view s, std::string_view literal, std::string_view with) {
    const auto p = s.find(literal);
    if (p == std::string_view::npos) {
        return std::string(s);
    }
    std::string out(s.substr(0, p));
    out += with;
    out += s.substr(p + literal.size());
    return out;
}

// String.replace(CharSequence, CharSequence): all occurrences.
[[nodiscard]] inline std::string replace_all(std::string_view s, std::string_view literal, std::string_view with) {
    std::string out;
    std::size_t i = 0;
    while (true) {
        const auto p = s.find(literal, i);
        if (p == std::string_view::npos) {
            out += s.substr(i);
            return out;
        }
        out += s.substr(i, p - i);
        out += with;
        i = p + literal.size();
    }
}

}  // namespace omnivar::util
