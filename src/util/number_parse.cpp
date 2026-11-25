#include "util/number_parse.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

namespace omnivar::util {

namespace {

std::string input_error(std::string_view s) { return "For input string: \"" + std::string(s) + "\""; }

// Digits in the given radix, accumulating as a negative number like Long.parseLong does so
// that MIN_VALUE parses. Returns false on any non-digit or overflow below `limit`.
bool accumulate(std::string_view digits, int radix, std::int64_t limit, std::int64_t& out) {
    if (digits.empty()) {
        return false;
    }
    const std::int64_t mult_min = limit / radix;
    std::int64_t result = 0;
    for (char c : digits) {
        int d;
        if (c >= '0' && c <= '9') {
            d = c - '0';
        } else if (c >= 'a' && c <= 'z') {
            d = c - 'a' + 10;
        } else if (c >= 'A' && c <= 'Z') {
            d = c - 'A' + 10;
        } else {
            return false;
        }
        if (d >= radix || result < mult_min) {
            return false;
        }
        result *= radix;
        if (result < limit + d) {
            return false;
        }
        result -= d;
    }
    out = result;
    return true;
}

// Long.parseLong(s, radix) / Integer.parseInt(s, radix) core.
std::expected<std::int64_t, std::string> parse_signed(std::string_view s, int radix, std::int64_t min_value,
                                                      std::int64_t max_value) {
    if (s.empty()) {
        return std::unexpected(input_error(s));
    }
    bool negative = false;
    std::string_view digits = s;
    if (digits[0] == '-' || digits[0] == '+') {
        negative = digits[0] == '-';
        digits.remove_prefix(1);
    }
    const std::int64_t limit = negative ? min_value : -max_value;
    std::int64_t neg = 0;
    if (!accumulate(digits, radix, limit, neg)) {
        return std::unexpected(input_error(s));
    }
    return negative ? neg : -neg;
}

}  // namespace

std::expected<std::int32_t, std::string> parse_int(std::string_view s) {
    auto r = parse_signed(s, 10, std::numeric_limits<std::int32_t>::min(), std::numeric_limits<std::int32_t>::max());
    if (!r) {
        return std::unexpected(r.error());
    }
    return static_cast<std::int32_t>(*r);
}

std::expected<std::int64_t, std::string> parse_long(std::string_view s) {
    return parse_signed(s, 10, std::numeric_limits<std::int64_t>::min(), std::numeric_limits<std::int64_t>::max());
}

std::expected<std::int32_t, std::string> decode_int(std::string_view nm) {
    if (nm.empty()) {
        return std::unexpected(std::string("Zero length string"));
    }
    std::size_t index = 0;
    bool negative = false;
    if (nm[0] == '-') {
        negative = true;
        ++index;
    } else if (nm[0] == '+') {
        ++index;
    }
    int radix = 10;
    std::string_view rest = nm.substr(index);
    if (rest.starts_with("0x") || rest.starts_with("0X")) {
        index += 2;
        radix = 16;
    } else if (rest.starts_with("#")) {
        index += 1;
        radix = 16;
    } else if (rest.starts_with("0") && rest.size() > 1) {
        index += 1;
        radix = 8;
    }
    std::string_view body = nm.substr(index);
    if (body.starts_with("-") || body.starts_with("+")) {
        return std::unexpected(std::string("Sign character in wrong position"));
    }
    // Parses the body (with the sign re-applied) in the radix; MIN_VALUE handled by parsing negative.
    std::int64_t neg = 0;
    const std::int64_t limit = negative ? std::numeric_limits<std::int32_t>::min()
                                        : -static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max());
    if (!accumulate(body, radix, limit, neg)) {
        return std::unexpected(input_error(body));
    }
    return static_cast<std::int32_t>(negative ? neg : -neg);
}

std::expected<double, std::string> parse_double(std::string_view in) {
    std::string_view s = in;
    while (!s.empty() && static_cast<unsigned char>(s.front()) <= ' ') {
        s.remove_prefix(1);
    }
    while (!s.empty() && static_cast<unsigned char>(s.back()) <= ' ') {
        s.remove_suffix(1);
    }
    if (s.empty()) {
        return std::unexpected(std::string("empty String"));
    }
    bool negative = false;
    std::string_view body = s;
    if (body[0] == '-' || body[0] == '+') {
        negative = body[0] == '-';
        body.remove_prefix(1);
    }
    if (body == "NaN") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (body == "Infinity") {
        return negative ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
    }
    // Optional type suffix.
    if (!body.empty() && (body.back() == 'f' || body.back() == 'F' || body.back() == 'd' || body.back() == 'D')) {
        // Not for hex digits: a hex literal must end with its binary exponent, so a trailing
        // d/f there is still a suffix (the exponent is decimal digits).
        body.remove_suffix(1);
    }
    if (body.empty()) {
        return std::unexpected(input_error(in));
    }
    const bool hex = body.starts_with("0x") || body.starts_with("0X");
    auto valid_decimal = [](std::string_view b) {
        // digits [. digits] [eE [+-] digits], at least one mantissa digit
        std::size_t i = 0;
        std::size_t mant = 0;
        while (i < b.size() && b[i] >= '0' && b[i] <= '9') {
            ++i;
            ++mant;
        }
        if (i < b.size() && b[i] == '.') {
            ++i;
            while (i < b.size() && b[i] >= '0' && b[i] <= '9') {
                ++i;
                ++mant;
            }
        }
        if (mant == 0) {
            return false;
        }
        if (i < b.size() && (b[i] == 'e' || b[i] == 'E')) {
            ++i;
            if (i < b.size() && (b[i] == '+' || b[i] == '-')) {
                ++i;
            }
            std::size_t exp = 0;
            while (i < b.size() && b[i] >= '0' && b[i] <= '9') {
                ++i;
                ++exp;
            }
            if (exp == 0) {
                return false;
            }
        }
        return i == b.size();
    };
    auto valid_hex = [](std::string_view b) {
        // hexdigits [. hexdigits] pP [+-] digits ; binary exponent required
        std::size_t i = 0;
        std::size_t mant = 0;
        auto is_hex = [](char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'); };
        while (i < b.size() && is_hex(b[i])) {
            ++i;
            ++mant;
        }
        if (i < b.size() && b[i] == '.') {
            ++i;
            while (i < b.size() && is_hex(b[i])) {
                ++i;
                ++mant;
            }
        }
        if (mant == 0 || i >= b.size() || (b[i] != 'p' && b[i] != 'P')) {
            return false;
        }
        ++i;
        if (i < b.size() && (b[i] == '+' || b[i] == '-')) {
            ++i;
        }
        std::size_t exp = 0;
        while (i < b.size() && b[i] >= '0' && b[i] <= '9') {
            ++i;
            ++exp;
        }
        return exp > 0 && i == b.size();
    };
    if (hex ? !valid_hex(body.substr(2)) : !valid_decimal(body)) {
        return std::unexpected(input_error(in));
    }
    // Floating-point std::from_chars is missing from libc++ < 20 (and AppleClang 16), so use
    // strtod: correctly rounded on glibc and macOS, accepts hex floats, and the program never
    // calls setlocale, so the decimal point is '.'. Overflow gives HUGE_VAL (= Infinity) and
    // underflow the correctly rounded subnormal or 0.
    const std::string copy(body);
    char* end = nullptr;
    const double value = std::strtod(copy.c_str(), &end);
    if (end != copy.c_str() + copy.size()) {
        return std::unexpected(input_error(in));
    }
    return negative ? -value : value;
}

std::int32_t d2i_cast(double d) {
    if (std::isnan(d)) {
        return 0;
    }
    if (d >= 2147483647.0) {
        return std::numeric_limits<std::int32_t>::max();
    }
    if (d <= -2147483648.0) {
        return std::numeric_limits<std::int32_t>::min();
    }
    return static_cast<std::int32_t>(d);
}

std::int64_t d2l_cast(double d) {
    if (std::isnan(d)) {
        return 0;
    }
    if (d >= 9223372036854775807.0) {
        return std::numeric_limits<std::int64_t>::max();
    }
    if (d <= -9223372036854775808.0) {
        return std::numeric_limits<std::int64_t>::min();
    }
    return static_cast<std::int64_t>(d);
}

std::int32_t ParsedNumber::int_value() const {
    if (const auto* l = std::get_if<std::int64_t>(&v_)) {
        return static_cast<std::int32_t>(static_cast<std::uint32_t>(static_cast<std::uint64_t>(*l)));
    }
    return d2i_cast(std::get<double>(v_));
}

double ParsedNumber::double_value() const {
    if (const auto* l = std::get_if<std::int64_t>(&v_)) {
        return static_cast<double>(*l);
    }
    return std::get<double>(v_);
}

std::expected<ParsedNumber, std::string> create_number(std::string_view s) {
    if (s.find('.') != std::string_view::npos) {
        auto d = parse_double(s);
        if (!d) {
            return std::unexpected(d.error());
        }
        return ParsedNumber(*d);
    }
    auto l = parse_long(s);
    if (!l) {
        return std::unexpected(l.error());
    }
    return ParsedNumber(*l);
}

}  // namespace omnivar::util
