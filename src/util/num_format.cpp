#include "util/num_format.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <limits>

#include "util/checked_string.hpp"
#include "util/number_parse.hpp"

namespace omnivar::util {

namespace {

constexpr std::string_view kInfinitySymbol = "\xE2\x88\x9E";  // U+221E, DecimalFormatSymbols.getInfinity()

// Shortest decimal digits of a finite positive double and its decimal exponent, such that
// value ~= 0.d1d2d3... x 10^exp10. Implements the JDK 19+ Double.toString digit selection.
struct Decimal {
    std::string digits;  // no leading or trailing zeros, at least one digit
    int exp10;           // value = 0.digits * 10^exp10
};

Decimal parse_scientific(std::string_view sci) {
    // "d.ddde[+-]xx" or "de[+-]xx"
    Decimal d;
    auto e = sci.find('e');
    std::string_view mant = sci.substr(0, e);
    int exp = 0;
    std::from_chars(sci.data() + e + 1 + (sci[e + 1] == '+' ? 1 : 0), sci.data() + sci.size(), exp);
    for (char c : mant) {
        if (c != '.') {
            d.digits.push_back(c);
        }
    }
    while (d.digits.size() > 1 && d.digits.back() == '0') {
        d.digits.pop_back();
    }
    d.exp10 = exp + 1;
    return d;
}

Decimal shortest_decimal(double v) {
    std::array<char, 64> buf{};
    auto [p, ec] = std::to_chars(buf.data(), buf.data() + buf.size(), v, std::chars_format::scientific);
    (void)ec;
    Decimal d = parse_scientific(std::string_view(buf.data(), p));
    if (d.digits.size() == 1) {
        // JDK 19+: when the shortest length is 1, choose the closest decimal of length 1 or 2.
        auto [p2, ec2] = std::to_chars(buf.data(), buf.data() + buf.size(), v, std::chars_format::scientific, 1);
        (void)ec2;
        d = parse_scientific(std::string_view(buf.data(), p2));
    }
    return d;
}

enum class RoundMode { half_up, half_even };

// Sign of (exact value of a) - (0.digits x 10^exp10), for finite positive a.
int compare_exact_to_digits(double a, const Decimal& d) {
    // Exact decimal expansion: a double has at most 767 significant decimal digits.
    std::array<char, 800> buf{};
    auto [p, ec] = std::to_chars(buf.data(), buf.data() + buf.size(), a, std::chars_format::scientific, 770);
    (void)ec;
    std::string_view sci(buf.data(), p);
    const auto e = sci.find('e');
    int exp = 0;
    std::from_chars(sci.data() + e + 1 + (sci[e + 1] == '+' ? 1 : 0), sci.data() + sci.size(), exp);
    if (exp + 1 != d.exp10) {
        return exp + 1 > d.exp10 ? 1 : -1;
    }
    std::string exact;
    for (char c : sci.substr(0, e)) {
        if (c != '.') {
            exact.push_back(c);
        }
    }
    for (std::size_t i = 0; i < exact.size(); ++i) {
        const char want = i < d.digits.size() ? d.digits[i] : '0';
        if (exact[i] != want) {
            return exact[i] > want ? 1 : -1;
        }
    }
    return d.digits.size() > exact.size() ? -1 : 0;
}

// Fixed notation of a finite positive double with `precision` fraction digits, rounding the
// shortest decimal digits (not the exact binary value) at the requested position.
std::string digits_fixed(double a, int precision, RoundMode mode) {
    if (a == 0.0) {
        std::string z = "0";
        if (precision > 0) {
            z += '.';
            z.append(static_cast<std::size_t>(precision), '0');
        }
        return z;
    }
    Decimal d = shortest_decimal(a);
    const long keep = static_cast<long>(d.exp10) + precision;  // digits kept
    std::string digits = d.digits;
    int exp10 = d.exp10;
    if (keep < 0) {
        digits = "0";
        exp10 = 1;
    } else if (static_cast<std::size_t>(keep) < digits.size()) {
        const auto k = static_cast<std::size_t>(keep);
        const char next = digits[k];
        bool round_up;
        if (mode == RoundMode::half_up) {
            round_up = next >= '5';
        } else if (next < '5') {
            round_up = false;
        } else if (next > '5' || k + 1 < digits.size()) {
            round_up = true;  // above 5, or 5 followed by more (non-zero) digits
        } else {
            // Known divergence (2 of ~8900 JVM vectors, none in the caller's ranges): 5.0E-4 and
            // 5.0E-5 round down on the JVM although the double lies above the tie.
            const int cmp = compare_exact_to_digits(a, d);
            if (cmp != 0) {
                round_up = cmp > 0;
            } else {
                round_up = k > 0 && ((digits[k - 1] - '0') % 2 != 0);
            }
        }
        digits.resize(k);
        if (round_up) {
            int i = static_cast<int>(digits.size()) - 1;
            while (i >= 0 && digits[static_cast<std::size_t>(i)] == '9') {
                digits[static_cast<std::size_t>(i)] = '0';
                --i;
            }
            if (i < 0) {
                digits.insert(digits.begin(), '1');
                ++exp10;
            } else {
                ++digits[static_cast<std::size_t>(i)];
            }
        }
        if (digits.empty()) {
            digits = "0";
            exp10 = 1;
        }
    }
    std::string int_part;
    std::string frac_part;
    if (exp10 <= 0) {
        int_part = "0";
        frac_part = std::string(static_cast<std::size_t>(-exp10), '0') + digits;
    } else if (static_cast<std::size_t>(exp10) >= digits.size()) {
        int_part = digits + std::string(static_cast<std::size_t>(exp10) - digits.size(), '0');
    } else {
        int_part = digits.substr(0, static_cast<std::size_t>(exp10));
        frac_part = digits.substr(static_cast<std::size_t>(exp10));
    }
    while (int_part.size() > 1 && int_part[0] == '0') {
        int_part.erase(0, 1);
    }
    if (frac_part.size() < static_cast<std::size_t>(precision)) {
        frac_part.append(static_cast<std::size_t>(precision) - frac_part.size(), '0');
    } else {
        frac_part.resize(static_cast<std::size_t>(precision));
    }
    std::string out = int_part;
    if (precision > 0) {
        out += '.';
        out += frac_part;
    }
    return out;
}

}  // namespace

int double_compare(double a, double b) {
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    // equal, or at least one NaN: compare raw bits as Double.compare does
    const bool an = std::isnan(a);
    const bool bn = std::isnan(b);
    if (an || bn) {
        return an == bn ? 0 : (an ? 1 : -1);
    }
    const bool as = std::signbit(a);
    const bool bs = std::signbit(b);
    return as == bs ? 0 : (as ? -1 : 1);
}

long long math_round(double value) {
    if (std::isnan(value)) {
        return 0;
    }
    double r = std::floor(value);
    if (value - r >= 0.5) {  // value - floor(value) is exact for doubles
        r += 1.0;
    }
    return d2l_cast(r);
}

std::string decimal_format(double value, int fraction_digits) {
    if (std::isnan(value)) {
        return "NaN";
    }
    const bool negative = value < 0.0 || (value == 0.0 && std::signbit(value));
    const double a = std::fabs(value);
    std::string out = negative ? "-" : "";
    if (std::isinf(a)) {
        out += kInfinitySymbol;
        return out;
    }
    out += digits_fixed(a, fraction_digits, RoundMode::half_even);
    return out;
}

double round_half_even(double value, int fraction_digits) {
    std::string s = decimal_format(value, fraction_digits);
    auto d = parse_double(s);
    if (!d) {
        throw RecoverableError("NumberFormatException: For input string: \"" + s + "\"");
    }
    return *d;
}

std::string rounded_value_to_print(double value, int fraction_digits) {
    if (value == static_cast<double>(math_round(value))) {
        return decimal_format(value, 0);
    }
    std::string s = decimal_format(value, fraction_digits);
    // .replaceAll("0+$", ""): may leave a trailing '.', e.g. 0.99999 with "0.0" -> "1." (the JVM does too)
    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }
    return s;
}

std::string printf_fixed(double value, int precision) {
    if (std::isnan(value)) {
        return "NaN";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Inf" : "-Inf";
    }
    char buf[400];
    auto r = std::to_chars(buf, buf + sizeof buf, value, std::chars_format::fixed, precision);
    return std::string(buf, r.ptr);
}

std::string string_format_fixed(double value, int precision) {
    if (std::isnan(value)) {
        return "NaN";
    }
    const bool negative = std::signbit(value);  // Double.compare(value, 0.0) == -1, includes -0.0
    const double a = std::fabs(value);
    if (std::isinf(a)) {
        return negative ? "-Infinity" : "Infinity";
    }
    return (negative ? "-" : "") + digits_fixed(a, precision, RoundMode::half_up);
}

std::string double_to_string(double value) {
    if (std::isnan(value)) {
        return "NaN";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Infinity" : "-Infinity";
    }
    std::string out = std::signbit(value) ? "-" : "";
    const double a = std::fabs(value);
    if (a == 0.0) {
        return out + "0.0";
    }
    Decimal d = shortest_decimal(a);
    const std::string& g = d.digits;
    if (a >= 1e-3 && a < 1e7) {
        if (d.exp10 <= 0) {
            out += "0.";
            out.append(static_cast<std::size_t>(-d.exp10), '0');
            out += g;
        } else if (static_cast<std::size_t>(d.exp10) >= g.size()) {
            out += g;
            out.append(static_cast<std::size_t>(d.exp10) - g.size(), '0');
            out += ".0";
        } else {
            out += g.substr(0, static_cast<std::size_t>(d.exp10));
            out += '.';
            out += g.substr(static_cast<std::size_t>(d.exp10));
        }
        return out;
    }
    out += g[0];
    out += '.';
    out += g.size() > 1 ? g.substr(1) : "0";
    out += 'E';
    out += std::to_string(d.exp10 - 1);
    return out;
}

}  // namespace omnivar::util
