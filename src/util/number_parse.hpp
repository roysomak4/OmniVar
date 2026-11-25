#pragma once

// Number parsing with the upstream and commons-cli 1.2 rules.

#include <cstdint>
#include <expected>
#include <string>
#include <string_view>
#include <variant>

namespace omnivar::util {

// Integer.parseInt: optional sign, decimal digits, must fit in int. No whitespace.
[[nodiscard]] std::expected<std::int32_t, std::string> parse_int(std::string_view s);

// Long.parseLong / Long.valueOf.
[[nodiscard]] std::expected<std::int64_t, std::string> parse_long(std::string_view s);

// Integer.decode: optional sign, then "0x"/"0X"/"#" (hex), "0" followed by digits (octal) or decimal.
[[nodiscard]] std::expected<std::int32_t, std::string> decode_int(std::string_view s);

// Double.parseDouble / Double.valueOf: leading/trailing chars <= ' ' trimmed, optional sign,
// decimal or hex floating literal with optional exponent and optional f/F/d/D suffix,
// "NaN", "Infinity".
[[nodiscard]] std::expected<double, std::string> parse_double(std::string_view s);

// A Number as produced by commons-cli 1.2 TypeHandler.createNumber: a Double when the
// text contains '.', otherwise a Long.
class ParsedNumber {
public:
    explicit ParsedNumber(std::int64_t v) : v_(v) {}
    explicit ParsedNumber(double v) : v_(v) {}

    // Number.intValue(): Long narrows by truncating to 32 bits; Double casts as an (int) cast does.
    [[nodiscard]] std::int32_t int_value() const;
    // Number.doubleValue().
    [[nodiscard]] double double_value() const;
    [[nodiscard]] bool is_long() const { return std::holds_alternative<std::int64_t>(v_); }

private:
    std::variant<std::int64_t, double> v_;
};

[[nodiscard]] std::expected<ParsedNumber, std::string> create_number(std::string_view s);

// (int) cast of a double: NaN -> 0, saturates at the int range, truncates toward zero.
[[nodiscard]] std::int32_t d2i_cast(double d);

// (long) cast of a double, same rules.
[[nodiscard]] std::int64_t d2l_cast(double d);

}  // namespace omnivar::util
