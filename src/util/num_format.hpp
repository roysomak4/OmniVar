#pragma once

// Number formatting that matches the upstream output (JVM rules). Every number that reaches output goes through here.
//
// - decimal_format(v, n):     new DecimalFormat("0.<n zeros>").format(v)  (HALF_EVEN, exact value)
// - round_half_even(v, n):    Utils.roundHalfEven("0.<n zeros>", v)
// - rounded_value_to_print:   Utils.getRoundedValueToPrint("0.<n zeros>", v)
// - string_format_fixed(v,n): String.format("%.<n>f", v)  (HALF_UP on Double.toString digits)
// - double_to_string(v):      Double.toString(v) / string concatenation of a double (JDK 19+)
//
// Verified against the JVM: tests/golden/num_format_vectors.tsv.

#include <string>
#include <string_view>

namespace omnivar::util {

// DecimalFormat with pattern "0" (fraction_digits == 0) or "0." followed by fraction_digits zeros.
// RoundingMode.HALF_EVEN on the shortest decimal digits of the double; a '5' at the rounding
// position is resolved by where the exact binary value lies (above -> up, below -> down,
// exact -> even), as JDK 8+ DigitList does. Negative values that round to zero keep the sign
// ("-0.0"). NaN -> "NaN", infinities -> "∞"/"-∞" (en_US DecimalFormatSymbols).
[[nodiscard]] std::string decimal_format(double value, int fraction_digits);

// Double.parseDouble(new DecimalFormat(pattern).format(value)). Throws util::RecoverableError for
// infinite values ("\u221E" does not parse).
[[nodiscard]] double round_half_even(double value, int fraction_digits);

[[nodiscard]] std::string rounded_value_to_print(double value, int fraction_digits);

// String.format("%.<precision>f", value) in the ROOT/English locale.
[[nodiscard]] std::string string_format_fixed(double value, int precision);

// Perl/C sprintf("%.<precision>f", value): the exact binary value correctly rounded (ties to
// even), as glibc does; Perl's own stringification of non-finite values ("Inf", "-Inf", "NaN").
// For numbers the Perl annotation scripts format (var2vcf_paired.pl), not caller output.
[[nodiscard]] std::string printf_fixed(double value, int precision);

// Double.toString(value), JDK 19+ algorithm (shortest decimal that rounds back; two digits
// minimum). Plain notation for 1e-3 <= |v| < 1e7, otherwise computerized scientific ("1.0E7").
[[nodiscard]] std::string double_to_string(double value);

// Double.compare(a, b): NaN is greater than everything (and equal to itself), -0.0 < 0.0.
[[nodiscard]] int double_compare(double a, double b);

// Math.round(double): floor(v + 0.5) with the JVM's special cases (NaN -> 0, saturation).
[[nodiscard]] long long math_round(double value);

}  // namespace omnivar::util
