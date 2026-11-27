#pragma once

// Helpers shared by the var2vcf ports: fields of a TSV row read the way the Perl scripts read
// them (internal to annotate/).

#include <string>
#include <string_view>
#include <vector>

#include "util/num_format.hpp"
#include "util/number_parse.hpp"

namespace senri::annotate::detail {

inline std::vector<std::string> split_tabs(std::string_view s) {
    // Perl split(/\t/) drops trailing empty fields.
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const auto p = s.find('\t', start);
        out.emplace_back(s.substr(start, p == std::string_view::npos ? std::string_view::npos : p - start));
        if (p == std::string_view::npos) {
            break;
        }
        start = p + 1;
    }
    while (!out.empty() && out.back().empty()) {
        out.pop_back();
    }
    return out;
}

// Perl numeric value of a field. Fields are numbers as the caller prints them; anything that does not
// parse (an empty or missing field, an infinity "∞") is 0, as in Perl.
inline double num(const std::vector<std::string>& f, std::size_t i) {
    if (i >= f.size()) {
        return 0;
    }
    auto v = util::parse_double(f[i]);
    return v ? *v : 0;
}

inline const std::string& str(const std::vector<std::string>& f, std::size_t i) {
    static const std::string empty;
    return i < f.size() ? f[i] : empty;
}

// Perl truth of a string: "" and "0" are false.
inline bool perl_true(const std::string& s) { return !s.empty() && s != "0"; }

inline std::string integer_text(double v) { return util::decimal_format(v, 0); }

}  // namespace senri::annotate::detail
