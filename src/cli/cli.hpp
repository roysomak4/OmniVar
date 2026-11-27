#pragma once

// Command line: the upstream option rules (short options, defaults and their interactions) on top
// of a small getopt-style parser that behaves like the upstream one.
//
// The long-form options (--version, --help, --tsv, --vcf-top-allele-only, --filter-*,
// --no-filter-pstd, --dev-tsv-to-vcf) are taken out of the argument list before those rules apply.

#include <expected>
#include <span>
#include <string>
#include <string_view>

#include "cli/config.hpp"

namespace senri::cli {

enum class Action {
    help,        // no options, -H or -?; also --help
    version,     // --version
    run,         // run the caller
    tsv_to_vcf,  // --dev-tsv-to-vcf FILE: VCF from a Simple TSV (no caller run)
};

struct ParseError {
    std::string message;
};

struct ParsedCommandLine {
    Action action = Action::run;
    Config config;
    // Text written to stderr while parsing (e.g. the missing -G notice, missing -b).
    std::string warnings;
    std::string tsv_path;  // Action::tsv_to_vcf ("-" = stdin)
};

// args excludes argv[0]. `available_processors` is used for -th without a value.
[[nodiscard]] std::expected<ParsedCommandLine, ParseError> parse_command_line(std::span<const std::string_view> args,
                                                                              int available_processors);

[[nodiscard]] std::string usage();
[[nodiscard]] std::string version_string();

}  // namespace senri::cli
