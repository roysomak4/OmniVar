#pragma once

// Command line: the upstream CmdParser on top of an emulation of commons-cli 1.2 BasicParser,
// so every upstream command line parses to the same Configuration (checked against the jar by
// omnivar-testkit cli-parity).
//
// Options added by this port are long-form only and are taken out of the argument list before
// the upstream rules apply: --version, --help, --tsv, --vcf-top-allele-only, --filter-*, --no-filter-pstd,
// --dev-dump-config, --dev-dump-inputs, --dev-tsv-to-vcf.

#include <expected>
#include <span>
#include <string>
#include <string_view>

#include "cli/config.hpp"

namespace omnivar::cli {

enum class Action {
    help,         // no options, -H or -?; also --help
    version,      // --version
    run,          // run the caller
    dump_config,  // --dev-dump-config: print the parsed Configuration (parity testing)
    dump_inputs,  // --dev-dump-inputs: print samples, regions, reference and read counts (parity testing)
    tsv_to_vcf,   // --dev-tsv-to-vcf FILE: VCF from a Simple TSV (no caller run)
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

// Configuration in the test oracle's format (tests/oracle; sorted "name=value" lines).
[[nodiscard]] std::string dump_config(const Config& config);

}  // namespace omnivar::cli
