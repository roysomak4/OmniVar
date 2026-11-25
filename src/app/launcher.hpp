#pragma once

#include <expected>
#include <optional>
#include <string>
#include <unordered_map>

#include "cli/config.hpp"
#include "core/scope.hpp"
#include "region/region.hpp"

namespace omnivar::app {

using RunContext = core::RunContext;

[[nodiscard]] std::expected<RunContext, std::string> init_resources(const Config& config);

[[nodiscard]] std::expected<std::pair<std::string, std::string>, std::string> sample_names(const Config& config);

// --dev-dump-inputs: the M2 parity dump.
int dump_inputs(const Config& config);

// --dev-tsv-to-vcf: VCF from a Simple-mode TSV (36 columns, optional CRISPR column),
// as teststrandbias.R | var2vcf_valid.pl would write it. `path` "-" reads stdin.
int tsv_to_vcf(const Config& config, const std::string& path);

// Normal run. Returns the process exit status.
int run(const Config& config);

}  // namespace omnivar::app
