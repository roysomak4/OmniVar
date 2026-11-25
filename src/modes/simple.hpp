#pragma once

// Simple (single-sample) mode output:

#include <string>

#include "core/scope.hpp"
#include "core/variant.hpp"

namespace omnivar::modes {

// SimpleOutputVariant.toString for one variant (nullptr = a null variant).
[[nodiscard]] std::string simple_output_line(const core::Variant* variant, const Region& region, const std::string& sv,
                                             int position, const Config& conf, const core::RunContext& run);

// -p (rejected variants): a rejected variant's output line is prefixed with this byte, its
// FILTER ID (core::reject_id) and a tab; the launcher strips the prefix and writes the reason
// (TSV: a last column, VCF: FILTER). Lines without it are ordinary calls.
inline constexpr char kRejectedLineMarker = '\x1f';

// SimpleMode.printHeader (empty unless -h)
[[nodiscard]] std::string simple_header(const Config& conf);

// SimplePostProcessModule.accept: appends output lines for one region to `out`.
void simple_post_process(core::AlignedVars& aligned, const core::RegionScope& scope, std::string& out);

}  // namespace omnivar::modes
