#pragma once

// SAMFileParser -> CigarParser -> VariationRealigner -> StructuralVariantsProcessor ->
// ToVarsBuilder for one region and one sample (Scope.bam, BAM paths joined with ':').

#include <string>
#include <unordered_set>

#include "core/scope.hpp"
#include "core/variant.hpp"
#include "io/record_preprocessor.hpp"
#include "io/reference_provider.hpp"

namespace senri::modes {

// What every pipeline run needs besides its region, BAMs and reference.
struct PipelineEnv {
    const Config& conf;
    const core::RunContext& run;
    const io::ReferenceProvider& provider;
    io::BamCache& cache;
    core::ExceptionCounter& exceptions;
};

// Scope<AlignedVarsData>: the variants and Scope.maxReadLength (CigarParser's final value).
struct PipelineResult {
    core::AlignedVars aligned;
    int max_read_length = 0;
};

// One pipeline run. `max_read_length` is Scope.maxReadLength on entry (0 for a first sample).
// Errors that stop the run are core::FatalError.
[[nodiscard]] PipelineResult run_pipeline(const PipelineEnv& env, const Region& region, const std::string& bam,
                                          int max_read_length, io::Reference& reference,
                                          std::unordered_set<std::string>& splice);

// "a.bam:b.bam" -> {"a.bam", "b.bam"} (String.split(":")).
[[nodiscard]] std::vector<std::string> split_bams(const std::string& bam);

}  // namespace senri::modes
