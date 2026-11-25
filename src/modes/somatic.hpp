#pragma once

// Somatic (tumour/normal) mode:

#include <string>
#include <unordered_set>

#include "modes/pipeline.hpp"

namespace omnivar::modes {

// SomaticMode.printHeader (empty unless -h)
[[nodiscard]] std::string somatic_header(const Config& conf);

// SomaticMode.processBothBamsInPipeline for one region: both pipelines (tumour = -b bam1,
// normal = bam2) and SomaticPostProcessModule.accept; appends the output lines to `out`.
void process_somatic_region(const PipelineEnv& env, const Region& region, std::string& out);

}  // namespace omnivar::modes
