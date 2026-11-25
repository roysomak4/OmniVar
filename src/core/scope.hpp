#pragma once

// What every pipeline stage sees: the immutable run configuration and context
// (GlobalReadOnlyScope), the region and its reference, plus the mutable per-region splice set
// (Scope.splice) and the run-wide exception counter (Configuration.exceptionCounter).

#include <atomic>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "cli/config.hpp"
#include "io/reference.hpp"
#include "region/region.hpp"

namespace omnivar::io {
class BamCache;
class ReferenceProvider;
}  // namespace omnivar::io

namespace omnivar::core {

// GlobalReadOnlyScope, minus conf.
struct RunContext {
    ChrLengths chr_lengths;
    std::vector<std::pair<std::string, long long>> contigs;  // header order (VCF output)
    std::string sample;
    std::string samplem;
    std::optional<std::string> amplicon_based_calling;
    std::unordered_map<std::string, int> adaptor_forward;
    std::unordered_map<std::string, int> adaptor_reverse;
    Segments segments;
};

// Thrown when more than Configuration.MAX_EXCEPTION_COUNT exceptions were continued past; the
// run then stops with exit status 1.
struct FatalError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Configuration.exceptionCounter (an AtomicInteger shared by all worker threads).
struct ExceptionCounter {
    std::atomic<int> count{0};
};

void print_exception_and_continue(const std::exception& e, const std::string& place, const std::string& place_def,
                                  const Region* region, ExceptionCounter& counter);

struct RegionScope {
    const Config& conf;
    const RunContext& run;
    Region region;
    // Scope.regionRef: the realigner extends it (ReferenceResource.getReference) before
    // re-reading reads outside the region.
    io::Reference& reference;
    std::unordered_set<std::string>& splice;
    ExceptionCounter& exceptions;
    // What VariationRealigner needs to re-read data: the reference source, the open BAMs and
    // Scope.bam (this sample's BAM paths joined with ':'). Null/empty in unit tests.
    const io::ReferenceProvider* provider = nullptr;
    io::BamCache* bam_cache = nullptr;
    std::string bam;
};

// RecordPreprocessor.getChrName: strips "chr" with -C.
[[nodiscard]] inline std::string chr_name(const Config& conf, const Region& region) {
    if (conf.chromosome_name_is_number && region.chr.starts_with("chr")) {
        return region.chr.substr(3);
    }
    return region.chr;
}

}  // namespace omnivar::core
