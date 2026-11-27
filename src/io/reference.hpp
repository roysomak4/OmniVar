#pragma once

// ReferenceResource.getReference: region bases plus the seed
// index of every 17-mer (SEED_1) and 12-mer (SEED_2) start position.

#include <expected>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/config.hpp"
#include "io/reference_provider.hpp"
#include "region/region.hpp"

namespace senri::io {

struct Reference {
    struct LoadedRegion {
        std::string chr;
        int sequence_start;
        int sequence_end;
        bool operator==(const LoadedRegion&) const = default;
    };
    std::vector<LoadedRegion> loaded_regions;
    std::unordered_map<int, char> reference_sequences;  // position -> upper-cased base
    // seed sequence -> start positions, in insertion order
    std::unordered_map<std::string, std::vector<int>> seed;
};

[[nodiscard]] bool is_loaded(const std::string& chr, int sequence_start, int sequence_end, const Reference& reference);

// ReferenceResource.getReference(region): default extension (-Y, referenceExtension).
[[nodiscard]] std::expected<Reference, ReferenceError> get_reference(const ReferenceProvider& provider,
                                                                     const Region& region, const Config& config,
                                                                     const ChrLengths& chr_lengths);

// ReferenceResource.getReference(region, extension, ref): adds the region (+ -x + extension) to `ref`.
[[nodiscard]] std::expected<void, ReferenceError> get_reference(const ReferenceProvider& provider, const Region& region,
                                                                int extension, Reference& ref, const Config& config,
                                                                const ChrLengths& chr_lengths);

// Opens the -G reference: a thin-FASTA CSV when the path ends in ".csv", else an indexed FASTA.
[[nodiscard]] std::expected<std::unique_ptr<ReferenceProvider>, ReferenceError> open_reference(const std::string& path);

}  // namespace senri::io
