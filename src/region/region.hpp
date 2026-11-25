#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace omnivar {

struct Region {
    std::string chr;
    int start = 0;
    int end = 0;
    std::string gene;
    int insert_start = 0;  // amplicon mode: where to look for variants
    int insert_end = 0;

    bool operator==(const Region&) const = default;

    // Region.printRegion: chr:start-end
    [[nodiscard]] std::string print_region() const {
        return chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    }

    // Region.newModifiedRegion
    [[nodiscard]] Region with_bounds(int new_start, int new_end) const {
        Region r = *this;
        r.start = new_start;
        r.end = new_end;
        return r;
    }
};

// Segments as built by RegionBuilder: one inner list per BED line (or per amplicon group).
using Segments = std::vector<std::vector<Region>>;

// Contig name -> length, from the BAM header.
using ChrLengths = std::unordered_map<std::string, int>;

}  // namespace omnivar
