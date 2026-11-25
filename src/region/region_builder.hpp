#pragma once

#include <expected>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cli/config.hpp"
#include "region/region.hpp"
#include "util/regex_engine.hpp"

namespace omnivar {

// Result of reading the BED file.
struct BedFile {
    std::optional<std::string> amplicon_parameters;  // set if -a, or if the BED looks like an amplicon BED
    std::optional<bool> zero_based;
    std::vector<std::string> lines;  // without #, browser and track lines
};

// Reads the BED from disk. Lines end at \n, \r or \r\n.
[[nodiscard]] std::expected<BedFile, std::string> read_bed_file(const Config& config);
// Same, from text (unit tests).
[[nodiscard]] std::expected<BedFile, std::string> read_bed_text(const Config& config, std::string_view text);

class RegionBuilder {
public:
    // Fails only if the -d delimiter is not a supported regular expression.
    [[nodiscard]] static std::expected<RegionBuilder, std::string> create(const ChrLengths& chr_lengths,
                                                                          const Config& config);

    [[nodiscard]] std::expected<Segments, std::string> build_regions(const std::vector<std::string>& seg_raws,
                                                                     std::optional<bool> zero_based) const;

    [[nodiscard]] std::expected<Segments, std::string> build_amp_regions(const std::vector<std::string>& seg_raws,
                                                                         bool zero_based) const;

    [[nodiscard]] std::expected<Segments, std::string> build_region_from_configuration() const;

    [[nodiscard]] std::string correct_chromosome(const std::string& chr) const;

    // Test hook replacing correctChromosome, like Mockito.doReturn in RegionBuilderTest.
    void set_chromosome_override(std::string chr) { chromosome_override_ = std::move(chr); }

private:
    RegionBuilder(const ChrLengths& chr_lengths, const Config& config, util::Regex delimiter)
        : chr_lengths_(&chr_lengths), config_(&config), delimiter_(std::move(delimiter)) {}

    const ChrLengths* chr_lengths_;
    const Config* config_;
    util::Regex delimiter_;
    std::optional<std::string> chromosome_override_;
};

}  // namespace omnivar
