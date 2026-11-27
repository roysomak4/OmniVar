#pragma once

// Counts totalReads and duplicateReads.

#include <deque>
#include <expected>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "cli/config.hpp"
#include "io/bam.hpp"
#include "region/region.hpp"

namespace senri::io {

// Opens each BAM once and hands out the open file (one reader per file per
// thread in SamView.threadLocalSAMReaders).
class BamCache {
public:
    [[nodiscard]] std::expected<const BamFile*, BamError> get(const std::string& path);

private:
    std::map<std::string, std::unique_ptr<BamFile>> files_;
};

class RecordPreprocessor {
public:
    // `bams`: one sample's BAM paths (Configuration.BamNames bam1 or bam2 split on ':').
    [[nodiscard]] static std::expected<RecordPreprocessor, BamError> create(const std::vector<std::string>& bams,
                                                                            const Region& region, const Config& config,
                                                                            BamCache& cache);

    // RecordPreprocessor.nextRecord: the next record that passes the filters; false at the end.
    // The file of the returned record is available from current_file().
    [[nodiscard]] std::expected<bool, BamError> next_record(BamRecord& out);

    [[nodiscard]] const BamFile& current_file() const { return *current_file_; }

    int total_reads = 0;
    int duplicate_reads = 0;

private:
    RecordPreprocessor(const Region& region, const Config& config, BamCache& cache)
        : region_(region), config_(&config), cache_(&cache) {}

    std::expected<void, BamError> next_reader();
    std::expected<bool, BamError> next_on_current_reader(BamRecord& out);
    bool preprocess_record(const BamRecord& record);

    Region region_;
    const Config* config_;
    BamCache* cache_;
    std::deque<std::string> bams_;
    const BamFile* current_file_ = nullptr;
    std::optional<SamView> current_reader_;
    std::unordered_set<std::string> duplicates_;
    int first_matching_position_ = -1;
    std::mt19937_64 rnd_;
};

}  // namespace senri::io
