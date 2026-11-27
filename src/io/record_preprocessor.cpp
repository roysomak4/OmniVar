#include "io/record_preprocessor.hpp"

#include <chrono>

namespace senri::io {

std::expected<const BamFile*, BamError> BamCache::get(const std::string& path) {
    auto it = files_.find(path);
    if (it == files_.end()) {
        auto f = BamFile::open(path);
        if (!f) {
            return std::unexpected(f.error());
        }
        it = files_.emplace(path, std::make_unique<BamFile>(std::move(*f))).first;
    }
    return it->second.get();
}

std::expected<RecordPreprocessor, BamError> RecordPreprocessor::create(const std::vector<std::string>& bams,
                                                                       const Region& region, const Config& config,
                                                                       BamCache& cache) {
    RecordPreprocessor p(region, config, cache);
    p.rnd_.seed(static_cast<std::uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()));
    // ArrayDeque.push (addFirst) of each BAM, then pollLast: BAMs are read in the given order.
    for (const auto& b : bams) {
        p.bams_.push_front(b);
    }
    if (auto r = p.next_reader(); !r) {
        return std::unexpected(r.error());
    }
    return p;
}

std::expected<bool, BamError> RecordPreprocessor::next_record(BamRecord& out) {
    auto r = next_on_current_reader(out);
    if (!r || *r) {
        return r;
    }
    if (auto n = next_reader(); !n) {
        return std::unexpected(n.error());
    }
    return next_on_current_reader(out);
}

std::expected<void, BamError> RecordPreprocessor::next_reader() {
    if (bams_.empty()) {
        return {};
    }
    const std::string path = bams_.back();
    bams_.pop_back();
    auto file = cache_->get(path);
    if (!file) {
        return std::unexpected(file.error());
    }
    auto view = SamView::query(**file, config_->samfilter, region_);
    if (!view) {
        return std::unexpected(view.error());
    }
    current_file_ = *file;
    current_reader_.emplace(std::move(*view));
    duplicates_.clear();
    first_matching_position_ = -1;
    return {};
}

std::expected<bool, BamError> RecordPreprocessor::next_on_current_reader(BamRecord& out) {
    do {
        if (!current_reader_->read(out)) {
            if (const auto& e = current_reader_->error()) {
                return std::unexpected(BamError{*e + " in " + current_file_->path()});
            }
            return false;
        }
    } while (!preprocess_record(out));
    return true;
}

bool RecordPreprocessor::preprocess_record(const BamRecord& record) {
    const Config& conf = *config_;
    if (conf.is_downsampling()) {
        // Random.nextDouble() <= downsampling
        const double r = std::uniform_real_distribution<double>(0.0, 1.0)(rnd_);
        if (r <= *conf.downsampling) {
            return false;
        }
    }
    // Ignore low mapping quality reads
    if (conf.has_mapping_quality() && record.mapping_quality() < *conf.mapping_quality) {
        return false;
    }
    // Skip not primary alignment reads
    if (record.secondary() && conf.samfilter != "0") {
        return false;
    }
    // Skip reads where sequence is not stored in read ("*")
    if (record.read_length() == 0) {
        return false;
    }
    ++total_reads;

    if (conf.remove_duplicated_reads) {
        const std::string_view mate_reference_name = current_file_->mate_reference_name_or_equals(record);
        if (record.alignment_start() != first_matching_position_) {
            duplicates_.clear();
        }
        if (record.mate_alignment_start() < 10) {
            // POS-RNEXT-PNEXT
            std::string key = std::to_string(record.alignment_start()) + "-" + std::string(mate_reference_name) + "-" +
                              std::to_string(record.mate_alignment_start());
            if (duplicates_.contains(key)) {
                ++duplicate_reads;
                return false;
            }
            duplicates_.insert(std::move(key));
            first_matching_position_ = record.alignment_start();
        } else if (record.read_paired() && record.mate_unmapped()) {
            // POS-CIGAR
            std::string key = std::to_string(record.alignment_start()) + "-" + record.cigar_string();
            if (duplicates_.contains(key)) {
                ++duplicate_reads;
                return false;
            }
            duplicates_.insert(std::move(key));
            first_matching_position_ = record.alignment_start();
        }
    }
    return true;
}

}  // namespace senri::io
