#include "io/reference.hpp"

#include <algorithm>
#include <cctype>

#include "io/csv_slice_provider.hpp"
#include "io/faidx_provider.hpp"

namespace omnivar::io {

namespace {

// Utils.substr(string, begin, len) for begin >= 0, len > 0 (ReferenceResource's only use).
std::string_view substr(std::string_view s, std::size_t begin, std::size_t len) {
    return s.substr(begin, std::min(len, s.size() - begin));
}

// ReferenceResource.addPositionsToSeedSequence
void add_seed(Reference& ref, int sequence_start, int i, std::string_view key) {
    ref.seed[std::string(key)].push_back(i + sequence_start);
}

}  // namespace

bool is_loaded(const std::string& chr, int sequence_start, int sequence_end, const Reference& reference) {
    for (const auto& r : reference.loaded_regions) {
        if (chr == r.chr && sequence_start >= r.sequence_start && sequence_end <= r.sequence_end) {
            return true;
        }
    }
    return false;
}

std::expected<Reference, ReferenceError> get_reference(const ReferenceProvider& provider, const Region& region,
                                                       const Config& config, const ChrLengths& chr_lengths) {
    Reference ref;
    if (auto r = get_reference(provider, region, config.reference_extension, ref, config, chr_lengths); !r) {
        return std::unexpected(r.error());
    }
    return ref;
}

std::expected<void, ReferenceError> get_reference(const ReferenceProvider& provider, const Region& region,
                                                  int extension, Reference& ref, const Config& config,
                                                  const ChrLengths& chr_lengths) {
    const int ext = config.number_nucleotide_to_extend;
    const int sequence_start = region.start - ext - extension < 1 ? 1 : region.start - ext - extension;
    const auto it = chr_lengths.find(region.chr);
    const int len = it != chr_lengths.end() ? it->second : 0;
    const int sequence_end = region.end + ext + extension > len ? len : region.end + ext + extension;

    auto sub = provider.sub_seq(region.chr, sequence_start, sequence_end);
    if (!sub) {
        return std::unexpected(sub.error());
    }
    std::string exon = std::move(*sub);
    // String.toUpperCase (bases are ASCII)
    std::transform(exon.begin(), exon.end(), exon.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (is_loaded(region.chr, sequence_start, sequence_end, ref)) {
        return {};
    }
    ref.loaded_regions.push_back(Reference::LoadedRegion{region.chr, sequence_start, sequence_end});

    const int exon_len = static_cast<int>(exon.size());
    // To process ends of chromosomes without decreasing by SEED_1
    const int site_end = len == sequence_end ? exon_len : exon_len - Config::kSeed1;
    for (int i = 0; i < site_end; ++i) {
        // don't process it more than once
        if (ref.reference_sequences.contains(i + sequence_start)) {
            continue;
        }
        ref.reference_sequences.emplace(i + sequence_start, exon[static_cast<std::size_t>(i)]);

        // Do not create adaptor sequences for the very end of chromosome
        if (len == sequence_end && i > exon_len - Config::kSeed1) {
            continue;
        }
        const auto ui = static_cast<std::size_t>(i);
        add_seed(ref, sequence_start, i, substr(exon, ui, Config::kSeed1));
        add_seed(ref, sequence_start, i, substr(exon, ui, Config::kSeed2));
    }
    return {};
}

std::expected<std::unique_ptr<ReferenceProvider>, ReferenceError> open_reference(const std::string& path) {
    if (is_csv_reference(path)) {
        auto p = CsvSliceProvider::open(path);
        if (!p) {
            return std::unexpected(p.error());
        }
        return std::make_unique<CsvSliceProvider>(std::move(*p));
    }
    auto p = FaidxProvider::open(path);
    if (!p) {
        return std::unexpected(p.error());
    }
    return std::make_unique<FaidxProvider>(std::move(*p));
}

}  // namespace omnivar::io
