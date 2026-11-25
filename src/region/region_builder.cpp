#include "region/region_builder.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

#include "util/hash_order.hpp"
#include "util/number_parse.hpp"

namespace omnivar {

namespace {

// BufferedReader.readLine: lines end at \n, \r or \r\n; no empty line after a final terminator.
std::vector<std::string_view> read_lines(std::string_view text) {
    std::vector<std::string_view> lines;
    std::size_t i = 0;
    while (i < text.size()) {
        std::size_t j = i;
        while (j < text.size() && text[j] != '\n' && text[j] != '\r') {
            ++j;
        }
        lines.push_back(text.substr(i, j - i));
        if (j < text.size() && text[j] == '\r' && j + 1 < text.size() && text[j + 1] == '\n') {
            ++j;
        }
        i = j + 1;
    }
    return lines;
}

// Utils.toInt = Integer.parseInt, with the upstream exception text.
std::expected<int, std::string> to_int(std::string_view s) {
    auto r = util::parse_int(s);
    if (!r) {
        return std::unexpected("NumberFormatException: " + r.error());
    }
    return *r;
}

std::string out_of_bounds(std::size_t index, std::size_t length) {
    return "ArrayIndexOutOfBoundsException: Index " + std::to_string(index) + " out of bounds for length " +
           std::to_string(length);
}

template <class V>
std::expected<std::string_view, std::string> column(const V& cols, int index) {
    if (index < 0 || static_cast<std::size_t>(index) >= cols.size()) {
        return std::unexpected(out_of_bounds(static_cast<std::size_t>(index), cols.size()));
    }
    return cols[static_cast<std::size_t>(index)];
}

const util::Regex& comma_regex() {
    static const util::Regex re = util::Regex::compile(",").value();
    return re;
}

}  // namespace

std::expected<BedFile, std::string> read_bed_file(const Config& config) {
    std::ifstream in(config.bed.value(), std::ios::binary);
    if (!in) {
        return std::unexpected("FileNotFoundException: " + config.bed.value() + " (No such file or directory)");
    }
    std::ostringstream buf;
    buf << in.rdbuf();
    return read_bed_text(config, buf.str());
}

std::expected<BedFile, std::string> read_bed_text(const Config& config, std::string_view text) {
    auto delim = util::Regex::compile(config.delimiter);
    if (!delim) {
        return std::unexpected("unsupported -d delimiter: " + delim.error());
    }
    static const util::Regex integer_only = util::Regex::compile("^\\d+$").value();  // Patterns.INTEGER_ONLY

    BedFile bed;
    bed.amplicon_parameters = config.amplicon_based_calling;
    bed.zero_based = config.zero_based;
    for (std::string_view line : read_lines(text)) {
        if (line.starts_with("#") || line.starts_with("browser") || line.starts_with("track")) {
            continue;
        }
        if (!bed.amplicon_parameters) {
            auto cols = util::regex_split(line, *delim);
            // Amplicon mode only if the BED has exactly 8 columns and columns 7 and 8 are numbers.
            if (cols.size() == 8 && integer_only.find(cols[6]) && integer_only.find(cols[7])) {
                auto a1 = to_int(cols[1]);
                auto a2 = to_int(cols[2]);
                auto a6 = to_int(cols[6]);
                auto a7 = to_int(cols[7]);
                if (!a1 || !a2 || !a6 || !a7) {
                    return std::unexpected(
                        "Incorrect format of BED file for amplicon mode. It must be 8 columns and 2, 3, 7 and 8 "
                        "columns must contain region and amplicon starts and ends.\n" +
                        (!a1   ? a1.error()
                         : !a2 ? a2.error()
                         : !a6 ? a6.error()
                               : a7.error()));
                }
                if (*a6 >= *a1 && *a7 <= *a2) {
                    bed.amplicon_parameters = Config::kDefaultAmpliconParameters;
                    if (!config.is_zero_based_defined()) {
                        bed.zero_based = true;
                    }
                }
            }
        }
        bed.lines.emplace_back(line);
    }
    return bed;
}

std::expected<RegionBuilder, std::string> RegionBuilder::create(const ChrLengths& chr_lengths, const Config& config) {
    auto delim = util::Regex::compile(config.delimiter);
    if (!delim) {
        return std::unexpected("unsupported -d delimiter: " + delim.error());
    }
    return RegionBuilder(chr_lengths, config, std::move(*delim));
}

std::string RegionBuilder::correct_chromosome(const std::string& chr) const {
    if (chromosome_override_) {
        return *chromosome_override_;
    }
    if (!chr_lengths_->contains(chr)) {
        if (chr.starts_with("chr")) {
            return chr.substr(3);
        }
        return "chr" + chr;
    }
    return chr;
}

std::expected<Segments, std::string> RegionBuilder::build_regions(const std::vector<std::string>& seg_raws,
                                                                  std::optional<bool> zero_based) const {
    const Config& config = *config_;
    bool is_zero_based = config.is_zero_based_defined() ? *config.zero_based : false;
    BedRowFormat format = config.bed_row_format;
    Segments segs;
    for (const std::string& seg : seg_raws) {
        auto cols = util::regex_split(seg, delimiter_);
        if (!config.is_column_for_chromosome_set() && cols.size() == 4) {
            auto a1 = to_int(cols[1]);
            auto a2 = to_int(cols[2]);
            if (!a1 || !a2) {
                return std::unexpected(
                    "Incorrect format of BED file. 2 and 3 columns must contain region start and end.\n" +
                    (!a1 ? a1.error() : a2.error()));
            }
            if (*a1 <= *a2) {
                format = kCustomBedRowFormat;
                if (!zero_based) {
                    is_zero_based = true;
                }
            }
        }
        auto chr_col = column(cols, format.chr_column);
        if (!chr_col) {
            return std::unexpected(chr_col.error());
        }
        const std::string chr = correct_chromosome(std::string(*chr_col));
        auto start_col = column(cols, format.start_column);
        if (!start_col) {
            return std::unexpected(start_col.error());
        }
        auto cds_start = to_int(*start_col);
        if (!cds_start) {
            return std::unexpected(cds_start.error());
        }
        auto end_col = column(cols, format.end_column);
        if (!end_col) {
            return std::unexpected(end_col.error());
        }
        auto cds_end = to_int(*end_col);
        if (!cds_end) {
            return std::unexpected(cds_end.error());
        }
        const std::string gene = (format.gene_column >= 0 && static_cast<std::size_t>(format.gene_column) < cols.size())
                                     ? std::string(cols[static_cast<std::size_t>(format.gene_column)])
                                     : chr;

        auto ts_col = column(cols, format.thick_start_column);
        if (!ts_col) {
            return std::unexpected(ts_col.error());
        }
        auto te_col = column(cols, format.thick_end_column);
        if (!te_col) {
            return std::unexpected(te_col.error());
        }
        const auto thick_starts = util::regex_split(*ts_col, comma_regex());
        const auto thick_ends = util::regex_split(*te_col, comma_regex());
        std::vector<Region> thick_regions;
        for (std::size_t i = 0; i < thick_starts.size(); ++i) {
            auto thick_start = to_int(thick_starts[i]);
            if (!thick_start) {
                return std::unexpected(thick_start.error());
            }
            if (i >= thick_ends.size()) {
                return std::unexpected(out_of_bounds(i, thick_ends.size()));
            }
            auto thick_end = to_int(thick_ends[i]);
            if (!thick_end) {
                return std::unexpected(thick_end.error());
            }
            int ts = *thick_start;
            int te = *thick_end;
            if (*cds_start > te) {
                continue;
            }
            if (*cds_end > te) {
                break;
            }
            if (ts < *cds_start) {
                ts = *cds_start;
            }
            if (te > *cds_end) {
                te = *cds_end;
            }
            ts -= config.number_nucleotide_to_extend;
            te += config.number_nucleotide_to_extend;
            if (is_zero_based && ts < te) {
                ++ts;
            }
            thick_regions.push_back(Region{chr, ts, te, gene, 0, 0});
        }
        segs.push_back(std::move(thick_regions));
    }
    return segs;
}

std::expected<Segments, std::string> RegionBuilder::build_amp_regions(const std::vector<std::string>& seg_raws,
                                                                      bool zero_based) const {
    // tsegs: HashMap<String, List<Region>>; keep first-insertion order of chromosomes to
    // reproduce the HashMap iteration order below.
    std::vector<std::string> chr_order;
    std::map<std::string, std::vector<Region>> tsegs;
    for (const std::string& line : seg_raws) {
        auto split = util::regex_split(line, delimiter_);
        auto get = [&](int idx) { return column(split, idx); };
        auto chr_col = get(kAmpBedRowFormat.chr_column);
        auto start_col = get(kAmpBedRowFormat.start_column);
        auto end_col = get(kAmpBedRowFormat.end_column);
        auto gene_col = get(kAmpBedRowFormat.gene_column);
        auto is_col = get(kAmpBedRowFormat.thick_start_column);
        auto ie_col = get(kAmpBedRowFormat.thick_end_column);
        for (const auto* c : {&chr_col, &start_col, &end_col, &gene_col, &is_col, &ie_col}) {
            if (!*c) {
                return std::unexpected(c->error());
            }
        }
        const std::string chr = correct_chromosome(std::string(*chr_col));
        auto start = to_int(*start_col);
        auto end = to_int(*end_col);
        auto insert_start = to_int(*is_col);
        auto insert_end = to_int(*ie_col);
        for (const auto* v : {&start, &end, &insert_start, &insert_end}) {
            if (!*v) {
                return std::unexpected(v->error());
            }
        }
        int s = *start;
        int is = *insert_start;
        if (zero_based && s < *end) {
            ++s;
            ++is;
        }
        auto [it, inserted] = tsegs.try_emplace(chr);
        if (inserted) {
            chr_order.push_back(chr);
        }
        it->second.push_back(Region{chr, s, *end, std::string(*gene_col), is, *insert_end});
    }

    std::vector<std::int32_t> hashes;
    hashes.reserve(chr_order.size());
    for (const auto& c : chr_order) {
        hashes.push_back(util::string_hash_code(c));
    }
    std::vector<std::string> iteration;
    if (auto order = util::hashmap_order(hashes)) {
        for (std::size_t idx : *order) {
            iteration.push_back(chr_order[idx]);
        }
    } else {
        iteration = chr_order;  // tree bins (> 48 colliding contigs): order not modelled
    }

    Segments segs;
    segs.emplace_back();
    int previous_end = -1;
    for (const std::string& chr : iteration) {
        auto& chr_regions = tsegs[chr];
        // List.sort is a stable merge sort
        std::stable_sort(chr_regions.begin(), chr_regions.end(),
                         [](const Region& a, const Region& b) { return a.insert_start < b.insert_start; });
        const std::string* previous_chr = nullptr;
        for (const Region& region : chr_regions) {
            if (previous_end != -1 &&
                (previous_chr == nullptr || region.chr != *previous_chr || region.insert_start > previous_end)) {
                segs.emplace_back();
            }
            segs.back().push_back(region);
            previous_chr = &region.chr;
            previous_end = region.insert_end;
        }
    }
    return segs;
}

std::expected<Segments, std::string> RegionBuilder::build_region_from_configuration() const {
    const Config& config = *config_;
    static const util::Regex colon = util::Regex::compile(":").value();
    static const util::Regex dash = util::Regex::compile("-").value();
    const std::string& roi = config.region_of_interest.value();
    auto split = util::regex_split(roi, colon);
    if (split.empty()) {
        return std::unexpected(out_of_bounds(0, 0));
    }
    const std::string chr = correct_chromosome(std::string(split[0]));
    const std::string gene = split.size() < 3 ? chr : std::string(split[2]);
    if (split.size() < 2) {
        return std::unexpected(out_of_bounds(1, split.size()));
    }
    auto range = util::regex_split(split[1], dash);
    if (range.empty()) {
        return std::unexpected(out_of_bounds(0, 0));
    }
    auto strip_commas = [](std::string_view s) {
        std::string out;
        for (char c : s) {
            if (c != ',') {
                out += c;
            }
        }
        return out;
    };
    auto start = to_int(strip_commas(range[0]));
    if (!start) {
        return std::unexpected(start.error());
    }
    int s = *start;
    int e = s;
    if (range.size() >= 2) {
        auto end = to_int(strip_commas(range[1]));
        if (!end) {
            return std::unexpected(end.error());
        }
        e = *end;
    }
    s -= config.number_nucleotide_to_extend;
    e += config.number_nucleotide_to_extend;
    const bool zero_based = config.is_zero_based_defined() ? *config.zero_based : false;
    if (zero_based && s < e) {
        ++s;
    }
    if (s > e) {
        s = e;
    }
    return Segments{{Region{chr, s, e, gene, 0, 0}}};
}

}  // namespace omnivar
