#include "core/cigar_parser.hpp"

#include <cmath>
#include <cstdio>
#include <print>

#include "core/cigar_modifier.hpp"
#include "core/patterns.hpp"
#include "util/checked_string.hpp"
#include "util/number_parse.hpp"
#include "util/regex_engine.hpp"
#include "util/sequence.hpp"

namespace senri::core {

using util::char_at;
using util::qual_at;
using util::substr;

namespace {

// BEGIN_dig_dig_S_ANY_dig_dig_S_END "^\d\dS.*\d\dS$" on Cigar.toString(): the first element is
// a soft clip of exactly two digits (10-99) and the last a soft clip of >= 10, distinct elements.
bool both_ends_soft_clipped(const Cigar& c) {
    return c.size() >= 2 && c.front().op == 'S' && c.front().length >= 10 && c.front().length <= 99 &&
           c.back().op == 'S' && c.back().length >= 10;
}

// SamRecord.getBaseQualityString via CigarParser.getBaseQualityString:
// qualities above SAMUtils.MAX_PHRED_SCORE (93) make htsjdk throw, and the fallback encodes
// (char)(q + 33) of the signed byte without actually capping.
std::string cigar_base_quality_string(const io::BamRecord& record) {
    const bam1_t* b = record.raw();
    const int n = b->core.l_qseq;
    const uint8_t* q = bam_get_qual(b);
    if (n == 0 || q[0] == 0xff) {
        return "*";
    }
    std::string s(static_cast<std::size_t>(n), '!');
    bool over = false;
    for (int i = 0; i < n; ++i) {
        over |= q[i] > 93;
        s[static_cast<std::size_t>(i)] = static_cast<char>(q[i] + 33);
    }
    if (over) {
        std::println(stderr,
                     "WARNING: Cannot get encode qualities for SAM entry at {}:{}, record name: '{}'. Qualities "
                     "capped at 93",
                     "?", record.alignment_start(), record.read_name());
    }
    return s;
}

// Utils.getReverseComplementedSequence(record, startIndex, length):
// Arrays.copyOfRange pads with 0 bytes past the end; a negative start throws.
std::string reverse_complemented_sequence(const std::string& bases, int start_index, int length) {
    if (start_index < 0) {
        start_index = static_cast<int>(bases.size()) + start_index;
    }
    if (start_index < 0 || start_index > static_cast<int>(bases.size()) || length < 0) {
        throw util::RecoverableError("ArrayIndexOutOfBoundsException: copyOfRange " + std::to_string(start_index));
    }
    std::string range(static_cast<std::size_t>(length), '\0');
    for (int i = 0; i < length; ++i) {
        const int j = start_index + i;
        if (j < static_cast<int>(bases.size())) {
            range[static_cast<std::size_t>(i)] = bases[static_cast<std::size_t>(j)];
        }
    }
    // SequenceUtil.reverseComplement
    std::string out(range.rbegin(), range.rend());
    for (char& c : out) {
        c = util::complement(c);
    }
    return out;
}

std::pair<int, std::string> adj_ins_pos(int bi, std::string ins, const RefMap& ref) {
    int n = 1;
    const int len = static_cast<int>(ins.size());
    while (true) {
        const char c = char_at(ins, static_cast<long>(ins.size()) - n);
        if (!is_equals(ref_get(ref, bi), c)) {
            break;
        }
        ++n;
        if (n > len) {
            n = 1;
        }
        --bi;
    }
    if (n > 1) {
        ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
    }
    return {bi, ins};
}

const util::Regex& comma() {
    static const util::Regex re = util::Regex::compile(",").value();
    return re;
}

}  // namespace

CigarParser::CigarParser(const RegionScope& scope, bool svflag, int max_read_length)
    : scope_(scope), conf_(scope.conf), ref_(scope.reference.reference_sequences), svflag_(svflag) {
    data_.max_read_length = max_read_length;
}

const CigarElement& CigarParser::element(int i) const {
    if (i < 0 || static_cast<std::size_t>(i) >= cigar_.size()) {
        throw util::RecoverableError("IndexOutOfBoundsException: Index " + std::to_string(i) +
                                     " out of bounds for length " + std::to_string(cigar_.size()));
    }
    return cigar_[static_cast<std::size_t>(i)];
}

namespace {

// Patterns.BEGIN_NUM_S_OR_BEGIN_NUM_H "^(\\d+)S|^\\d+H" find, group 1: the leading soft clip.
std::optional<int> begin_num_s(std::string_view cigar) {
    std::size_t d = patterns::digit_run(cigar, 0);
    if (d > 0 && d < cigar.size() && cigar[d] == 'S') {
        return util::parse_int(cigar.substr(0, d)).value_or(0);
    }
    return std::nullopt;
}

// Patterns.END_NUM_S_OR_NUM_H "(\\d+)S$|H$" find, group 1: the trailing soft clip.
std::optional<int> end_num_s(std::string_view cigar) {
    if (cigar.size() < 2 || cigar.back() != 'S' || !patterns::is_digit(cigar[cigar.size() - 2])) {
        return std::nullopt;
    }
    std::size_t b = cigar.size() - 1;
    while (b > 0 && patterns::is_digit(cigar[b - 1])) {
        --b;
    }
    return util::parse_int(cigar.substr(b, cigar.size() - 1 - b)).value_or(0);
}

// Patterns.MC_Z_NUM_S_ANY_NUM_S "\\d+S\\S*\\d+S" find: a soft clip followed (later) by another.
bool mc_two_soft_clips(std::string_view mc) {
    std::size_t first = std::string_view::npos;
    for (std::size_t i = 1; i < mc.size(); ++i) {
        if (mc[i] == 'S' && patterns::is_digit(mc[i - 1])) {
            first = i;
            break;
        }
    }
    if (first == std::string_view::npos) {
        return false;
    }
    for (std::size_t i = first + 2; i < mc.size(); ++i) {
        if (mc[i] == 'S' && patterns::is_digit(mc[i - 1])) {
            return true;
        }
    }
    return false;
}

Sclip& last_sv(std::vector<Sclip>& v) { return v.back(); }

}  // namespace

void CigarParser::prepare_sv_structures_for_analysis(const io::BamRecord& record, const io::BamFile& file,
                                                     const std::string& query_quality, int number_of_mismatches,
                                                     bool read_direction, bool mate_direction, int start,
                                                     int total_length_including_soft_clipped) {
    SVStructures& sv = data_.sv_structures;
    const int mate_start = record.mate_alignment_start();
    const int mend = mate_start + total_length_including_soft_clipped;
    const std::string cigar_text = cigar_to_string(cigar_);
    int end = start;
    for (const CigarElement& e : cigar_) {  // sum(globalFind(ALIGNED_LENGTH_MND, cigar))
        if (e.op == 'M' || e.op == 'N' || e.op == 'D') {
            end += e.length;
        }
    }
    const int max_read_length = data_.max_read_length;

    int soft5 = 0;
    if (auto tt = begin_num_s(cigar_text)) {
        if (*tt != 0 && util::qual_at(query_quality, *tt - 1) > conf_.goodq) {
            soft5 = start;
        }
    }
    int soft3 = 0;
    if (auto tt = end_num_s(cigar_text)) {
        if (*tt != 0 && util::qual_at(query_quality, static_cast<long>(query_quality.size()) - *tt) > conf_.goodq) {
            soft3 = end;
        }
    }
    // Mate direction is true if it is forward, read direction is true if it is reverse.
    const int read_dir_num = read_direction ? -1 : 1;
    const int mate_dir_num = mate_direction ? 1 : -1;
    const int MIN_D = 75;  // Minimal distance from start to end of SV
    const double min_svc_dist = Config::kMinSvcDist * max_read_length;
    const int quals_len = static_cast<int>(query_quality.size());
    const auto qual_at_minmapbase = [&] {
        return static_cast<double>(util::qual_at(query_quality, Config::kMinMapBase));
    };
    const int mq = record.mapping_quality();
    const auto mc = record.string_attribute("MC");
    const auto mate_mq = record.integer_attribute("MQ");
    const bool mc_clipped = mc && mc_two_soft_clips(*mc);
    const bool mate_low_mq = mate_mq && *mate_mq < 15;
    const auto disc_if = [](std::vector<Sclip>& v, bool cond) {
        if (!v.empty() && cond) {
            last_sv(v).disc++;
        }
    };
    // New cluster when there is none yet or the read starts far from the previous one.
    const auto cluster = [&](std::vector<Sclip>& v, int last_end) -> Sclip& {
        if (v.empty() || start - last_end > min_svc_dist) {
            v.emplace_back();
        }
        return last_sv(v);
    };

    if (file.mate_reference_name_or_equals(record) == "=") {
        int mlen = record.inferred_insert_size();
        if (mc_clipped || mate_low_mq) {
            // Ignore those with mates mapped with softclipping at both ends, or with mate mapping
            // quality less than 15
        } else if (read_dir_num * mate_dir_num == -1 && (mlen * read_dir_num) > 0 && quals_len > Config::kMinMapBase) {
            // deletion candidate
            mlen = mate_start > start ? mend - start : end - mate_start;
            if (std::abs(mlen) > conf_.ins_size + conf_.ins_std_amt * conf_.ins_std) {
                if (read_dir_num == 1) {
                    Sclip& c = cluster(sv.svfdel, sv.svdelfend);
                    add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, mlen,
                           soft3, max_read_length / 2.0, qual_at_minmapbase(), mq, number_of_mismatches);
                    sv.svdelfend = end;
                } else {
                    Sclip& c = cluster(sv.svrdel, sv.svdelrend);
                    add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, mlen,
                           soft5, max_read_length / 2.0, qual_at_minmapbase(), mq, number_of_mismatches);
                    sv.svdelrend = end;
                }
                disc_if(sv.svfdel, std::abs(start - sv.svdelfend) <= min_svc_dist);
                disc_if(sv.svrdel, std::abs(start - sv.svdelrend) <= min_svc_dist);
                disc_if(sv.svfdup, std::abs(start - sv.svdupfend) <= MIN_D);
                disc_if(sv.svrdup, std::abs(start - sv.svduprend) <= MIN_D);
                disc_if(sv.svfinv5, std::abs(start - sv.svinvfend5) <= MIN_D);
                disc_if(sv.svrinv5, std::abs(start - sv.svinvrend5) <= MIN_D);
                disc_if(sv.svfinv3, std::abs(start - sv.svinvfend3) <= MIN_D);
                disc_if(sv.svrinv3, std::abs(start - sv.svinvrend3) <= MIN_D);
            }
        } else if (read_dir_num * mate_dir_num == -1 && read_dir_num * mlen < 0 && quals_len > Config::kMinMapBase) {
            // duplication
            if (read_dir_num == 1) {
                Sclip& c = cluster(sv.svfdup, sv.svdupfend);
                add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, mlen, soft3,
                       max_read_length / 2.0, qual_at_minmapbase(), mq, number_of_mismatches);
                sv.svdupfend = end;
            } else {
                Sclip& c = cluster(sv.svrdup, sv.svduprend);
                add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, mlen, soft5,
                       max_read_length / 2.0, qual_at_minmapbase(), mq, number_of_mismatches);
                sv.svduprend = end;
            }
            disc_if(sv.svfdup, std::abs(start - sv.svdupfend) <= min_svc_dist);
            disc_if(sv.svrdup, std::abs(start - sv.svduprend) <= min_svc_dist);
            disc_if(sv.svfdel, std::abs(start - sv.svdelfend) <= MIN_D);
            disc_if(sv.svrdel, std::abs(start - sv.svdelrend) <= MIN_D);
            disc_if(sv.svfinv5, std::abs(start - sv.svinvfend5) <= MIN_D);
            disc_if(sv.svrinv5, std::abs(start - sv.svinvrend5) <= MIN_D);
            disc_if(sv.svfinv3, std::abs(start - sv.svinvfend3) <= MIN_D);
            disc_if(sv.svrinv3, std::abs(start - sv.svinvrend3) <= MIN_D);
        } else if (read_dir_num * mate_dir_num == 1 && quals_len > Config::kMinMapBase) {  // Inversion
            const auto inv = [&](std::vector<Sclip>& v, int& last_end, int softp) {
                Sclip& c = cluster(v, last_end);
                add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, mlen, softp,
                       max_read_length / 2.0, qual_at_minmapbase(), mq, number_of_mismatches);
                last_end = end;
                last_sv(v).disc++;
            };
            if (read_dir_num == 1 && mlen != 0) {
                if (mlen < -3 * max_read_length) {
                    inv(sv.svfinv3, sv.svinvfend3, soft3);
                } else if (mlen > 3 * max_read_length) {
                    inv(sv.svfinv5, sv.svinvfend5, soft3);
                }
            } else if (mlen != 0) {
                if (mlen < -3 * max_read_length) {
                    inv(sv.svrinv3, sv.svinvrend3, soft5);
                } else if (mlen > 3 * max_read_length) {
                    inv(sv.svrinv5, sv.svinvrend5, soft5);
                }
            }
            if (mlen != 0) {
                disc_if(sv.svfdel, (start - sv.svdelfend) <= MIN_D);
                disc_if(sv.svrdel, (start - sv.svdelrend) <= MIN_D);
                disc_if(sv.svfdup, (start - sv.svdupfend) <= MIN_D);
                disc_if(sv.svrdup, (start - sv.svduprend) <= MIN_D);
            }
        }
    } else if (quals_len > Config::kMinMapBase) {  // Inter-chr translocation
        const std::string mchr(file.mate_reference_name_or_equals(record));
        if (mc_clipped || mate_low_mq) {
            // Ignore those with mates mapped with softclipping at both ends, or with mate mapping
            // quality less than 15
        } else {
            auto& fus = read_dir_num == 1 ? sv.svffus : sv.svrfus;
            auto& fus_end = read_dir_num == 1 ? sv.svfusfend : sv.svfusrend;
            std::vector<Sclip>* list = fus.get(mchr);
            bool new_cluster = list == nullptr;
            if (!new_cluster) {
                const int* last_end = fus_end.get(mchr);
                if (last_end == nullptr) {
                    util::throw_npe("svfusfend.get(mchr)");
                }
                new_cluster = start - *last_end > min_svc_dist;
            }
            if (new_cluster) {
                fus.get_or_create(mchr).emplace_back();
            }
            Sclip& c = fus.get(mchr)->back();
            add_sv(c, start, end, mate_start, mend, read_dir_num, total_length_including_soft_clipped, 0,
                   read_dir_num == 1 ? soft3 : soft5, max_read_length / 2.0, qual_at_minmapbase(), mq,
                   number_of_mismatches);
            fus_end.put(mchr, end);
            c.disc++;
        }
        disc_if(sv.svfdel, start - sv.svdelfend <= Config::kMinSvPos);
        disc_if(sv.svrdel, start - sv.svdelrend <= Config::kMinSvPos);
        disc_if(sv.svfdup, start - sv.svdupfend <= Config::kMinSvPos);
        disc_if(sv.svrdup, start - sv.svduprend <= Config::kMinSvPos);
        disc_if(sv.svfinv5, start - sv.svinvfend5 <= Config::kMinSvPos);
        disc_if(sv.svrinv5, start - sv.svinvrend5 <= Config::kMinSvPos);
        disc_if(sv.svfinv3, start - sv.svinvfend3 <= Config::kMinSvPos);
        disc_if(sv.svrinv3, start - sv.svinvrend3 <= Config::kMinSvPos);
    }
}

void CigarParser::add_sv(Sclip& sdref, int start_s, int end_e, int mate_start_ms, int mate_end_me, int dir, int rlen,
                         int mlen, int softp, double pmean_rp, double qmean, double big_qmean, double nm) const {
    sdref.vars_count++;
    sdref.inc_dir(dir != 1);
    if (qmean >= conf_.goodq) {
        sdref.high_quality_reads_count++;
    } else {
        sdref.low_quality_reads_count++;
    }
    if (sdref.start == 0 || sdref.start >= start_s) {
        sdref.start = start_s;
    }
    if (sdref.end == 0 || sdref.end <= end_e) {
        sdref.end = end_e;
    }
    sdref.mates.push_back(Mate{mate_start_ms, mate_end_me, mlen, start_s, end_e, pmean_rp, qmean, big_qmean, nm});
    if (sdref.mstart == 0 || sdref.mstart >= mate_start_ms) {
        sdref.mstart = mate_start_ms;
    }
    if (sdref.mend == 0 || sdref.mend <= mate_end_me) {
        sdref.mend = mate_start_ms + rlen;
    }
    if (softp != 0) {
        const int anchor = dir == 1 ? sdref.end : sdref.start;
        if (std::abs(softp - anchor) < 10) {
            const int* count = sdref.soft.get(softp);
            sdref.soft.put(softp, (count ? *count : 0) + 1);
        }
    }
}

CigarParserResult CigarParser::process(io::RecordPreprocessor& processor) {
    const int total_reads = 0;
    const int duplicate_reads = 0;

    const std::string chr = chr_name(conf_, scope_.region);
    io::BamRecord record;
    while (true) {
        auto r = processor.next_record(record);
        if (!r) {
            throw FatalError(r.error().message);
        }
        if (!*r) {
            break;
        }
        try {
            parse_cigar(chr, record, processor.current_file());
        } catch (const util::RecoverableError& e) {
            print_exception_and_continue(e, "record", std::string(record.read_name()), &scope_.region,
                                         scope_.exceptions);
        }
    }

    CigarParserResult result;
    if (conf_.output_splicing) {
        data_.splice_count.for_each([&](const std::string& intron, int count) {
            result.splicing_output +=
                scope_.run.sample + "\t" + scope_.region.chr + "\t" + intron + "\t" + std::to_string(count) + "\n";
        });
        result.data.max_read_length = data_.max_read_length;
        return result;
    }

    if (svflag_) {
        // Integer division, as upstream: the ratio is 0 or at least 1, so "> 0.5" means ">= 1".
        data_.duprate = (discordant_count_ + 1) / (total_reads - duplicate_reads + 1) >= 1 ? 0.0 : 1.0;
    } else {
        data_.duprate = 0.0;  // removeDuplicatedReads && totalReads != 0 is never true (see above)
    }
    result.data = std::move(data_);
    return result;
}

namespace {

void swap_initial_data(VariationData& a, VariationData& b) {
    std::swap(a.non_insertion_variants, b.non_insertion_variants);
    std::swap(a.insertion_variants, b.insertion_variants);
    std::swap(a.ref_coverage, b.ref_coverage);
    std::swap(a.soft_clips_3end, b.soft_clips_3end);
    std::swap(a.soft_clips_5end, b.soft_clips_5end);
}

}  // namespace

void CigarParser::process_into(io::RecordPreprocessor& processor, VariationData& shared) {
    swap_initial_data(data_, shared);
    CigarParserResult result;
    try {
        result = process(processor);
    } catch (...) {
        swap_initial_data(data_, shared);
        throw;
    }
    // Splicing output (-i) never reaches the realigner; there process() leaves data_ in place.
    swap_initial_data(conf_.output_splicing ? data_ : result.data, shared);
}

void CigarParser::parse_cigar(const std::string& chr_name, const io::BamRecord& record, const io::BamFile& file) {
    std::string query_sequence = record.read_string();
    const int mapping_quality = record.mapping_quality();
    cigar_ = cigar_of(record);
    const RefMap& ref = ref_;

    const int insertion_deletion_length = core::insertion_deletion_length(cigar_);

    // Number of mismatches from NM tag (or nM for STAR).
    int total_number_of_mismatches = 0;
    auto nm_tag = record.integer_attribute("NM");
    if (auto star = record.integer_attribute("nM")) {
        nm_tag = star;
    }
    if (nm_tag) {
        total_number_of_mismatches = static_cast<int>(*nm_tag) - insertion_deletion_length;
        if (total_number_of_mismatches > conf_.mismatch) {
            return;
        }
    } else {
        const std::string cigar_string = record.cigar_string();
        if (conf_.y && cigar_string != "*") {
            std::println(stderr, "No NM tag for mismatches. {}", record.read_name());
        }
        if (record.read_unmapped() || cigar_string == "*") {
            return;
        }
    }

    std::string query_quality = cigar_base_quality_string(record);
    const bool is_mate_reference_name_equal =
        file.reference_name(record.reference_index()) == file.reference_name(record.mate_reference_index());
    const int number_of_mismatches = total_number_of_mismatches;
    const bool direction = record.read_negative_strand();

    if (scope_.run.amplicon_based_calling) {
        if (parse_cigar_with_amp_case(record, is_mate_reference_name_equal)) {
            return;
        }
    }

    int position;
    read_position_including_soft_clipped_ = 0;
    read_position_excluding_soft_clipped_ = 0;

    if (conf_.perform_local_realignment) {
        ModifiedCigar mc = modify_cigar(record.alignment_start(), record.cigar_string(), query_sequence, query_quality,
                                        scope_.reference, insertion_deletion_length, scope_.region,
                                        data_.max_read_length, conf_, scope_.exceptions);
        position = mc.position;
        cigar_ = parse_cigar_text(mc.cigar);
        query_sequence = std::move(mc.query_sequence);
        query_quality = std::move(mc.query_quality);
    } else {
        position = record.alignment_start();
        cigar_ = cigar_of(record);
    }

    // cleanupCigar(record) replaces the record's own CIGAR, not the
    // `cigar` being parsed; it only matters through record.getCigar() below, whose reference
    // length it does not change. So it has no effect here.

    start_ = position;
    offset_ = 0;

    // determine discordant reads
    if (file.mate_reference_name_or_equals(record) != "=") {
        ++discordant_count_;
    }

    // Ignore reads that are softclipped at both ends and both greater than 10 bp
    if (both_ends_soft_clipped(cigar_)) {
        return;
    }
    // Only match and insertion counts toward read length
    const int read_length_include_matching_and_insertions = match_insertion_length(cigar_);
    if (conf_.minmatch != 0 && read_length_include_matching_and_insertions < conf_.minmatch) {
        return;
    }
    // The total length, including soft-clipped bases
    const int total_length_including_soft_clipped = soft_clipped_length(cigar_);
    if (total_length_including_soft_clipped > data_.max_read_length) {
        data_.max_read_length = total_length_including_soft_clipped;
    }

    // If supplementary alignment is present
    if (conf_.samfilter != "0" && record.supplementary()) {
        return;
    }

    // Skip sites that are not in region of interest in CRISPR mode
    if (skip_sites_out_region_of_interest()) {
        return;
    }

    // true if mate is in forward direction
    const bool mate_direction = (record.flags() & 0x20) == 0;

    if (record.read_paired() && record.mate_unmapped()) {
        // Mate unmapped, potential insertion: not implemented upstream either
    } else if (record.mapping_quality() > 10 && !conf_.disable_sv) {
        // Consider high mapping quality mates only
        prepare_sv_structures_for_analysis(record, file, query_quality, number_of_mismatches, direction, mate_direction,
                                           position, total_length_including_soft_clipped);
    }
    const int mate_alignment_start = record.mate_alignment_start();

    const int n_elements = static_cast<int>(cigar_.size());
    bool stop_processing = false;
    // Loop over CIGAR records
    for (int ci = 0; ci < n_elements && !stop_processing; ci++) {
        if (skip_overlapping_reads(record, file, position, direction, mate_alignment_start)) {
            break;
        }
        cigar_element_length_ = element(ci).length;

        // CigarParser.getCigarOperator: insertions at the edges are soft clips
        char op = element(ci).op;
        if ((ci == 0 || ci == n_elements - 1) && op == 'I') {
            op = 'S';
        }
        switch (op) {
            case 'N': process_not_matched(); continue;
            case 'S':
                process_soft_clip(chr_name, record, file, query_sequence, mapping_quality, query_quality,
                                  number_of_mismatches, direction, position, total_length_including_soft_clipped, ci);
                continue;
            case 'H': offset_ = 0; continue;
            case 'I':
                offset_ = 0;
                ci = process_insertion(query_sequence, mapping_quality, query_quality, number_of_mismatches, direction,
                                       position, read_length_include_matching_and_insertions, ci);
                continue;
            case 'D':
                offset_ = 0;
                ci = process_deletion(query_sequence, mapping_quality, query_quality, number_of_mismatches, direction,
                                      read_length_include_matching_and_insertions, ci);
                continue;
            default: break;
        }

        // Now dealing with matching part
        int nmoff = 0;
        int moffset = 0;
        for (int i = offset_; i < cigar_element_length_; i++) {
            const bool trim = is_trim_at_opt_t_bases(direction, total_length_including_soft_clipped);

            const char ch1 = char_at(query_sequence, read_position_including_soft_clipped_);
            std::string s(1, ch1);
            bool start_with_deletion = false;
            // skip if base is unknown
            if (ch1 == 'N') {
                if (conf_.include_n_in_total_depth) {
                    inc_cnt(data_.ref_coverage, start_, 1);
                }
                start_++;
                read_position_including_soft_clipped_++;
                read_position_excluding_soft_clipped_++;
                continue;
            }

            double q = qual_at(query_quality, read_position_including_soft_clipped_);
            int qbases = 1;
            int qibases = 0;
            std::string ss;

            while ((start_ + 1) >= scope_.region.start && (start_ + 1) <= scope_.region.end &&
                   (i + 1) < cigar_element_length_ && q >= conf_.goodq &&
                   is_has_and_not_equals(ref, start_, query_sequence, read_position_including_soft_clipped_) &&
                   is_not_equals('N', ref_get(ref, start_))) {
                // Require higher quality for MNV
                if (qual_at(query_quality, read_position_including_soft_clipped_ + 1) < conf_.goodq + 5) {
                    break;
                }
                const char nuc = char_at(query_sequence, read_position_including_soft_clipped_ + 1);
                if (nuc == 'N') {
                    break;
                }
                if (is_has_and_equals('N', ref, start_ + 1)) {
                    break;
                }
                if (is_not_equals(ref_get(ref, start_ + 1), nuc)) {
                    ss += nuc;
                    q += qual_at(query_quality, read_position_including_soft_clipped_ + 1);
                    qbases++;
                    read_position_including_soft_clipped_++;
                    read_position_excluding_soft_clipped_++;
                    i++;
                    start_++;
                    nmoff++;
                } else {
                    // look ahead to see whether there're more mismatches, if yes, add reference to grow MNV
                    int ssn = 0;
                    for (int ssi = 1; ssi <= conf_.vext; ssi++) {
                        if (i + 1 + ssi >= cigar_element_length_) {
                            break;
                        }
                        if (read_position_including_soft_clipped_ + 1 + ssi < static_cast<int>(query_sequence.size())) {
                            const char c = char_at(query_sequence, read_position_including_soft_clipped_ + 1 + ssi);
                            if (is_has_and_not_equals(c, ref, start_ + 1 + ssi)) {
                                ssn = ssi + 1;
                                break;
                            }
                        }
                    }
                    if (ssn == 0) {
                        break;
                    }
                    if (qual_at(query_quality, read_position_including_soft_clipped_ + ssn) < conf_.goodq + 5) {
                        break;
                    }
                    for (int ssi = 1; ssi <= ssn; ssi++) {
                        ss += char_at(query_sequence, read_position_including_soft_clipped_ + ssi);
                        q += qual_at(query_quality, read_position_including_soft_clipped_ + ssi);
                        qbases++;
                    }
                    read_position_including_soft_clipped_ += ssn;
                    read_position_excluding_soft_clipped_ += ssn;
                    i += ssn;
                    start_ += ssn;
                }
            }

            if (!ss.empty()) {
                s += "&" + ss;
            }
            int ddlen = 0;

            if (is_closer_then_vext_and_good_base(query_sequence, query_quality, ci, i, ss, 'D')) {
                while (i + 1 < cigar_element_length_) {
                    s += char_at(query_sequence, read_position_including_soft_clipped_ + 1);
                    q += qual_at(query_quality, read_position_including_soft_clipped_ + 1);
                    qbases++;
                    i++;
                    read_position_including_soft_clipped_++;
                    read_position_excluding_soft_clipped_++;
                    start_++;
                }
                s = util::replace_first(s, "&", "");
                s = "-" + std::to_string(element(ci + 1).length) + "&" + s;
                start_with_deletion = true;
                ddlen = element(ci + 1).length;
                ci += 1;

                if (two_insertions_ahead(ci)) {
                    s +=
                        "^" + substr(query_sequence, read_position_including_soft_clipped_ + 1, element(ci + 1).length);
                    const int next_len = element(ci + 1).length;
                    for (int qi = 1; qi <= next_len; qi++) {
                        q += qual_at(query_quality, read_position_including_soft_clipped_ + 1 + qi);
                        qibases++;
                    }
                    read_position_including_soft_clipped_ += next_len;
                    read_position_excluding_soft_clipped_ += next_len;
                    ci += 1;
                }
                if (next_after_num_matched(ci, 1)) {
                    Offset tpl = find_offset(start_ + ddlen + 1, read_position_including_soft_clipped_ + 1,
                                             element(ci + 1).length, query_sequence, query_quality);
                    if (tpl.offset != 0) {
                        moffset = tpl.offset;
                        nmoff += tpl.offset_number_of_mismatches;
                        s += "&" + tpl.sequence;
                        for (char tq : tpl.quality_sequence) {
                            q += static_cast<unsigned char>(tq) - 33;
                            qibases++;
                        }
                    }
                }
            } else if (is_closer_then_vext_and_good_base(query_sequence, query_quality, ci, i, ss, 'I')) {
                while (i + 1 < cigar_element_length_) {
                    s += char_at(query_sequence, read_position_including_soft_clipped_ + 1);
                    q += qual_at(query_quality, read_position_including_soft_clipped_ + 1);
                    qbases++;
                    i++;
                    read_position_including_soft_clipped_++;
                    read_position_excluding_soft_clipped_++;
                    start_++;
                }
                s = util::replace_first(s, "&", "");
                const int next_len = element(ci + 1).length;
                s += substr(query_sequence, read_position_including_soft_clipped_ + 1, next_len);
                s = substr(s, 0, next_len) + "&" + util::substring(s, next_len);
                s = "+" + s;
                for (int qi = 1; qi <= next_len; qi++) {
                    q += qual_at(query_quality, read_position_including_soft_clipped_ + 1 + qi);
                    qibases++;
                }
                read_position_including_soft_clipped_ += next_len;
                read_position_excluding_soft_clipped_ += next_len;
                ci += 1;
                qibases--;
                qbases++;  // need to add to set the correction insertion position
            }
            if (!trim) {
                const int pos = start_ - qbases + 1;
                if (pos >= scope_.region.start && pos <= scope_.region.end && !s.contains('N')) {
                    add_variation_for_matching_part(mapping_quality, number_of_mismatches, direction,
                                                    read_length_include_matching_and_insertions, nmoff, s,
                                                    start_with_deletion, q, qbases, qibases, ddlen, pos);
                }
            }

            if (start_with_deletion) {
                start_ += ddlen;
            }
            if (op != 'I') {
                start_++;
            }
            if (op != 'D') {
                read_position_including_soft_clipped_++;
                read_position_excluding_soft_clipped_++;
            }
            if (skip_overlapping_reads(record, file, position, direction, mate_alignment_start)) {
                stop_processing = true;  // break processCigar
                break;
            }
        }
        if (stop_processing) {
            break;
        }
        if (moffset != 0) {
            offset_ = moffset;
            read_position_including_soft_clipped_ += moffset;
            start_ += moffset;
            read_position_excluding_soft_clipped_ += moffset;
        }
        if (start_ > scope_.region.end) {
            break;
        }
    }
}

bool CigarParser::is_trim_at_opt_t_bases(bool direction, int total_length_including_soft_clipped) const {
    if (conf_.trim_bases_after != 0) {
        if (!direction) {
            return read_position_including_soft_clipped_ > conf_.trim_bases_after;
        }
        return total_length_including_soft_clipped - read_position_including_soft_clipped_ > conf_.trim_bases_after;
    }
    return false;
}

bool CigarParser::skip_sites_out_region_of_interest() const {
    const int cut_site = conf_.crispr_cutting_site;
    const int filter_bp = conf_.crispr_filtering_bp;
    if (cut_site != 0) {
        // ALIGNED_LENGTH_MD "(\d+)[MD=X]": total aligned length excluding clips and insertions
        int rlen3 = 0;
        for (const auto& e : cigar_) {
            if (e.op == 'M' || e.op == 'D' || e.op == '=' || e.op == 'X') {
                rlen3 += e.length;
            }
        }
        if (filter_bp != 0) {
            return !(cut_site - start_ > filter_bp && start_ + rlen3 - cut_site > filter_bp);
        }
    }
    return false;
}

int CigarParser::process_deletion(const std::string& query_sequence, int mapping_quality,
                                  const std::string& query_quality, int number_of_mismatches, bool direction,
                                  int read_length_include_matching_and_insertions, int ci) {
    const RefMap& ref = ref_;
    // Ignore deletions right after introns at exon edge in RNA-seq
    if (skip_indel_next_to_intron(ci)) {
        read_position_excluding_soft_clipped_ += cigar_element_length_;
        return ci;
    }

    std::string desc = "-" + std::to_string(cigar_element_length_);
    std::string sequence_to_append;
    const char quality_of_last_segment_before_del = char_at(query_quality, read_position_including_soft_clipped_ - 1);
    std::string quality_of_segment;

    int multoffs = 0;
    int multoffp = 0;
    int nmoff = 0;

    if (insertion_or_deletion_with_next_matched(ci)) {
        const int m_len = element(ci + 1).length;
        const int indel_len = element(ci + 2).length;
        const int begin = read_position_including_soft_clipped_;
        append_segments(query_sequence, query_quality, ci, desc, quality_of_segment, m_len, indel_len, begin, false);

        multoffs += m_len + (element(ci + 2).op == 'D' ? indel_len : 0);
        multoffp += m_len + (element(ci + 2).op == 'I' ? indel_len : 0);
        if (next_after_num_matched(ci, 3)) {
            int vsn = 0;
            const int tn = read_position_including_soft_clipped_ + multoffp;
            const int ts = start_ + multoffs + cigar_element_length_;
            for (int vi = 0; vsn <= conf_.vext && vi < element(ci + 3).length; vi++) {
                if (char_at(query_sequence, tn + vi) == 'N') {
                    break;
                }
                if (qual_at(query_quality, tn + vi) < conf_.goodq) {
                    break;
                }
                if (is_has_and_equals('N', ref, ts + vi)) {
                    break;
                }
                if (auto ref_ch = ref_get(ref, ts + vi)) {
                    if (is_not_equals(char_at(query_sequence, tn + vi), ref_ch)) {
                        offset_ = vi + 1;
                        nmoff++;
                        vsn = 0;
                    } else {
                        vsn++;
                    }
                }
            }
            if (offset_ != 0) {
                sequence_to_append += substr(query_sequence, tn, offset_);
                quality_of_segment += substr(query_quality, tn, offset_);
            }
        }
        // skip next 2 CIGAR segments
        ci += 2;
    } else if (next_insertion(ci)) {
        const int ins_len = element(ci + 1).length;
        desc += "^" + substr(query_sequence, read_position_including_soft_clipped_, ins_len);
        quality_of_segment += substr(query_quality, read_position_including_soft_clipped_, ins_len);
        multoffp += ins_len;
        if (next_after_num_matched(ci, 2)) {
            const int m_len = element(ci + 2).length;
            int vsn = 0;
            const int tn = read_position_including_soft_clipped_ + multoffp;
            const int ts = start_ + cigar_element_length_;
            for (int vi = 0; vsn <= conf_.vext && vi < m_len; vi++) {
                const char seq_ch = char_at(query_sequence, tn + vi);
                if (seq_ch == 'N') {
                    break;
                }
                if (qual_at(query_quality, tn + vi) < conf_.goodq) {
                    break;
                }
                if (auto ref_ch = ref_get(ref, ts + vi)) {
                    if (is_equals('N', ref_ch)) {
                        break;
                    }
                    if (is_not_equals(seq_ch, ref_ch)) {
                        offset_ = vi + 1;
                        nmoff++;
                        vsn = 0;
                    } else {
                        vsn++;
                    }
                }
            }
            if (offset_ != 0) {
                sequence_to_append += substr(query_sequence, tn, offset_);
                quality_of_segment += substr(query_quality, tn, offset_);
            }
        }
        ci += 1;
    } else {
        if (next_matched(ci)) {
            const int m_len = element(ci + 1).length;
            int vsn = 0;
            for (int vi = 0; vsn <= conf_.vext && vi < m_len; vi++) {
                const char seq_ch = char_at(query_sequence, read_position_including_soft_clipped_ + vi);
                if (seq_ch == 'N') {
                    break;
                }
                if (qual_at(query_quality, read_position_including_soft_clipped_ + vi) < conf_.goodq) {
                    break;
                }
                if (auto ref_ch = ref_get(ref, start_ + cigar_element_length_ + vi)) {
                    if (is_equals('N', ref_ch)) {
                        break;
                    }
                    if (is_not_equals(seq_ch, ref_ch)) {
                        offset_ = vi + 1;
                        nmoff++;
                        vsn = 0;
                    } else {
                        vsn++;
                    }
                }
            }
            if (offset_ != 0) {
                sequence_to_append += substr(query_sequence, read_position_including_soft_clipped_, offset_);
                quality_of_segment += substr(query_quality, read_position_including_soft_clipped_, offset_);
            }
        }
    }
    if (offset_ > 0) {
        desc += "&" + sequence_to_append;
    }

    // quality of first matched base after deletion: best of q1 and q2
    if (read_position_including_soft_clipped_ + offset_ >= static_cast<int>(query_quality.size())) {
        quality_of_segment += quality_of_last_segment_before_del;
    } else {
        const char q2 = char_at(query_quality, read_position_including_soft_clipped_ + offset_);
        quality_of_segment +=
            static_cast<unsigned char>(quality_of_last_segment_before_del) > static_cast<unsigned char>(q2)
                ? quality_of_last_segment_before_del
                : q2;
    }

    if (start_ >= scope_.region.start && start_ <= scope_.region.end) {
        add_variation_for_deletion(mapping_quality, number_of_mismatches, direction,
                                   read_length_include_matching_and_insertions, desc, quality_of_segment, nmoff);
    }

    start_ += cigar_element_length_ + offset_ + multoffs;
    read_position_including_soft_clipped_ += offset_ + multoffp;
    read_position_excluding_soft_clipped_ += offset_ + multoffp;
    return ci;
}

bool CigarParser::next_after_num_matched(int ci, int number) const {
    return static_cast<int>(cigar_.size()) > ci + number && element(ci + number).op == 'M';
}
bool CigarParser::two_insertions_ahead(int ci) const {
    return static_cast<int>(cigar_.size()) > ci + 1 && element(ci + 1).op == 'I';
}
bool CigarParser::next_insertion(int ci) const {
    return conf_.perform_local_realignment && static_cast<int>(cigar_.size()) > ci + 1 && element(ci + 1).op == 'I';
}
bool CigarParser::next_matched(int ci) const {
    return conf_.perform_local_realignment && static_cast<int>(cigar_.size()) > ci + 1 && element(ci + 1).op == 'M';
}
bool CigarParser::insertion_or_deletion_with_next_matched(int ci) const {
    return conf_.perform_local_realignment && static_cast<int>(cigar_.size()) > ci + 2 &&
           element(ci + 1).length <= conf_.vext && element(ci + 1).op == 'M' &&
           (element(ci + 2).op == 'I' || element(ci + 2).op == 'D') && element(ci + 3).op != 'I' &&
           element(ci + 3).op != 'D';
}

bool CigarParser::is_closer_then_vext_and_good_base(const std::string& query_sequence, const std::string& query_quality,
                                                    int ci, int i, const std::string& ss, char op) const {
    // Do not adjust complex if we have hard-clips after insertion/deletion
    if (static_cast<int>(cigar_.size()) > ci + 2 && element(ci + 2).op == 'H') {
        return false;
    }
    if (!(conf_.perform_local_realignment && cigar_element_length_ - i <= conf_.vext &&
          static_cast<int>(cigar_.size()) > ci + 1 && element(ci + 1).op == op && ref_.contains(start_))) {
        return false;
    }
    if (!ss.empty() ||
        is_not_equals(char_at(query_sequence, read_position_including_soft_clipped_), ref_get(ref_, start_))) {
        return qual_at(query_quality, read_position_including_soft_clipped_) >= conf_.goodq;
    }
    return false;
}

int CigarParser::process_insertion(const std::string& query_sequence, int mapping_quality,
                                   const std::string& query_quality, int number_of_mismatches, bool direction,
                                   int position, int read_length_include_matching_and_insertions, int ci) {
    const RefMap& ref = ref_;
    if (skip_indel_next_to_intron(ci)) {
        read_position_including_soft_clipped_ += cigar_element_length_;
        return ci;
    }

    std::string desc = substr(query_sequence, read_position_including_soft_clipped_, cigar_element_length_);
    std::string quality_string = substr(query_quality, read_position_including_soft_clipped_, cigar_element_length_);
    std::string ss;

    int multoffs = 0;
    int multoffp = 0;
    int nmoff = 0;

    if (insertion_or_deletion_with_next_matched(ci)) {
        const int m_len = element(ci + 1).length;
        const int indel_len = element(ci + 2).length;
        const int begin = read_position_including_soft_clipped_ + cigar_element_length_;
        append_segments(query_sequence, query_quality, ci, desc, quality_string, m_len, indel_len, begin, true);
        multoffs += m_len + (element(ci + 2).op == 'D' ? indel_len : 0);
        multoffp += m_len + (element(ci + 2).op == 'I' ? indel_len : 0);

        const int ci6 = static_cast<int>(cigar_.size()) > ci + 3 ? element(ci + 3).length : 0;
        if (ci6 != 0 && element(ci + 3).op == 'M') {
            Offset tpl =
                find_offset(start_ + multoffs, read_position_including_soft_clipped_ + cigar_element_length_ + multoffp,
                            ci6, query_sequence, query_quality);
            offset_ = tpl.offset;
            ss = tpl.sequence;
            quality_string += tpl.quality_sequence;
        }
        ci += 2;
    } else {
        if (next_matched(ci)) {
            int vsn = 0;
            for (int vi = 0; vsn <= conf_.vext && vi < element(ci + 1).length; vi++) {
                const int rp = read_position_including_soft_clipped_ + cigar_element_length_ + vi;
                if (char_at(query_sequence, rp) == 'N') {
                    break;
                }
                if (qual_at(query_quality, rp) < conf_.goodq) {
                    break;
                }
                if (ref_.contains(start_ + vi)) {
                    if (is_not_equals(char_at(query_sequence, rp), ref_get(ref, start_ + vi))) {
                        offset_ = vi + 1;
                        nmoff++;
                        vsn = 0;
                    } else {
                        vsn++;
                    }
                }
            }
            if (offset_ != 0) {
                ss += substr(query_sequence, read_position_including_soft_clipped_ + cigar_element_length_, offset_);
                quality_string +=
                    substr(query_quality, read_position_including_soft_clipped_ + cigar_element_length_, offset_);
                for (int osi = 0; osi < offset_; osi++) {
                    inc_cnt(data_.ref_coverage, start_ + osi, 1);
                }
            }
        }
    }

    if (offset_ > 0) {
        desc += "&" + ss;
    }

    if (start_ - 1 >= scope_.region.start && start_ - 1 <= scope_.region.end && !desc.contains('N')) {
        int insertion_position = start_ - 1;
        if (patterns::begin_atgc_end(desc)) {
            auto [adjusted_position, adjusted_desc] = adj_ins_pos(start_ - 1, desc, ref);
            // check for the index, if it is negative, there can be deletions, will not be adjusted.
            if (read_position_including_soft_clipped_ - 1 - (start_ - 1 - adjusted_position) > 0) {
                insertion_position = adjusted_position;
                desc = adjusted_desc;
            }
        }
        inc_cnt(data_.position_to_insertion_count.get_or_create(insertion_position), "+" + desc, 1);

        Variation& hv = get_variation(data_.insertion_variants, insertion_position, "+" + desc);
        hv.inc_dir(direction);
        hv.vars_count++;
        const int rpes = read_position_excluding_soft_clipped_;
        const int tp = rpes < read_length_include_matching_and_insertions - rpes
                           ? rpes + 1
                           : read_length_include_matching_and_insertions - rpes;

        double tmpq = 0;
        for (char c : quality_string) {
            tmpq += static_cast<unsigned char>(c) - 33;
        }
        tmpq = tmpq / static_cast<double>(quality_string.size());

        if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
            hv.pstd = true;
        }
        if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
            hv.qstd = true;
        }
        hv.mean_position += tp;
        hv.mean_quality += tmpq;
        hv.mean_mapping_quality += mapping_quality;
        hv.pp = tp;
        hv.pq = tmpq;
        if (tmpq >= conf_.goodq) {
            hv.high_quality_reads_count++;
        } else {
            hv.low_quality_reads_count++;
        }
        hv.number_of_mismatches += number_of_mismatches - nmoff;

        // Adjust the reference count for insertion reads
        const int index_in_query = read_position_including_soft_clipped_ - 1 - (start_ - 1 - insertion_position);
        if (insertion_position > position) {
            const char c = char_at(query_sequence, index_in_query);
            if (is_has_and_equals(c, ref, insertion_position)) {
                Variation* tv = get_variation_maybe(data_.non_insertion_variants, insertion_position,
                                                    char_at(query_sequence, index_in_query));
                if (tv != nullptr) {
                    sub_cnt(*tv, direction, tp, qual_at(query_quality, index_in_query), mapping_quality,
                            number_of_mismatches - nmoff);
                }
            }
        }
        // Adjust count if the insertion is at the edge so that the AF won't > 1
        if (ci == 1 && (element(0).op == 'S' || element(0).op == 'H')) {
            // Variation hv is not touched below, but get_variation may insert into another map;
            // take the values needed first (references into HashOrderMap values stay valid).
            const bool hv_pstd = hv.pstd;
            const bool hv_qstd = hv.qstd;
            Variation& ttref = get_variation(data_.non_insertion_variants, insertion_position,
                                             std::string(1, ref_at(ref, insertion_position)));
            ttref.inc_dir(direction);
            ttref.vars_count++;
            ttref.pstd = hv_pstd;
            ttref.qstd = hv_qstd;
            ttref.mean_position += tp;
            ttref.mean_quality += tmpq;
            ttref.mean_mapping_quality += mapping_quality;
            ttref.pp = tp;
            ttref.pq = tmpq;
            ttref.number_of_mismatches += number_of_mismatches - nmoff;
            inc_cnt(data_.ref_coverage, insertion_position, 1);
        }
    }

    read_position_including_soft_clipped_ += cigar_element_length_ + offset_ + multoffp;
    read_position_excluding_soft_clipped_ += cigar_element_length_ + offset_ + multoffp;
    start_ += offset_ + multoffs;
    return ci;
}

void CigarParser::process_soft_clip(const std::string& chr_name, const io::BamRecord& record, const io::BamFile& file,
                                    const std::string& query_sequence, int mapping_quality,
                                    const std::string& query_quality, int number_of_mismatches, bool direction,
                                    int position, int total_length_including_soft_clipped, int ci) {
    const RefMap& ref = ref_;
    if (ci == 0) {  // 5' soft clipped
        if (!conf_.chimeric) {
            auto sa = record.string_attribute("SA");
            if (cigar_element_length_ >= 20 && sa) {
                if (is_read_chimeric_with_sa(record, file, position, std::string(*sa), direction, true)) {
                    read_position_including_soft_clipped_ += cigar_element_length_;
                    offset_ = 0;
                    start_ = position;
                    return;
                }
            } else if (cigar_element_length_ >= Config::kSeed1) {
                const std::string sequence = reverse_complemented_sequence(query_sequence, 0, cigar_element_length_);
                const std::string seed = util::substring(sequence, 0, Config::kSeed1);
                auto it = scope_.reference.seed.find(seed);
                if (it != scope_.reference.seed.end()) {
                    const auto& positions = it->second;
                    if (positions.size() == 1 && std::abs(start_ - positions[0]) < 2 * data_.max_read_length) {
                        read_position_including_soft_clipped_ += cigar_element_length_;
                        offset_ = 0;
                        start_ = position;
                        if (conf_.y) {
                            std::println(stderr, "{} at 5' is a chimeric at {} by SEED {}", sequence, start_,
                                         Config::kSeed1);
                        }
                        return;
                    }
                }
            }
        }
        // Align softclipped but matched sequences due to mis-softclipping
        while (cigar_element_length_ - 1 >= 0 && start_ - 1 > 0) {
            auto len_it = scope_.run.chr_lengths.find(chr_name);
            if (len_it == scope_.run.chr_lengths.end()) {
                util::throw_npe("chrLengths.get(" + chr_name + ")");
            }
            if (!(start_ - 1 <= len_it->second)) {
                break;
            }
            const char c = char_at(query_sequence, cigar_element_length_ - 1);
            if (!is_has_and_equals(c, ref, start_ - 1)) {
                break;
            }
            if (!(qual_at(query_quality, cigar_element_length_ - 1) > 10)) {
                break;
            }
            Variation& variation =
                get_variation(data_.non_insertion_variants, start_ - 1, std::string(1, ref_at(ref, start_ - 1)));
            add_cnt(variation, direction, cigar_element_length_, qual_at(query_quality, cigar_element_length_ - 1),
                    mapping_quality, number_of_mismatches);
            inc_cnt(data_.ref_coverage, start_ - 1, 1);
            start_--;
            cigar_element_length_--;
        }
        if (cigar_element_length_ > 0) {
            int sum_of_read_qualities = 0;
            int number_of_high_quality_bases = 0;
            int number_of_low_quality_bases = 0;
            for (int si = cigar_element_length_ - 1; si >= 0; si--) {
                if (char_at(query_sequence, si) == 'N') {
                    break;
                }
                const int base_quality = qual_at(query_quality, si);
                if (base_quality <= 12) {
                    number_of_low_quality_bases++;
                }
                if (number_of_low_quality_bases > 1) {
                    break;
                }
                sum_of_read_qualities += base_quality;
                number_of_high_quality_bases++;
            }
            sclip5_high_quality_processing(query_sequence, mapping_quality, query_quality, number_of_mismatches,
                                           direction, sum_of_read_qualities, number_of_high_quality_bases,
                                           number_of_low_quality_bases);
        }
        cigar_element_length_ = element(ci).length;
    } else if (ci == static_cast<int>(cigar_.size()) - 1) {  // 3' soft clipped
        if (!conf_.chimeric) {
            auto sa = record.string_attribute("SA");
            if (cigar_element_length_ >= 20 && sa) {
                if (is_read_chimeric_with_sa(record, file, position, std::string(*sa), direction, false)) {
                    read_position_including_soft_clipped_ += cigar_element_length_;
                    offset_ = 0;
                    start_ = position;
                    return;
                }
            } else if (cigar_element_length_ >= Config::kSeed1) {
                const std::string sequence =
                    reverse_complemented_sequence(query_sequence, -cigar_element_length_, cigar_element_length_);
                const std::string seed = substr(sequence, -Config::kSeed1, Config::kSeed1);
                auto it = scope_.reference.seed.find(seed);
                if (it != scope_.reference.seed.end()) {
                    const auto& positions = it->second;
                    if (positions.size() == 1 && std::abs(start_ - positions[0]) < 2 * data_.max_read_length) {
                        read_position_including_soft_clipped_ += cigar_element_length_;
                        offset_ = 0;
                        start_ = position;
                        if (conf_.y) {
                            std::println(stderr, "{} at 3' is a chimeric at {} by SEED {}", sequence, start_,
                                         Config::kSeed1);
                        }
                        return;
                    }
                }
            }
        }
        while (read_position_including_soft_clipped_ < static_cast<int>(query_sequence.size())) {
            const char c = char_at(query_sequence, read_position_including_soft_clipped_);
            if (!is_has_and_equals(c, ref, start_)) {
                break;
            }
            if (!(qual_at(query_quality, read_position_including_soft_clipped_) > 10)) {
                break;
            }
            Variation& variation =
                get_variation(data_.non_insertion_variants, start_, std::string(1, ref_at(ref, start_)));
            add_cnt(variation, direction, total_length_including_soft_clipped - read_position_excluding_soft_clipped_,
                    qual_at(query_quality, read_position_including_soft_clipped_), mapping_quality,
                    number_of_mismatches);
            inc_cnt(data_.ref_coverage, start_, 1);
            read_position_including_soft_clipped_++;
            start_++;
            cigar_element_length_--;
            read_position_excluding_soft_clipped_++;
        }
        if (static_cast<int>(query_sequence.size()) - read_position_including_soft_clipped_ > 0) {
            int sum_of_read_qualities = 0;
            int number_of_high_quality_bases = 0;
            int number_of_low_quality_bases = 0;
            for (int si = 0; si < cigar_element_length_; si++) {
                if (char_at(query_sequence, read_position_including_soft_clipped_ + si) == 'N') {
                    break;
                }
                const int base_quality = qual_at(query_quality, read_position_including_soft_clipped_ + si);
                if (base_quality <= 12) {
                    number_of_low_quality_bases++;
                }
                if (number_of_low_quality_bases > 1) {
                    break;
                }
                sum_of_read_qualities += base_quality;
                number_of_high_quality_bases++;
            }
            sclip3_high_quality_processing(query_sequence, mapping_quality, query_quality, number_of_mismatches,
                                           direction, sum_of_read_qualities, number_of_high_quality_bases,
                                           number_of_low_quality_bases);
        }
    }
    read_position_including_soft_clipped_ += cigar_element_length_;
    offset_ = 0;
    start_ = position;  // Had to reset the start due to softclipping adjustment
}

void CigarParser::process_not_matched() {
    const std::string key = std::to_string(start_ - 1) + "-" + std::to_string(start_ + cigar_element_length_ - 1);
    scope_.splice.insert(key);
    inc_cnt(data_.splice_count, key, 1);
    start_ += cigar_element_length_;
    offset_ = 0;
}

bool CigarParser::parse_cigar_with_amp_case(const io::BamRecord& record, bool is_mate_reference_name_equal) const {
    static const util::Regex colon = util::Regex::compile(":").value();
    auto split = util::regex_split(*scope_.run.amplicon_based_calling, colon);
    int distance_to_amplicon = 10;
    double overlap_fraction = 0.95;
    {
        auto d = split.empty() ? std::expected<int, std::string>(std::unexpected("")) : util::parse_int(split[0]);
        auto f = util::parse_double(split.size() > 1 ? split[1] : std::string_view{});
        if (d && f) {
            distance_to_amplicon = *d;
            overlap_fraction = *f;
        }
    }
    const int read_length_include_matched_and_deleted = aligned_length(cigar_);
    int segstart = record.alignment_start();
    int segend = segstart + read_length_include_matched_and_deleted - 1;
    const Region& region = scope_.region;

    if (element(0).op == 'S') {
        const int ts1 = segstart > region.start ? segstart : region.start;
        const int te1 = segend < region.end ? segend : region.end;
        if ((std::abs(ts1 - te1) / static_cast<double>(segend - segstart) > overlap_fraction) == false) {
            return true;
        }
    } else if (element(static_cast<int>(cigar_.size()) - 1).op == 'S') {
        const int ts1 = segstart > region.start ? segstart : region.start;
        const int te1 = segend < region.end ? segend : region.end;
        if ((std::abs(te1 - ts1) / static_cast<double>(segend - segstart) > overlap_fraction) == false) {
            return true;
        }
    } else {
        const int isize = static_cast<int>(record.raw()->core.isize);
        if (is_mate_reference_name_equal && isize != 0) {
            if (isize > 0) {
                segend = segstart + isize - 1;
            } else {
                segstart = record.mate_alignment_start();
                segend = record.mate_alignment_start() - isize - 1;
            }
        }
        const int ts1 = segstart > region.start ? segstart : region.start;
        const int te1 = segend < region.end ? segend : region.end;
        if ((std::abs(segstart - region.start) > distance_to_amplicon ||
             std::abs(segend - region.end) > distance_to_amplicon) ||
            std::abs((ts1 - te1) / static_cast<double>(segend - segstart)) <= overlap_fraction) {
            return true;
        }
    }
    return false;
}

void CigarParser::sub_cnt(Variation& v, bool dir, int read_position, double base_quality, int mapping_quality,
                          int nm) const {
    v.vars_count--;
    v.dec_dir(dir);
    v.mean_position -= read_position;
    v.mean_quality -= base_quality;
    v.mean_mapping_quality -= mapping_quality;
    v.number_of_mismatches -= nm;
    if (base_quality >= conf_.goodq) {
        v.high_quality_reads_count--;
    } else {
        v.low_quality_reads_count--;
    }
}

void CigarParser::add_cnt(Variation& v, bool dir, int read_position, double base_quality, int mapping_quality,
                          int nm) const {
    v.vars_count++;
    v.inc_dir(dir);
    v.mean_position += read_position;
    v.mean_quality += base_quality;
    v.mean_mapping_quality += mapping_quality;
    v.number_of_mismatches += nm;
    if (base_quality >= conf_.goodq) {
        v.high_quality_reads_count++;
    } else {
        v.low_quality_reads_count++;
    }
}

CigarParser::Offset CigarParser::find_offset(int reference_position, int read_position, int cigar_length,
                                             const std::string& query_sequence, const std::string& query_quality) {
    int offset = 0;
    std::string ss;
    std::string q;
    int tnm = 0;
    int vsn = 0;
    for (int vi = 0; vsn <= conf_.vext && vi < cigar_length; vi++) {
        if (char_at(query_sequence, read_position + vi) == 'N') {
            break;
        }
        if (qual_at(query_quality, read_position + vi) < conf_.goodq) {
            break;
        }
        if (auto ref_ch = ref_get(ref_, reference_position + vi)) {
            const char ch = char_at(query_sequence, read_position + vi);
            if (is_not_equals(ch, ref_ch)) {
                offset = vi + 1;
                tnm++;
                vsn = 0;
            } else {
                vsn++;
            }
        }
    }
    if (offset > 0) {
        ss = substr(query_sequence, read_position, offset);
        q = substr(query_quality, read_position, offset);
        for (int osi = 0; osi < offset; osi++) {
            inc_cnt(data_.ref_coverage, reference_position + osi, 1);
        }
    }
    return Offset{offset, ss, q, tnm};
}

bool CigarParser::skip_indel_next_to_intron(int ci) const {
    return (static_cast<int>(cigar_.size()) > ci + 1 && element(ci + 1).op == 'N') ||
           (ci > 0 && element(ci - 1).op == 'N');
}

void CigarParser::add_variation_for_matching_part(int mapping_quality, int nm, bool dir, int rlen1, int nmoff,
                                                  const std::string& s, bool start_with_deletion, double q, int qbases,
                                                  int qibases, int ddlen, int pos) {
    Variation* hv_ptr;
    if (s.starts_with("+")) {
        inc_cnt(data_.position_to_insertion_count.get_or_create(pos), s, 1);
        hv_ptr = &get_variation(data_.insertion_variants, pos, s);
    } else {
        // isBEGIN_ATGC_AMP_ATGCs_END: one base, '&', one or more bases
        if (s.size() > 2 && s[1] == '&' && patterns::is_atgc(s[0]) &&
            patterns::begin_atgc_end(std::string_view(s).substr(2))) {
            inc_cnt(data_.mnp.get_or_create(pos), s, 1);
        }
        hv_ptr = &get_variation(data_.non_insertion_variants, pos, s);
    }
    Variation& hv = *hv_ptr;
    hv.inc_dir(dir);
    hv.vars_count++;
    const int rpes = read_position_excluding_soft_clipped_;
    const int tp = rpes < rlen1 - rpes ? rpes + 1 : rlen1 - rpes;
    q = q / (qbases + qibases);
    if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
        hv.pstd = true;
    }
    if (!hv.qstd && hv.pq != 0 && q != hv.pq) {
        hv.qstd = true;
    }
    hv.mean_position += tp;
    hv.mean_quality += q;
    hv.mean_mapping_quality += mapping_quality;
    hv.pp = tp;
    hv.pq = q;
    hv.number_of_mismatches += nm - nmoff;
    if (q >= conf_.goodq) {
        hv.high_quality_reads_count++;
    } else {
        hv.low_quality_reads_count++;
    }
    const int shift = (s.starts_with("+") && s.contains('&')) ? 1 : 0;
    for (int qi = 1; qi <= qbases - shift; qi++) {
        inc_cnt(data_.ref_coverage, start_ - qi + 1, 1);
    }
    if (start_with_deletion) {
        inc_cnt(data_.position_to_deletion_count.get_or_create(pos), s, 1);
        for (int qi = 1; qi < ddlen; qi++) {
            inc_cnt(data_.ref_coverage, start_ + qi, 1);
        }
    }
}

void CigarParser::add_variation_for_deletion(int mapping_quality, int nm, bool dir, int rlen1, const std::string& desc,
                                             const std::string& quality_of_segment, int nmoff) {
    Variation& hv = get_variation(data_.non_insertion_variants, start_, desc);
    inc_cnt(data_.position_to_deletion_count.get_or_create(start_), desc, 1);
    hv.inc_dir(dir);
    hv.vars_count++;
    const int rpes = read_position_excluding_soft_clipped_;
    const int tp = rpes < rlen1 - rpes ? rpes + 1 : rlen1 - rpes;
    double tmpq = 0;
    for (char c : quality_of_segment) {
        tmpq += static_cast<unsigned char>(c) - 33;
    }
    tmpq = tmpq / static_cast<double>(quality_of_segment.size());
    if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
        hv.pstd = true;
    }
    if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
        hv.qstd = true;
    }
    hv.mean_position += tp;
    hv.mean_quality += tmpq;
    hv.mean_mapping_quality += mapping_quality;
    hv.pp = tp;
    hv.pq = tmpq;
    hv.number_of_mismatches += nm - nmoff;
    if (tmpq >= conf_.goodq) {
        hv.high_quality_reads_count++;
    } else {
        hv.low_quality_reads_count++;
    }
    for (int i = 0; i < cigar_element_length_; i++) {
        inc_cnt(data_.ref_coverage, start_ + i, 1);
    }
}

bool CigarParser::is_read_chimeric_with_sa(const io::BamRecord& record, const io::BamFile& file, int position,
                                           const std::string& sa_tag, bool dir, bool is5side) const {
    auto f = util::regex_split(sa_tag, comma());
    auto field = [&](std::size_t i) -> std::string {
        if (i >= f.size()) {
            throw util::RecoverableError("ArrayIndexOutOfBoundsException: Index " + std::to_string(i) +
                                         " out of bounds for length " + std::to_string(f.size()));
        }
        return std::string(f[i]);
    };
    const std::string sa_chromosome = field(0);
    auto sa_pos = util::parse_int(field(1));
    if (!sa_pos) {
        throw util::RecoverableError("NumberFormatException: " + sa_pos.error());
    }
    const std::string sa_direction = field(2);
    const std::string sa_cigar = field(3);
    const bool sa_direction_is_forward = sa_direction == "+";
    bool clip_match;
    if (is5side) {
        // SA_CIGAR_D_S_5clip "^\d\d+S"
        std::size_t d = 0;
        while (d < sa_cigar.size() && patterns::is_digit(sa_cigar[d])) {
            ++d;
        }
        clip_match = d >= 2 && d < sa_cigar.size() && sa_cigar[d] == 'S';
    } else {
        // SA_CIGAR_D_S_3clip "\d\dS$"
        const std::size_t n = sa_cigar.size();
        clip_match = n >= 3 && sa_cigar[n - 1] == 'S' && patterns::is_digit(sa_cigar[n - 2]) &&
                     patterns::is_digit(sa_cigar[n - 3]);
    }
    const bool chimeric = ((dir && sa_direction_is_forward) || (!dir && !sa_direction_is_forward)) &&
                          sa_chromosome == file.reference_name(record.reference_index()) &&
                          std::abs(*sa_pos - position) < 2 * data_.max_read_length && clip_match;
    if (conf_.y && chimeric) {
        std::println(stderr, "{} {} {} {} {} is ignored as chimeric with SA: {},{},{}", record.read_name(),
                     file.reference_name(record.reference_index()), position, record.mapping_quality(),
                     cigar_to_string(cigar_), *sa_pos, sa_direction, sa_cigar);
    }
    return chimeric;
}

void CigarParser::append_segments(const std::string& query_sequence, const std::string& query_quality, int ci,
                                  std::string& desc, std::string& quality_segment, int m_len, int indel_len, int begin,
                                  bool is_insertion) {
    desc += "#" + substr(query_sequence, begin, m_len);
    quality_segment += substr(query_quality, begin, m_len);
    desc += '^';
    desc += element(ci + 2).op == 'I' ? substr(query_sequence, begin + m_len, indel_len) : std::to_string(indel_len);
    if (is_insertion) {
        if (element(ci + 2).op == 'I') {
            quality_segment += substr(query_quality, begin + m_len, indel_len);
        } else {
            quality_segment += char_at(query_quality, begin + m_len);
        }
    } else {
        if (element(ci + 2).op == 'I') {
            quality_segment += substr(query_quality, begin + m_len, indel_len);
        }
    }
}

bool CigarParser::skip_overlapping_reads(const io::BamRecord& record, const io::BamFile& file, int position, bool dir,
                                         int mate_alignment_start) const {
    auto paired_same_chr = [&] { return record.read_paired() && file.mate_reference_name_or_equals(record) == "="; };
    if (conf_.unique_mode_alignment_enabled && paired_same_chr() && !dir && start_ >= mate_alignment_start) {
        return true;
    }
    if (conf_.unique_mode_second_in_pair_enabled) {
        // htsjdk getSecondOfPairFlag throws IllegalStateException for unpaired reads
        if (!record.read_paired()) {
            throw util::RecoverableError("IllegalStateException: Inappropriate call if not paired read");
        }
        if ((record.flags() & BAM_FREAD2) != 0 && paired_same_chr()) {
            bool overlap;
            if (position >= mate_alignment_start) {
                overlap =
                    start_ >= mate_alignment_start && start_ <= (mate_alignment_start + record.reference_length() - 1);
            } else {
                overlap = start_ >= mate_alignment_start && record.mate_alignment_start() <= record.alignment_end();
            }
            if (overlap) {
                return true;
            }
        }
    }
    return false;
}

void CigarParser::sclip3_high_quality_processing(const std::string& query_sequence, int mapping_quality,
                                                 const std::string& query_quality, int number_of_mismatches, bool dir,
                                                 int sum_of_read_qualities, int number_of_high_quality_bases,
                                                 int number_of_low_quality_bases) {
    if (number_of_high_quality_bases >= 1 && number_of_high_quality_bases > number_of_low_quality_bases &&
        start_ >= scope_.region.start && start_ <= scope_.region.end) {
        Sclip& sclip = data_.soft_clips_3end.get_or_create(start_);
        for (int si = 0; si < number_of_high_quality_bases; si++) {
            const char ch = char_at(query_sequence, read_position_including_soft_clipped_ + si);
            const int idx = si;
            inc_cnt_tree(sclip.nt[idx], ch);
            Variation& variation = get_variation_from_seq(sclip, idx, ch);
            add_cnt(variation, dir, number_of_high_quality_bases - si,
                    qual_at(query_quality, read_position_including_soft_clipped_ + si), mapping_quality,
                    number_of_mismatches);
        }
        add_cnt(sclip, dir, cigar_element_length_,
                sum_of_read_qualities / static_cast<double>(number_of_high_quality_bases), mapping_quality,
                number_of_mismatches);
    }
}

void CigarParser::sclip5_high_quality_processing(const std::string& query_sequence, int mapping_quality,
                                                 const std::string& query_quality, int number_of_mismatches, bool dir,
                                                 int sum_of_read_qualities, int number_of_high_quality_bases,
                                                 int number_of_low_quality_bases) {
    if (number_of_high_quality_bases >= 1 && number_of_high_quality_bases > number_of_low_quality_bases &&
        start_ >= scope_.region.start && start_ <= scope_.region.end) {
        Sclip& sclip = data_.soft_clips_5end.get_or_create(start_);
        for (int si = cigar_element_length_ - 1; cigar_element_length_ - si <= number_of_high_quality_bases; si--) {
            const char ch = char_at(query_sequence, si);
            const int idx = cigar_element_length_ - 1 - si;
            inc_cnt_tree(sclip.nt[idx], ch);
            Variation& variation = get_variation_from_seq(sclip, idx, ch);
            add_cnt(variation, dir, si - (cigar_element_length_ - number_of_high_quality_bases),
                    qual_at(query_quality, si), mapping_quality, number_of_mismatches);
        }
        add_cnt(sclip, dir, cigar_element_length_,
                sum_of_read_qualities / static_cast<double>(number_of_high_quality_bases), mapping_quality,
                number_of_mismatches);
    }
}

}  // namespace senri::core
