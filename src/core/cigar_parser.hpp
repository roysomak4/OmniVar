#pragma once

#include <string>

#include "core/cigar.hpp"
#include "core/scope.hpp"
#include "core/variation.hpp"
#include "core/variation_utils.hpp"
#include "io/record_preprocessor.hpp"

namespace senri::core {

struct CigarParserResult {
    VariationData data;
    // Splicing mode (-i): "sample\tchr\tintron\tcount" lines, in spliceCount iteration order.
    std::string splicing_output;
};

class CigarParser {
public:
    // svflag: the partial pipeline variant (duprate computed differently).
    CigarParser(const RegionScope& scope, bool svflag, int max_read_length);

    // CigarParser.process: consumes all reads of the preprocessor.
    [[nodiscard]] CigarParserResult process(io::RecordPreprocessor& processor);

    // AbstractMode.partialPipeline: CigarParser(true).process on a scope
    // whose InitialData is `shared`'s nonInsertionVariants, insertionVariants, refCoverage,
    // softClips3End and softClips5End. The reads are added to those maps (objects already in
    // them keep their addresses); everything else the parser builds is discarded.
    void process_into(io::RecordPreprocessor& processor, VariationData& shared);

    // CigarParser.parseCigar for one record (exposed for tests).
    void parse_cigar(const std::string& chr_name, const io::BamRecord& record, const io::BamFile& file);

private:
    struct Offset {
        int offset;
        std::string sequence;
        std::string quality_sequence;
        int offset_number_of_mismatches;
    };

    void process_not_matched();
    void prepare_sv_structures_for_analysis(const io::BamRecord& record, const io::BamFile& file,
                                            const std::string& query_quality, int number_of_mismatches,
                                            bool read_direction, bool mate_direction, int start,
                                            int total_length_including_soft_clipped);
    void add_sv(Sclip& sdref, int start_s, int end_e, int mate_start_ms, int mate_end_me, int dir, int rlen, int mlen,
                int softp, double pmean_rp, double qmean, double big_qmean, double nm) const;
    void process_soft_clip(const std::string& chr_name, const io::BamRecord& record, const io::BamFile& file,
                           const std::string& query_sequence, int mapping_quality, const std::string& query_quality,
                           int number_of_mismatches, bool direction, int position,
                           int total_length_including_soft_clipped, int ci);
    int process_insertion(const std::string& query_sequence, int mapping_quality, const std::string& query_quality,
                          int number_of_mismatches, bool direction, int position,
                          int read_length_include_matching_and_insertions, int ci);
    int process_deletion(const std::string& query_sequence, int mapping_quality, const std::string& query_quality,
                         int number_of_mismatches, bool direction, int read_length_include_matching_and_insertions,
                         int ci);
    void add_variation_for_matching_part(int mapping_quality, int nm, bool dir, int rlen1, int nmoff,
                                         const std::string& s, bool start_with_deletion, double q, int qbases,
                                         int qibases, int ddlen, int pos);
    void add_variation_for_deletion(int mapping_quality, int nm, bool dir, int rlen1, const std::string& desc,
                                    const std::string& quality_of_segment, int nmoff);
    void sclip3_high_quality_processing(const std::string& query_sequence, int mapping_quality,
                                        const std::string& query_quality, int number_of_mismatches, bool dir,
                                        int sum_of_read_qualities, int number_of_high_quality_bases,
                                        int number_of_low_quality_bases);
    void sclip5_high_quality_processing(const std::string& query_sequence, int mapping_quality,
                                        const std::string& query_quality, int number_of_mismatches, bool dir,
                                        int sum_of_read_qualities, int number_of_high_quality_bases,
                                        int number_of_low_quality_bases);
    void append_segments(const std::string& query_sequence, const std::string& query_quality, int ci, std::string& desc,
                         std::string& quality_segment, int m_len, int indel_len, int begin, bool is_insertion);
    Offset find_offset(int reference_position, int read_position, int cigar_length, const std::string& query_sequence,
                       const std::string& query_quality);

    bool is_read_chimeric_with_sa(const io::BamRecord& record, const io::BamFile& file, int position,
                                  const std::string& sa_tag, bool dir, bool is5side) const;
    bool skip_overlapping_reads(const io::BamRecord& record, const io::BamFile& file, int position, bool dir,
                                int mate_alignment_start) const;
    bool skip_sites_out_region_of_interest() const;
    bool is_trim_at_opt_t_bases(bool direction, int total_length_including_soft_clipped) const;
    bool parse_cigar_with_amp_case(const io::BamRecord& record, bool is_mate_reference_name_equal) const;
    bool is_closer_then_vext_and_good_base(const std::string& query_sequence, const std::string& query_quality, int ci,
                                           int i, const std::string& ss, char op) const;

    void add_cnt(Variation& v, bool dir, int read_position, double base_quality, int mapping_quality, int nm) const;
    void sub_cnt(Variation& v, bool dir, int read_position, double base_quality, int mapping_quality, int nm) const;

    // cigar.getCigarElement(i) with the upstream IndexOutOfBoundsException
    const CigarElement& element(int i) const;
    bool next_after_num_matched(int ci, int number) const;
    bool two_insertions_ahead(int ci) const;
    bool next_insertion(int ci) const;
    bool next_matched(int ci) const;
    bool insertion_or_deletion_with_next_matched(int ci) const;
    bool skip_indel_next_to_intron(int ci) const;

    const RegionScope& scope_;
    const Config& conf_;
    const RefMap& ref_;
    VariationData data_;
    bool svflag_;

    Cigar cigar_;
    int start_ = 0;
    int discordant_count_ = 0;
    int offset_ = 0;
    int cigar_element_length_ = 0;
    int read_position_including_soft_clipped_ = 0;
    int read_position_excluding_soft_clipped_ = 0;
};

}  // namespace senri::core
