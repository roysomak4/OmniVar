#pragma once

// BAM access through htslib, with accessors that behave like htsjdk's SAMRecord, and a port of

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdint>
#include <expected>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#include "region/region.hpp"

namespace omnivar::io {

struct BamError {
    std::string message;
};

// One alignment record. Accessor names follow htsjdk SAMRecord.
class BamRecord {
public:
    BamRecord();

    bam1_t* raw() { return b_.get(); }
    [[nodiscard]] const bam1_t* raw() const { return b_.get(); }

    [[nodiscard]] std::string_view read_name() const { return bam_get_qname(b_.get()); }
    [[nodiscard]] int flags() const { return b_->core.flag; }
    [[nodiscard]] int mapping_quality() const { return b_->core.qual; }
    [[nodiscard]] int reference_index() const { return b_->core.tid; }
    [[nodiscard]] int mate_reference_index() const { return b_->core.mtid; }
    // 1-based; 0 when unset (SAMRecord.NO_ALIGNMENT_START)
    [[nodiscard]] int alignment_start() const { return static_cast<int>(b_->core.pos + 1); }
    [[nodiscard]] int mate_alignment_start() const { return static_cast<int>(b_->core.mpos + 1); }
    // SAMRecord.getInferredInsertSize (TLEN)
    [[nodiscard]] int inferred_insert_size() const { return static_cast<int>(b_->core.isize); }
    // SAMRecord.getAlignmentEnd: start + reference length - 1 (start - 1 for an empty CIGAR); 0 if unmapped.
    [[nodiscard]] int alignment_end() const;
    // Reference length of the CIGAR (M, D, N, =, X).
    [[nodiscard]] int reference_length() const;

    [[nodiscard]] bool read_paired() const { return (flags() & BAM_FPAIRED) != 0; }
    [[nodiscard]] bool read_unmapped() const { return (flags() & BAM_FUNMAP) != 0; }
    [[nodiscard]] bool mate_unmapped() const { return (flags() & BAM_FMUNMAP) != 0; }
    [[nodiscard]] bool read_negative_strand() const { return (flags() & BAM_FREVERSE) != 0; }
    [[nodiscard]] bool secondary() const { return (flags() & BAM_FSECONDARY) != 0; }
    [[nodiscard]] bool supplementary() const { return (flags() & BAM_FSUPPLEMENTARY) != 0; }

    [[nodiscard]] int read_length() const { return b_->core.l_qseq; }
    // SAMRecord.getReadString: "*" when SEQ is absent.
    [[nodiscard]] std::string read_string() const;
    // SAMRecord.getBaseQualityString: phred+33, "*" when absent.
    [[nodiscard]] std::string base_quality_string() const;
    // SAMRecord.getCigarString: "*" when there is no CIGAR.
    [[nodiscard]] std::string cigar_string() const;
    // SAMRecord.getIntegerAttribute: nullopt when absent or not an integer type.
    [[nodiscard]] std::optional<std::int64_t> integer_attribute(const char tag[2]) const;
    // SAMRecord.getStringAttribute: nullopt when absent or not a string ('Z') type.
    [[nodiscard]] std::optional<std::string_view> string_attribute(const char tag[2]) const;

private:
    struct Free {
        void operator()(bam1_t* b) const { bam_destroy1(b); }
    };
    std::unique_ptr<bam1_t, Free> b_;
};

// An open, indexed BAM/CRAM file.
class BamFile {
public:
    [[nodiscard]] static std::expected<BamFile, BamError> open(const std::string& path);

    [[nodiscard]] const sam_hdr_t* header() const { return hdr_.get(); }
    [[nodiscard]] const std::string& path() const { return path_; }

    // Contig names and lengths from the header.
    [[nodiscard]] ChrLengths chr_lengths() const;
    // The same, in header order (VCF ##contig lines).
    [[nodiscard]] std::vector<std::pair<std::string, long long>> contigs() const;

    // SAMRecord.getReferenceName / getMateReferenceName: "*" when the index is -1.
    [[nodiscard]] std::string_view reference_name(int tid) const;

    // CigarParser.getMateReferenceName: "=" if same as the read's contig.
    [[nodiscard]] std::string_view mate_reference_name_or_equals(const BamRecord& r) const;

    struct Impl;

private:
    struct HdrFree {
        void operator()(sam_hdr_t* h) const { sam_hdr_destroy(h); }
    };
    std::string path_;
    std::shared_ptr<Impl> impl_;  // htsFile + index, shared by queries
    std::shared_ptr<sam_hdr_t> hdr_;

    friend class SamView;
};

class SamView {
public:
    [[nodiscard]] static std::expected<SamView, BamError> query(const BamFile& file, const std::string& samfilter,
                                                                const Region& region);

    // SamView.read: next record, or false at the end. Errors from htslib end the iteration with
    // an error message available from error().
    bool read(BamRecord& out);

    [[nodiscard]] const std::optional<std::string>& error() const { return error_; }

private:
    struct ItrFree {
        void operator()(hts_itr_t* it) const { hts_itr_destroy(it); }
    };
    std::shared_ptr<BamFile::Impl> impl_;
    std::unique_ptr<hts_itr_t, ItrFree> itr_;
    int filter_ = 0;
    int qstart_ = 0;
    int qend_ = 0;
    bool empty_ = false;
    std::optional<std::string> error_;
};

}  // namespace omnivar::io
