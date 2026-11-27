#include "io/bam.hpp"

#include <limits>

#include "util/number_parse.hpp"

namespace senri::io {

struct BamFile::Impl {
    htsFile* fp = nullptr;
    hts_idx_t* idx = nullptr;
    ~Impl() {
        if (idx != nullptr) {
            hts_idx_destroy(idx);
        }
        if (fp != nullptr) {
            hts_close(fp);
        }
    }
};

BamRecord::BamRecord() : b_(bam_init1()) {}

int BamRecord::reference_length() const {
    const auto* b = b_.get();
    return static_cast<int>(bam_cigar2rlen(static_cast<int>(b->core.n_cigar), bam_get_cigar(b)));
}

int BamRecord::alignment_end() const {
    if (read_unmapped()) {
        return 0;
    }
    return alignment_start() + reference_length() - 1;
}

std::string BamRecord::read_string() const {
    const auto* b = b_.get();
    const int n = b->core.l_qseq;
    if (n == 0) {
        return "*";
    }
    std::string s(static_cast<std::size_t>(n), 'N');
    const uint8_t* seq = bam_get_seq(b);
    for (int i = 0; i < n; ++i) {
        s[static_cast<std::size_t>(i)] = seq_nt16_str[bam_seqi(seq, i)];
    }
    return s;
}

std::string BamRecord::base_quality_string() const {
    const auto* b = b_.get();
    const int n = b->core.l_qseq;
    const uint8_t* q = bam_get_qual(b);
    if (n == 0 || q[0] == 0xff) {
        return "*";
    }
    std::string s(static_cast<std::size_t>(n), '!');
    for (int i = 0; i < n; ++i) {
        s[static_cast<std::size_t>(i)] = static_cast<char>(q[i] + 33);
    }
    return s;
}

std::string BamRecord::cigar_string() const {
    const auto* b = b_.get();
    const uint32_t n = b->core.n_cigar;
    if (n == 0) {
        return "*";
    }
    std::string s;
    const uint32_t* cig = bam_get_cigar(b);
    for (uint32_t i = 0; i < n; ++i) {
        s += std::to_string(bam_cigar_oplen(cig[i]));
        s += bam_cigar_opchr(cig[i]);
    }
    return s;
}

std::optional<std::int64_t> BamRecord::integer_attribute(const char tag[2]) const {
    const uint8_t* aux = bam_aux_get(b_.get(), tag);
    if (aux == nullptr) {
        return std::nullopt;
    }
    switch (*aux) {
        case 'c':
        case 'C':
        case 's':
        case 'S':
        case 'i':
        case 'I': return bam_aux2i(aux);
        default: return std::nullopt;
    }
}

std::optional<std::string_view> BamRecord::string_attribute(const char tag[2]) const {
    const uint8_t* aux = bam_aux_get(b_.get(), tag);
    if (aux == nullptr || *aux != 'Z') {
        return std::nullopt;
    }
    return std::string_view(bam_aux2Z(aux));
}

std::expected<BamFile, BamError> BamFile::open(const std::string& path) {
    auto impl = std::make_shared<Impl>();
    impl->fp = hts_open(path.c_str(), "r");
    if (impl->fp == nullptr) {
        return std::unexpected(BamError{"FileNotFoundException: " + path + " (No such file or directory)"});
    }
    sam_hdr_t* h = sam_hdr_read(impl->fp);
    if (h == nullptr) {
        return std::unexpected(BamError{"Cannot read BAM header: " + path});
    }
    BamFile f;
    f.path_ = path;
    f.hdr_ = std::shared_ptr<sam_hdr_t>(h, HdrFree{});
    impl->idx = sam_index_load(impl->fp, path.c_str());
    // A missing index is reported when a query is made (htsjdk: "Cannot query stream ... index").
    f.impl_ = std::move(impl);
    return f;
}

ChrLengths BamFile::chr_lengths() const {
    ChrLengths m;
    const int n = sam_hdr_nref(hdr_.get());
    for (int i = 0; i < n; ++i) {
        const hts_pos_t len = sam_hdr_tid2len(hdr_.get(), i);
        m[sam_hdr_tid2name(hdr_.get(), i)] = static_cast<int>(len);
    }
    return m;
}

std::vector<std::pair<std::string, long long>> BamFile::contigs() const {
    std::vector<std::pair<std::string, long long>> out;
    const int n = sam_hdr_nref(hdr_.get());
    for (int i = 0; i < n; ++i) {
        out.emplace_back(sam_hdr_tid2name(hdr_.get(), i), static_cast<long long>(sam_hdr_tid2len(hdr_.get(), i)));
    }
    return out;
}

std::string_view BamFile::reference_name(int tid) const {
    if (tid < 0 || tid >= sam_hdr_nref(hdr_.get())) {
        return "*";
    }
    return sam_hdr_tid2name(hdr_.get(), tid);
}

std::string_view BamFile::mate_reference_name_or_equals(const BamRecord& r) const {
    // htsjdk never returns null for getMateReferenceName ("*" when unset), so the null branch
    // of CigarParser.getMateReferenceName is unreachable.
    const std::string_view mate = reference_name(r.mate_reference_index());
    if (reference_name(r.reference_index()) == mate) {
        return "=";
    }
    return mate;
}

std::expected<SamView, BamError> SamView::query(const BamFile& file, const std::string& samfilter,
                                                const Region& region) {
    SamView v;
    v.impl_ = file.impl_;
    auto filter = util::decode_int(samfilter);
    if (!filter) {
        return std::unexpected(BamError{"NumberFormatException: " + filter.error()});
    }
    v.filter_ = *filter;
    v.qstart_ = region.start;
    // htsjdk QueryInterval: an end <= 0 means "to the end of the contig".
    v.qend_ = region.end <= 0 ? std::numeric_limits<int>::max() : region.end;
    if (v.impl_->idx == nullptr) {
        return std::unexpected(BamError{"Cannot query stream: no index for " + file.path()});
    }
    const int tid = sam_hdr_name2tid(file.hdr_.get(), region.chr.c_str());
    if (tid < 0) {
        // htsjdk: a contig absent from the header gives an empty iterator.
        v.empty_ = true;
        return v;
    }
    // htslib takes a 0-based half-open interval. Ask for [start-2, end) so that records whose
    // htsjdk end is start-1 can be seen and filtered below; the htsjdk rule decides.
    const hts_pos_t beg = std::max<hts_pos_t>(0, static_cast<hts_pos_t>(region.start) - 2);
    const hts_pos_t end = v.qend_ == std::numeric_limits<int>::max() ? HTS_POS_MAX : static_cast<hts_pos_t>(v.qend_);
    v.itr_.reset(sam_itr_queryi(v.impl_->idx, tid, beg, end));
    if (!v.itr_) {
        return std::unexpected(BamError{"Cannot query " + region.print_region() + " in " + file.path()});
    }
    return v;
}

bool SamView::read(BamRecord& out) {
    if (empty_ || error_) {
        return false;
    }
    while (true) {
        const int r = sam_itr_next(impl_->fp, itr_.get(), out.raw());
        if (r == -1) {
            return false;
        }
        if (r < -1) {
            error_ = "htslib error " + std::to_string(r) + " reading records";
            return false;
        }
        // htsjdk BAMQueryMultipleIntervalsIteratorFilter (OVERLAPPING).
        const int rec_start = out.alignment_start();
        if (rec_start > qend_) {
            continue;  // htslib may return records past the region; htsjdk stops before them
        }
        const int rec_end = out.read_unmapped() ? rec_start : out.alignment_end();
        if (rec_end < qstart_) {
            continue;
        }
        // SamView.read: -F filter
        if (filter_ != 0 && (out.flags() & filter_) != 0) {
            continue;
        }
        return true;
    }
}

}  // namespace senri::io
