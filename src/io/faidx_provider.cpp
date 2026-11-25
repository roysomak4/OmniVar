#include "io/faidx_provider.hpp"

#include <htslib/faidx.h>

#include <cstdlib>

namespace omnivar::io {

void FaidxProvider::Free::operator()(faidx_t* f) const { fai_destroy(f); }

std::expected<FaidxProvider, ReferenceError> FaidxProvider::open(const std::string& path) {
    faidx_t* f = fai_load(path.c_str());
    if (f == nullptr) {
        return std::unexpected(ReferenceError{ReferenceError::Kind::io, "Couldn't open reference file: " + path});
    }
    FaidxProvider p;
    p.path_ = path;
    p.fai_ = std::shared_ptr<faidx_t>(f, Free{});
    return p;
}

std::expected<std::string, ReferenceError> FaidxProvider::sub_seq(std::string_view chr_view, int start, int end) const {
    const std::string chr(chr_view);
    if (faidx_has_seq(fai_.get(), chr.c_str()) == 0) {
        // htsjdk: "Unable to find entry for contig" -> WrongFastaOrBamException
        return std::unexpected(ReferenceError{ReferenceError::Kind::missing_contig,
                                              "Unable to find entry for contig: " + chr + " in " + path_});
    }
    const hts_pos_t len = faidx_seq_len64(fai_.get(), chr.c_str());
    if (end > len) {
        return std::unexpected(ReferenceError{ReferenceError::Kind::out_of_bounds,
                                              "Query asks for data past end of contig. Query contig " + chr +
                                                  " start:" + std::to_string(start) + " stop:" + std::to_string(end) +
                                                  " contigLength:" + std::to_string(len)});
    }
    if (start > end + 1) {
        return std::unexpected(ReferenceError{
            ReferenceError::Kind::out_of_bounds,
            "Malformed query; start point " + std::to_string(start) + " lies after end point " + std::to_string(end)});
    }
    if (start == end + 1) {
        return std::string{};
    }
    hts_pos_t got = 0;
    char* s = faidx_fetch_seq64(fai_.get(), chr.c_str(), start - 1, end - 1, &got);
    if (s == nullptr || got != static_cast<hts_pos_t>(end) - start + 1) {
        std::free(s);
        return std::unexpected(ReferenceError{
            ReferenceError::Kind::io,
            "Error reading " + chr + ":" + std::to_string(start) + "-" + std::to_string(end) + " from " + path_});
    }
    std::string out(s, static_cast<std::size_t>(got));
    std::free(s);
    return out;
}

}  // namespace omnivar::io
