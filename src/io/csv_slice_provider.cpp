#include "io/csv_slice_provider.hpp"

#include <charconv>
#include <fstream>
#include <sstream>

namespace omnivar::io {

namespace {

ReferenceError parse_error(const std::string& path, std::size_t line_no, std::string_view what) {
    return {ReferenceError::Kind::io,
            "Malformed reference CSV " + path + " line " + std::to_string(line_no) + ": " + std::string(what)};
}

}  // namespace

bool is_csv_reference(std::string_view path) { return path.ends_with(".csv"); }

std::expected<CsvSliceProvider, ReferenceError> CsvSliceProvider::open(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return std::unexpected(ReferenceError{ReferenceError::Kind::io, "Couldn't open reference file: " + path});
    }
    std::ostringstream buf;
    buf << in.rdbuf();
    return from_text(buf.str(), path);
}

std::expected<CsvSliceProvider, ReferenceError> CsvSliceProvider::from_text(std::string_view text, std::string path) {
    CsvSliceProvider p;
    p.path_ = std::move(path);
    std::size_t line_no = 0;
    while (!text.empty()) {
        auto nl = text.find('\n');
        std::string_view line = text.substr(0, nl);
        text = nl == std::string_view::npos ? std::string_view{} : text.substr(nl + 1);
        ++line_no;
        if (line.ends_with('\r')) {  // \r\n line ends too
            line.remove_suffix(1);
        }
        if (line.empty()) {
            continue;
        }
        // Fields: chr, start, end (ignored), sequence.
        std::string_view f[4];
        std::size_t n = 0;
        while (n < 3) {
            auto comma = line.find(',');
            if (comma == std::string_view::npos) {
                break;
            }
            f[n++] = line.substr(0, comma);
            line.remove_prefix(comma + 1);
        }
        if (n < 3) {
            return std::unexpected(parse_error(p.path_, line_no, "expected chr,start,end,sequence"));
        }
        f[3] = line.substr(0, line.find(','));
        int start = 0;
        auto [ptr, ec] = std::from_chars(f[1].data(), f[1].data() + f[1].size(), start);
        if (ec != std::errc{} || ptr != f[1].data() + f[1].size()) {
            return std::unexpected(parse_error(p.path_, line_no, "bad start"));
        }
        p.slices_[std::string(f[0])][start] = std::string(f[3]);
    }
    return p;
}

std::expected<std::string, ReferenceError> CsvSliceProvider::sub_seq(std::string_view chr, int start, int end) const {
    if (auto it = slices_.find(std::string(chr)); it != slices_.end()) {
        const auto& by_start = it->second;
        auto e = by_start.upper_bound(start);  // floorEntry(start)
        if (e != by_start.begin()) {
            --e;
            const int csv_start = e->first;
            const std::string& seq = e->second;
            const long long csv_end = static_cast<long long>(csv_start) + static_cast<long long>(seq.size()) - 1;
            if (csv_end >= end) {
                if (end < start) {
                    return std::string{};
                }
                return seq.substr(static_cast<std::size_t>(start - csv_start),
                                  static_cast<std::size_t>(end - start + 1));
            }
        }
    }
    return std::unexpected(
        ReferenceError{ReferenceError::Kind::no_slice, "No region: >" + std::string(chr) + ":" + std::to_string(start) +
                                                           "-" + std::to_string(end) + " in file " + path_});
}

}  // namespace omnivar::io
