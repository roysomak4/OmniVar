#pragma once

// Reference served from a "thin FASTA" CSV (docs/golden/fastas/*.csv), one slice per line:
//   chr,start,end,sequence
// Same lookup as the upstream test harness and the test oracle (tests/oracle): take the slice
// with the greatest
// start <= query start (a later line with the same chr and start replaces an earlier one),
// and serve the query only if that slice covers the query end. The CSV end column is ignored;
// the slice end is start + len(sequence) - 1.

#include <map>
#include <string>
#include <string_view>
#include <unordered_map>

#include "io/reference_provider.hpp"

namespace omnivar::io {

class CsvSliceProvider final : public ReferenceProvider {
public:
    [[nodiscard]] static std::expected<CsvSliceProvider, ReferenceError> open(const std::string& path);
    [[nodiscard]] static std::expected<CsvSliceProvider, ReferenceError> from_text(std::string_view text,
                                                                                   std::string path);

    [[nodiscard]] std::expected<std::string, ReferenceError> sub_seq(std::string_view chr, int start,
                                                                     int end) const override;

private:
    std::string path_;
    // chr -> slice start -> sequence
    std::unordered_map<std::string, std::map<int, std::string>> slices_;
};

// True when the -G path names a thin-FASTA CSV rather than an indexed FASTA.
[[nodiscard]] bool is_csv_reference(std::string_view path);

}  // namespace omnivar::io
