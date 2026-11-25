#pragma once

// Indexed FASTA through htslib faidx, with the bounds checks of htsjdk
// IndexedFastaSequenceFile.getSubsequenceAt (used by ReferenceResource.retrieveSubSeq).

#include <memory>
#include <string>

#include "io/reference_provider.hpp"

struct faidx_t;

namespace omnivar::io {

class FaidxProvider final : public ReferenceProvider {
public:
    [[nodiscard]] static std::expected<FaidxProvider, ReferenceError> open(const std::string& path);

    [[nodiscard]] std::expected<std::string, ReferenceError> sub_seq(std::string_view chr, int start,
                                                                     int end) const override;

private:
    struct Free {
        void operator()(faidx_t* f) const;
    };
    std::string path_;
    std::shared_ptr<faidx_t> fai_;
};

}  // namespace omnivar::io
