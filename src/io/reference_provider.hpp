#pragma once

// Source of reference bases. Mirrors ReferenceResource.retrieveSubSeq
// Callers clamp coordinates (ReferenceResource.getReference clamps to [1, contig length]);
// providers only look the range up.

#include <expected>
#include <string>
#include <string_view>

namespace senri::io {

struct ReferenceError {
    enum class Kind {
        missing_contig,  // upstream: WrongFastaOrBamException
        out_of_bounds,   // upstream: RegionBoundariesException
        no_slice,        // thin-FASTA CSV has no slice covering the range
        io,              // file could not be opened or read
    };
    Kind kind;
    std::string message;
};

class ReferenceProvider {
public:
    virtual ~ReferenceProvider() = default;

    // Bases for chr:start-end, 1-based inclusive, as stored in the file (case preserved;
    // the caller upper-cases).
    [[nodiscard]] virtual std::expected<std::string, ReferenceError> sub_seq(std::string_view chr, int start,
                                                                             int end) const = 0;
};

}  // namespace senri::io
