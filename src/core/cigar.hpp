#pragma once

// CIGAR as htsjdk's Cigar: a list of (length, operator) with operator characters
// M I D N S H P = X. Helpers port the static CigarParser length functions.

#include <string>
#include <string_view>
#include <vector>

#include "io/bam.hpp"

namespace senri::core {

struct CigarElement {
    int length;
    char op;
};

using Cigar = std::vector<CigarElement>;

[[nodiscard]] Cigar cigar_of(const io::BamRecord& record);

// Cigar.toString (TextCigarCodec.encode): "*" when empty.
[[nodiscard]] std::string cigar_to_string(const Cigar& c);

// TextCigarCodec.decode for well-formed text ("*" -> empty).
[[nodiscard]] Cigar parse_cigar_text(std::string_view text);

// CigarOperator.consumesReadBases / consumesReferenceBases
[[nodiscard]] constexpr bool consumes_read(char op) {
    return op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X';
}
[[nodiscard]] constexpr bool consumes_reference(char op) {
    return op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X';
}

// CigarParser.getAlignedLength: M + D
[[nodiscard]] int aligned_length(const Cigar& c);
// CigarParser.getSoftClippedLength: M + I + S
[[nodiscard]] int soft_clipped_length(const Cigar& c);
// CigarParser.getMatchInsertionLength: M + I
[[nodiscard]] int match_insertion_length(const Cigar& c);
// CigarParser.getInsertionDeletionLength: I + D
[[nodiscard]] int insertion_deletion_length(const Cigar& c);
// Cigar.getReferenceLength: M D N = X
[[nodiscard]] int reference_length(const Cigar& c);

}  // namespace senri::core
