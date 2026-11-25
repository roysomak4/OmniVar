#pragma once

#include <string>

#include "cli/config.hpp"
#include "core/scope.hpp"
#include "io/reference.hpp"
#include "region/region.hpp"

namespace omnivar::core {

struct ModifiedCigar {
    int position;
    std::string cigar;
    std::string query_sequence;
    std::string query_quality;
};

// new CigarModifier(position, cigarStr, querySequence, queryQuality, reference, indel, region,
// maxReadLength).modifyCigar(). An exception inside (any exception) is reported with
// Utils.printExceptionAndContinue and the values as modified so far are returned.
[[nodiscard]] ModifiedCigar modify_cigar(int position, const std::string& cigar, const std::string& query_sequence,
                                         const std::string& query_quality, const io::Reference& reference,
                                         int insertion_deletion_length, const Region& region, int max_read_length,
                                         const Config& conf, ExceptionCounter& exceptions);

}  // namespace omnivar::core
