#pragma once

// VCF comparison shared by `vcf` and `panel`: records matched on CHROM, POS, REF, ALT, END;
// QUAL, FILTER, GT and every non-numeric value identical; numeric INFO and FORMAT values equal as
// numbers (relative 1e-12), because R's write.table and Perl reprint numbers while the C++ VCF
// keeps the upstream text.

#include <format>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace testkit {

using RecordKey = std::tuple<std::string, std::string, std::string, std::string, std::optional<std::string>>;
using Record = std::vector<std::string>;

struct Records {
    std::map<RecordKey, std::vector<Record>> by_key;
    std::vector<RecordKey> order;
    [[nodiscard]] std::size_t count() const { return order.size(); }
};

Records parse_records(std::string_view vcf);
bool same_value(std::string_view x, std::string_view y);
// Differences between two records with the same key; accepted R artefacts go to `known`.
std::vector<std::string> record_diffs(const Record& want, const Record& got, std::vector<std::string>& known);
// Every difference between two record sets ("chr:pos ref>alt: ...").
std::vector<std::string> compare_records(const Records& want, const Records& got, std::vector<std::string>& known);
std::string key_label(const RecordKey& k);

// `bcftools view` parse check: nothing when clean (or bcftools is missing), else the first message.
std::optional<std::string> bcftools_check(const std::string& vcf);

// TSV | Rscript teststrandbias.R | perl var2vcf_valid.pl <perl_args>   (paired: testsomatic.R,
// var2vcf_paired.pl). Throws std::runtime_error when a step fails.
std::string chain(std::string_view tsv, int timeout, std::string_view perl_args, bool paired);

}  // namespace testkit
