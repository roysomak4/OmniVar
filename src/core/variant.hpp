#pragma once

// Vars.variants and Vars.varDescriptionStringToVariants, hence shared_ptr.

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cli/config.hpp"
#include "util/ordered_maps.hpp"

namespace senri::core {

// Why isGoodVar rejects a variant: one value per rule, in the order the rules are checked.
// -p writes rejected variants with FILTER = reject_id(reason).
enum class RejectReason {
    no_ref,       // empty reference allele
    af,           // frequency < -f
    min_reads,    // hicnt < -r
    pmean,        // meanPosition < -P
    qual,         // meanQuality < -q
    ref_mapq,     // AF < 0.25 and the reference allele has much better mean MAPQ
    splice,       // deletion that is a splice junction
    qratio,       // highQualityToLowQualityRatio < -o
    mapq,         // AF <= 0.30 and meanMappingQuality < -O
    msi_monomer,  // AF <= 0.30, monomer MSI >= 15 and AF <= -mfreq
    msi,          // AF <= 0.30, non-monomer MSI >= 12 and AF <= -nmfreq
    strand_bias,  // AF < 0.20, strand bias flag 2;1, SNV or short indel
};

// VCF FILTER ID of a reason (REJ_*; distinct from the var2vcf annotation tags).
[[nodiscard]] std::string_view reject_id(RejectReason reason);
// Every reason with a description of its rule and the thresholds from `conf` (##FILTER lines).
[[nodiscard]] std::vector<std::pair<std::string, std::string>> reject_filter_descriptions(const Config& conf);

struct Variant {
    std::string description_string;  // descriptionString
    int position_coverage = 0;       // positionCoverage
    int vars_count_on_forward = 0;
    int vars_count_on_reverse = 0;
    std::string strand_bias_flag = "0";
    double frequency = 0;
    double mean_position = 0;
    bool is_at_least_at_2_positions = false;
    double mean_quality = 0;
    bool has_at_least_2_diff_qualities = false;
    double mean_mapping_quality = 0;
    double high_quality_to_low_quality_ratio = 0;
    double high_quality_reads_frequency = 0;
    double extra_frequency = 0;
    int shift3 = 0;
    double msi = 0;
    int msint = 0;
    double number_of_mismatches = 0;
    int hicnt = 0;
    int hicov = 0;
    std::string leftseq;
    std::string rightseq;
    int start_position = 0;
    int end_position = 0;
    int ref_reverse_coverage = 0;
    int ref_forward_coverage = 0;
    int total_pos_coverage = 0;
    double duprate = 0;
    std::optional<std::string> genotype;  // null until set
    std::string varallele;
    std::string refallele;
    std::string vartype;
    std::string debug;
    int crispr = 0;

    bool is_noise(const Config& conf);
    void adj_complex();
    [[nodiscard]] std::string var_type() const;
    [[nodiscard]] bool is_good_var(const Variant* reference_var, std::string type,
                                   const std::unordered_set<std::string>& splice, const Config& conf) const;
    // isGoodVar's first failing rule (std::nullopt: the variant is good), for -p (rejected variants).
    [[nodiscard]] std::optional<RejectReason> reject_reason(const Variant* reference_var, std::string type,
                                                            const std::unordered_set<std::string>& splice,
                                                            const Config& conf) const;
};

using VariantPtr = std::shared_ptr<Variant>;

struct Vars {
    VariantPtr reference_variant;
    std::vector<VariantPtr> variants;
    util::HashOrderMap<std::string, VariantPtr> var_description_string_to_variants;
    std::string sv;
};

using AlignedVars = util::HashOrderMap<int, Vars>;

// VariationUtils.getOrPutVars
inline Vars& get_or_put_vars(AlignedVars& map, int position) { return map.get_or_create(position); }

}  // namespace senri::core
