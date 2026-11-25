#include "core/variant.hpp"

#include <utility>

#include "core/patterns.hpp"
#include "util/checked_string.hpp"
#include "util/num_format.hpp"

namespace omnivar::core {

bool Variant::is_noise(const Config& conf) {
    const double qual = mean_quality;
    if (((qual < 4.5 || (qual < 12 && !has_at_least_2_diff_qualities)) && position_coverage <= 3) ||
        (qual < conf.goodq && frequency < 2 * conf.lofreq && position_coverage <= 1)) {
        total_pos_coverage -= position_coverage;
        position_coverage = 0;
        vars_count_on_forward = 0;
        vars_count_on_reverse = 0;
        frequency = 0;
        high_quality_reads_frequency = 0;
        return true;
    }
    return false;
}

void Variant::adj_complex() {
    std::string ref_allele = refallele;
    std::string var_allele = varallele;
    if (util::char_at(var_allele, 0) == '<') {
        return;
    }
    long n = 0;
    while (static_cast<long>(ref_allele.size()) - n > 1 && static_cast<long>(var_allele.size()) - n > 1 &&
           util::char_at(ref_allele, n) == util::char_at(var_allele, n)) {
        ++n;
    }
    if (n > 0) {
        start_position += static_cast<int>(n);
        refallele = util::substr(ref_allele, n);
        varallele = util::substr(var_allele, n);
        leftseq += util::substr(ref_allele, 0, n);
        leftseq = util::substr(leftseq, n);
    }
    ref_allele = refallele;
    var_allele = varallele;
    n = 1;
    while (static_cast<long>(ref_allele.size()) - n > 0 && static_cast<long>(var_allele.size()) - n > 0 &&
           util::substr(ref_allele, -n, 1) == util::substr(var_allele, -n, 1)) {
        ++n;
    }
    if (n > 1) {
        end_position -= static_cast<int>(n - 1);
        refallele = util::substr(ref_allele, 0, 1 - n);
        varallele = util::substr(var_allele, 0, 1 - n);
        rightseq = util::substr(ref_allele, 1 - n, n - 1) + util::substr(rightseq, 0, 1 - n);
    }
}

std::string Variant::var_type() const {
    if (refallele == varallele && refallele.size() == 1) {
        return "";
    }
    if (refallele.size() == 1 && varallele.size() == 1) {
        return "SNV";
    }
    if (auto sv = patterns::any_sv(varallele)) {
        return *sv;
    }
    if (refallele.empty() || varallele.empty()) {
        return "Complex";
    }
    if (refallele[0] != varallele[0]) {
        return "Complex";
    }
    if (refallele.size() == 1 && varallele.size() > 1 && varallele.starts_with(refallele)) {
        return "Insertion";
    }
    if (refallele.size() > 1 && varallele.size() == 1 && refallele.starts_with(varallele)) {
        return "Deletion";
    }
    return "Complex";
}

std::optional<RejectReason> Variant::reject_reason(const Variant* reference_var, std::string type,
                                                   const std::unordered_set<std::string>& splice,
                                                   const Config& conf) const {
    if (refallele.empty()) {
        return RejectReason::no_ref;
    }
    if (type.empty()) {
        type = var_type();
    }
    // Upstream: one `||` chain; its first true operand is the reason.
    if (frequency < conf.freq) {
        return RejectReason::af;
    }
    if (hicnt < conf.minr) {
        return RejectReason::min_reads;
    }
    if (mean_position < conf.read_pos_filter) {
        return RejectReason::pmean;
    }
    if (mean_quality < conf.goodq) {
        return RejectReason::qual;
    }
    if (reference_var != nullptr && reference_var->hicnt > conf.minr && frequency < 0.25) {
        // The reference allele has much better mean mapq than var allele, thus likely false positives
        const double d =
            mean_mapping_quality + static_cast<double>(refallele.size()) + static_cast<double>(varallele.size());
        const double f = (1 + d) / (reference_var->mean_mapping_quality + 1);
        if ((d - 2 < 5 && reference_var->mean_mapping_quality > 20) || f < 0.25) {
            return RejectReason::ref_mapq;
        }
    }
    if (type == "Deletion" && splice.contains(std::to_string(start_position) + "-" + std::to_string(end_position))) {
        return RejectReason::splice;
    }
    if (high_quality_to_low_quality_ratio < conf.qratio) {
        return RejectReason::qratio;
    }
    if (frequency > 0.30) {
        return std::nullopt;
    }
    if (mean_mapping_quality < conf.mapq) {
        return RejectReason::mapq;
    }
    if (msi >= 15 && frequency <= conf.monomer_msi_frequency && msint == 1) {
        return RejectReason::msi_monomer;
    }
    if (msi >= 12 && frequency <= conf.non_monomer_msi_frequency && msint > 1) {
        return RejectReason::msi;
    }
    if (strand_bias_flag == "2;1" && frequency < 0.20) {
        if (type == "SNV" || (refallele.size() < 3 && varallele.size() < 3)) {
            return RejectReason::strand_bias;
        }
    }
    return std::nullopt;
}

bool Variant::is_good_var(const Variant* reference_var, std::string type, const std::unordered_set<std::string>& splice,
                          const Config& conf) const {
    return !reject_reason(reference_var, std::move(type), splice, conf).has_value();
}

std::string_view reject_id(RejectReason reason) {
    switch (reason) {
        case RejectReason::no_ref: return "REJ_NoRef";
        case RejectReason::af: return "REJ_AF";
        case RejectReason::min_reads: return "REJ_MinReads";
        case RejectReason::pmean: return "REJ_PMean";
        case RejectReason::qual: return "REJ_Qual";
        case RejectReason::ref_mapq: return "REJ_RefMQ";
        case RejectReason::splice: return "REJ_Splice";
        case RejectReason::qratio: return "REJ_QRatio";
        case RejectReason::mapq: return "REJ_MQ";
        case RejectReason::msi_monomer: return "REJ_MSI1";
        case RejectReason::msi: return "REJ_MSI";
        case RejectReason::strand_bias: return "REJ_Bias";
    }
    std::unreachable();
}

std::vector<std::pair<std::string, std::string>> reject_filter_descriptions(const Config& conf) {
    using util::double_to_string;
    const std::string pre = "Rejected by the caller (-p): ";
    const auto row = [&](RejectReason r, const std::string& text) {
        return std::pair<std::string, std::string>{std::string(reject_id(r)), pre + text};
    };
    return {
        row(RejectReason::af, "allele frequency < " + double_to_string(conf.freq) + " (-f)"),
        row(RejectReason::min_reads, "high-quality variant reads < " + std::to_string(conf.minr) + " (-r)"),
        row(RejectReason::pmean, "mean position in reads < " + std::to_string(conf.read_pos_filter) + " (-P)"),
        row(RejectReason::qual, "mean base quality < " + double_to_string(conf.goodq) + " (-q)"),
        row(RejectReason::ref_mapq,
            "AF < 0.25 and the variant's mean mapping quality is much lower than the "
            "reference allele's"),
        row(RejectReason::splice, "deletion matching a splice junction"),
        row(RejectReason::qratio, "high- to low-quality read ratio < " + double_to_string(conf.qratio) + " (-o)"),
        row(RejectReason::mapq, "AF <= 0.3 and mean mapping quality < " + double_to_string(conf.mapq) + " (-O)"),
        row(RejectReason::msi_monomer, "monomer microsatellite (MSI >= 15) and AF <= " +
                                           double_to_string(conf.monomer_msi_frequency) + " (-mfreq)"),
        row(RejectReason::msi, "non-monomer microsatellite (MSI >= 12) and AF <= " +
                                   double_to_string(conf.non_monomer_msi_frequency) + " (-nmfreq)"),
        row(RejectReason::strand_bias, "AF < 0.2 and strand bias 2;1 (SNV or indel shorter than 3 bp)"),
    };
}

}  // namespace omnivar::core
