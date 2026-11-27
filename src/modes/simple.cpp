#include "modes/simple.hpp"

#include <initializer_list>
#include <optional>
#include <utility>
#include <vector>

#include "stats/fisher.hpp"
#include "util/checked_string.hpp"
#include "util/num_format.hpp"

namespace senri::modes {

namespace {

std::string join(std::initializer_list<std::string> parts) {
    std::string s;
    bool first = true;
    for (const auto& p : parts) {
        if (!first) {
            s += '\t';
        }
        s += p;
        first = false;
    }
    return s;
}

std::string i2s(int v) { return std::to_string(v); }

// `x == 0 ? 0 : new DecimalFormat(pattern).format(x)` (prints the Integer 0)
std::string df_or_zero(double v, int digits) { return v == 0 ? "0" : util::decimal_format(v, digits); }

}  // namespace

std::string simple_header(const Config& conf) {
    if (!conf.print_header) {
        return "";
    }
    std::string header =
        "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\tDepth\tAltDepth\tRefFwdReads\tRefRevReads\tAltFwdReads\t"
        "AltRevReads\tGenotype\tAF\tBias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\tHiAF\tExtraAF\tshift3\tMSI\t"
        "MSI_NT\tNM\tHiCnt\tHiCov\t5pFlankSeq\t3pFlankSeq\tSeg\tVarType\tDuprate\tSV_info";
    if (conf.crispr_cutting_site != 0) {
        header += "\tCRISPR";
    }
    return header + "\n";
}

std::string simple_output_line(const core::Variant* v, const Region& region, const std::string& sv_in, int position,
                               const Config& conf, const core::RunContext& run) {
    const core::Variant empty;
    const bool null_variant = v == nullptr;
    const core::Variant& x = null_variant ? empty : *v;

    const int start_position = null_variant ? position : x.start_position;
    const int end_position = null_variant ? position : x.end_position;
    const std::string ref_allele = null_variant ? "" : x.refallele;
    const std::string var_allele = null_variant ? "" : x.varallele;
    const std::string genotype = null_variant ? "" : x.genotype.value_or("0");
    const std::string bias = null_variant ? "0;0" : x.strand_bias_flag;
    const std::string left = null_variant ? "" : (x.leftseq.empty() ? "0" : x.leftseq);
    const std::string right = null_variant ? "" : (x.rightseq.empty() ? "0" : x.rightseq);
    const std::string var_type = null_variant ? "" : x.vartype;
    const std::string region_str = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
    const std::string sv = sv_in.empty() ? "0" : sv_in;

    std::string line;
    if (!conf.fisher) {
        // create_simple_variant_36columns
        line = join({run.sample,
                     region.gene,
                     region.chr,
                     i2s(start_position),
                     i2s(end_position),
                     ref_allele,
                     var_allele,
                     i2s(x.total_pos_coverage),
                     i2s(x.position_coverage),
                     i2s(x.ref_forward_coverage),
                     i2s(x.ref_reverse_coverage),
                     i2s(x.vars_count_on_forward),
                     i2s(x.vars_count_on_reverse),
                     genotype,
                     df_or_zero(x.frequency, 4),
                     bias,
                     df_or_zero(x.mean_position, 1),
                     i2s(x.is_at_least_at_2_positions ? 1 : 0),
                     df_or_zero(x.mean_quality, 1),
                     i2s(x.has_at_least_2_diff_qualities ? 1 : 0),
                     df_or_zero(x.mean_mapping_quality, 1),
                     df_or_zero(x.high_quality_to_low_quality_ratio, 3),
                     df_or_zero(x.high_quality_reads_frequency, 4),
                     df_or_zero(x.extra_frequency, 4),
                     i2s(x.shift3),
                     df_or_zero(x.msi, 3),
                     i2s(x.msint),
                     x.number_of_mismatches > 0 ? util::decimal_format(x.number_of_mismatches, 1) : "0",
                     i2s(x.hicnt),
                     i2s(x.hicov),
                     left,
                     right,
                     region_str,
                     var_type,
                     df_or_zero(x.duprate, 1),
                     sv});
    } else {
        // create_simple_variant_38columns
        const stats::FisherExact fisher = null_variant
                                              ? stats::FisherExact(0, 0, 0, 0)
                                              : stats::FisherExact(x.ref_forward_coverage, x.ref_reverse_coverage,
                                                                   x.vars_count_on_forward, x.vars_count_on_reverse);
        const double pvalue = fisher.p_value();
        const std::string oddratio = fisher.odd_ratio();
        const std::string hifreq_f =
            x.high_quality_reads_frequency == 0 ? "0" : util::decimal_format(x.high_quality_reads_frequency, 4);
        const double nm = x.number_of_mismatches > 0 ? x.number_of_mismatches : 0;
        const std::string nm_f = nm == 0 ? "0" : util::decimal_format(nm, 1);
        line = join({run.sample,
                     region.gene,
                     region.chr,
                     i2s(start_position),
                     i2s(end_position),
                     ref_allele,
                     var_allele,
                     i2s(x.total_pos_coverage),
                     i2s(x.position_coverage),
                     i2s(x.ref_forward_coverage),
                     i2s(x.ref_reverse_coverage),
                     i2s(x.vars_count_on_forward),
                     i2s(x.vars_count_on_reverse),
                     genotype,
                     util::rounded_value_to_print(x.frequency, 4),
                     bias,
                     util::rounded_value_to_print(x.mean_position, 1),
                     i2s(x.is_at_least_at_2_positions ? 1 : 0),
                     util::rounded_value_to_print(x.mean_quality, 1),
                     i2s(x.has_at_least_2_diff_qualities ? 1 : 0),
                     util::rounded_value_to_print(pvalue, 5),
                     oddratio,
                     util::rounded_value_to_print(x.mean_mapping_quality, 1),
                     util::rounded_value_to_print(x.high_quality_to_low_quality_ratio, 3),
                     hifreq_f,
                     util::rounded_value_to_print(x.extra_frequency, 4),
                     i2s(x.shift3),
                     util::rounded_value_to_print(x.msi, 3),
                     i2s(x.msint),
                     nm_f,
                     i2s(x.hicnt),
                     i2s(x.hicov),
                     left,
                     right,
                     region_str,
                     var_type,
                     util::rounded_value_to_print(x.duprate, 2),
                     sv});
    }
    if (conf.crispr_cutting_site != 0) {
        line += "\t" + i2s(x.crispr);
    }
    if (conf.debug) {
        line += "\t" + x.debug;
    }
    return line + "\n";
}

void simple_post_process(core::AlignedVars& aligned, const core::RegionScope& scope, std::string& out) {
    const Config& conf = scope.conf;
    int last_position = 0;
    aligned.for_each([&](const int& position, core::Vars& variants_on_position) {
        try {
            last_position = position;
            // (variant, reject reason): rejected variants keep their place among the calls.
            std::vector<std::pair<core::VariantPtr, std::optional<core::RejectReason>>> vrefs;
            if (variants_on_position.sv.empty()) {
                if (position < scope.region.start || position > scope.region.end) {
                    return;
                }
            }
            if (variants_on_position.variants.empty()) {
                return;
            }
            for (const core::VariantPtr& vref : variants_on_position.variants) {
                if (vref->refallele.contains('N')) {
                    continue;
                }
                if (vref->refallele == vref->varallele) {
                    continue;
                }
                vref->vartype = vref->var_type();
                auto reason = vref->reject_reason(variants_on_position.reference_variant.get(), vref->vartype,
                                                  scope.splice, conf);
                if (reason) {
                    if (!conf.do_pileup || *reason == core::RejectReason::no_ref) {
                        continue;
                    }
                }
                vrefs.emplace_back(vref, reason);
            }
            for (const auto& [vref, reason] : vrefs) {
                if (reason) {
                    // -p: printed from a copy, so a rejected variant is never changed (upstream never
                    // trims rejected alleles; with an empty ALT adjComplex would throw).
                    core::Variant c = *vref;
                    if (c.vartype == "Complex" && !c.varallele.empty()) {
                        c.adj_complex();
                    }
                    if (conf.crispr_cutting_site == 0) {
                        c.crispr = 0;
                    }
                    out += kRejectedLineMarker;
                    out += core::reject_id(*reason);
                    out += '\t';
                    out += simple_output_line(&c, scope.region, variants_on_position.sv, position, conf, scope.run);
                    continue;
                }
                if (vref->vartype == "Complex") {
                    vref->adj_complex();
                }
                if (conf.crispr_cutting_site == 0) {
                    vref->crispr = 0;
                }
                out += simple_output_line(vref.get(), scope.region, variants_on_position.sv, position, conf, scope.run);
            }
        } catch (const util::RecoverableError& e) {
            core::print_exception_and_continue(e, "position", std::to_string(last_position), &scope.region,
                                               scope.exceptions);
        }
    });
}

}  // namespace senri::modes
