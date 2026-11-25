#include "modes/somatic.hpp"

#include <algorithm>
#include <cstdio>
#include <initializer_list>
#include <optional>
#include <print>
#include <set>
#include <utility>
#include <vector>

#include "core/variation_utils.hpp"
#include "io/reference.hpp"
#include "modes/simple.hpp"
#include "stats/fisher.hpp"
#include "util/checked_string.hpp"
#include "util/num_format.hpp"

namespace omnivar::modes {

namespace {

using core::Variant;
using core::VariantPtr;
using core::Vars;

const std::string STRONG_SOMATIC = "StrongSomatic";
const std::string SAMPLE_SPECIFIC = "SampleSpecific";
const std::string DELETION = "Deletion";
const std::string LIKELY_LOH = "LikelyLOH";
const std::string GERMLINE = "Germline";
const std::string STRONG_LOH = "StrongLOH";
const std::string LIKELY_SOMATIC = "LikelySomatic";
const std::string AF_DIFF = "AFDiff";
const std::string FALSE_TYPE = "FALSE";
const std::string REJECTED = "Rejected";  // -p: STATUS of rejected records (decision 2026-11-02)

// -p FILTER ID of a reason in one sample: REJ_<reason>_T or _N.
std::string reject_id(core::RejectReason reason, bool tumour) {
    return std::string(core::reject_id(reason)) + (tumour ? "_T" : "_N");
}

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
// `x == 0 ? 0 : new DecimalFormat(pattern).format(x)`
std::string df_or_zero(double v, int digits) { return v == 0 ? "0" : util::decimal_format(v, digits); }
std::string rv(double v, int digits) { return util::rounded_value_to_print(v, digits); }

// Patterns.MINUS_NUM_NUM "-\d\d" find
bool minus_num_num(const std::string& s) {
    for (std::size_t i = 0; i + 2 < s.size(); ++i) {
        if (s[i] == '-' && std::isdigit(static_cast<unsigned char>(s[i + 1])) &&
            std::isdigit(static_cast<unsigned char>(s[i + 2]))) {
            return true;
        }
    }
    return false;
}

// VariationUtils.getVarMaybe(vars, varn, key)
Variant* get_varn(Vars& vars, const std::string& key) {
    VariantPtr* p = vars.var_description_string_to_variants.get(key);
    return p ? p->get() : nullptr;
}
Variant* get_var(Vars& vars, std::size_t i) { return i < vars.variants.size() ? vars.variants[i].get() : nullptr; }

// One sample's columns of SomaticOutputVariant.
struct SampleCols {
    int total_coverage = 0;
    int variant_coverage = 0;
    int ref_forward = 0;
    int ref_reverse = 0;
    int var_forward = 0;
    int var_reverse = 0;
    std::string genotype = "0";
    double frequency = 0;
    std::string strand_bias_flag = "0";
    double mean_position = 0;
    int at_least_2_positions = 0;
    double mean_quality = 0;
    int at_least_2_qualities = 0;
    double mean_mapping_quality = 0;
    double hi_lo_ratio = 0;
    double hi_freq = 0;
    double extra_frequency = 0;
    double nm = 0;
    double duprate = 0;

    static SampleCols of(const Variant* v) {
        SampleCols c;
        if (v == nullptr) {
            return c;
        }
        c.total_coverage = v->total_pos_coverage;
        c.variant_coverage = v->position_coverage;
        c.ref_forward = v->ref_forward_coverage;
        c.ref_reverse = v->ref_reverse_coverage;
        c.var_forward = v->vars_count_on_forward;
        c.var_reverse = v->vars_count_on_reverse;
        c.genotype = v->genotype.value_or("0");
        c.frequency = v->frequency;
        c.strand_bias_flag = v->strand_bias_flag;
        c.mean_position = v->mean_position;
        c.at_least_2_positions = v->is_at_least_at_2_positions ? 1 : 0;
        c.mean_quality = v->mean_quality;
        c.at_least_2_qualities = v->has_at_least_2_diff_qualities ? 1 : 0;
        c.mean_mapping_quality = v->mean_mapping_quality;
        c.hi_lo_ratio = v->high_quality_to_low_quality_ratio;
        c.hi_freq = v->high_quality_reads_frequency;
        c.extra_frequency = v->extra_frequency;
        c.nm = v->number_of_mismatches;
        c.duprate = v->duprate;
        return c;
    }

    // create_somatic_variant_55columns, one sample
    [[nodiscard]] std::string columns55() const {
        return join({i2s(total_coverage), i2s(variant_coverage), i2s(ref_forward), i2s(ref_reverse), i2s(var_forward),
                     i2s(var_reverse), genotype, df_or_zero(frequency, 4), strand_bias_flag,
                     df_or_zero(mean_position, 1), i2s(at_least_2_positions), df_or_zero(mean_quality, 1),
                     i2s(at_least_2_qualities), df_or_zero(mean_mapping_quality, 1), df_or_zero(hi_lo_ratio, 3),
                     df_or_zero(hi_freq, 4), df_or_zero(extra_frequency, 4),
                     nm > 0 ? util::decimal_format(nm, 1) : "0"});
    }
    // create_somatic_variant_61columns, one sample (nm already clamped at 0)
    [[nodiscard]] std::string columns61(double pvalue, const std::string& oddratio) const {
        return join({i2s(total_coverage),
                     i2s(variant_coverage),
                     i2s(ref_forward),
                     i2s(ref_reverse),
                     i2s(var_forward),
                     i2s(var_reverse),
                     genotype,
                     rv(frequency, 4),
                     strand_bias_flag,
                     rv(mean_position, 1),
                     i2s(at_least_2_positions),
                     rv(mean_quality, 1),
                     i2s(at_least_2_qualities),
                     rv(mean_mapping_quality, 1),
                     rv(hi_lo_ratio, 3),
                     rv(hi_freq, 4),
                     rv(extra_frequency, 4),
                     rv(nm > 0 ? nm : 0, 1),
                     rv(pvalue, 5),
                     oddratio});
    }
};

std::string somatic_line(const Variant* begin, const Variant* end, const Variant* tumor, const Variant* normal,
                         const Region& region, const std::string& sv1, const std::string& sv2,
                         const std::string& var_label, const Config& conf, const core::RunContext& run) {
    int start_position = 0;
    int end_position = 0;
    std::string ref_allele;
    std::string var_allele;
    std::string var_type;
    std::string debug;
    if (begin != nullptr) {
        start_position = begin->start_position;
        end_position = begin->end_position;
        ref_allele = begin->refallele;
        var_allele = begin->varallele;
        var_type = begin->vartype;
        debug = begin->debug;
    }
    int shift3 = 0;
    double msi = 0;
    int msint = 0;
    std::string left;
    std::string right;
    if (end != nullptr) {
        shift3 = end->shift3;
        msi = end->msi;
        msint = end->msint;
        left = end->leftseq.empty() ? "0" : end->leftseq;
        right = end->rightseq.empty() ? "0" : end->rightseq;
    }
    const SampleCols t = SampleCols::of(tumor);
    const SampleCols n = SampleCols::of(normal);
    const std::string region_str = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
    const std::string var1sv = sv1.empty() ? "0" : sv1;
    const std::string var2sv = sv2.empty() ? "0" : sv2;
    const std::string head =
        join({run.sample, region.gene, region.chr, i2s(start_position), i2s(end_position), ref_allele, var_allele});
    std::string line;
    if (!conf.fisher) {
        line =
            join({head, t.columns55(), n.columns55(), i2s(shift3), df_or_zero(msi, 3), i2s(msint), left, right,
                  region_str, var_label, var_type, df_or_zero(t.duprate, 1), var1sv, df_or_zero(n.duprate, 1), var2sv});
    } else {
        // calculateFisherSomatic
        const auto strand = [](const Variant* v) {
            return v ? stats::FisherExact(v->ref_forward_coverage, v->ref_reverse_coverage, v->vars_count_on_forward,
                                          v->vars_count_on_reverse)
                     : stats::FisherExact(0, 0, 0, 0);
        };
        const stats::FisherExact f1 = strand(tumor);
        const stats::FisherExact f2 = strand(normal);
        const int tref = std::max(0, t.total_coverage - t.variant_coverage);
        const int rref = std::max(0, n.total_coverage - n.variant_coverage);
        const stats::FisherExact f(t.variant_coverage, tref, n.variant_coverage, rref);
        const double pvalue_greater = f.p_value_greater();
        const double pvalue_less = f.p_value_less();
        const double pvalue = pvalue_less < pvalue_greater ? pvalue_less : pvalue_greater;
        const std::string msi_f = msi == 0 ? "0" : util::decimal_format(msi, 3);
        line = join({head, t.columns61(f1.p_value(), f1.odd_ratio()), n.columns61(f2.p_value(), f2.odd_ratio()),
                     i2s(shift3), msi_f, i2s(msint), left, right, region_str, var_label, var_type, rv(t.duprate, 1),
                     var1sv, rv(n.duprate, 1), var2sv, rv(pvalue, 5), f.odd_ratio()});
    }
    if (conf.debug) {
        line += "\t" + debug;
    }
    return line + "\n";
}

struct CombineAnalysisData {
    int max_read_length;
    std::string type;
};

class SomaticPostProcess {
public:
    SomaticPostProcess(const PipelineEnv& env, const Region& region, std::unordered_set<std::string>& splice,
                       std::string& out)
        : env_(env), conf_(env.conf), region_(region), splice_(splice), out_(out) {}

    void accept(PipelineResult& bam1, PipelineResult& bam2);

private:
    void calling_for_one_sample(Vars& variants, bool is_first_cover, const std::string& var_label);
    void calling_for_both_samples(int position, Vars& v1, Vars& v2);
    void print_variations_from_first_sample(int position, Vars& v1, Vars& v2);
    void print_variations_from_second_sample(int position, Vars& v1, Vars& v2);
    std::string determinate_type(Vars& variants, const Variant& standard, Variant& to_compare);
    CombineAnalysisData combine_analysis(const Variant& variant1, Variant& variant2, int position,
                                         const std::string& description_string, int max_read_length);
    void print(const Variant* b, const Variant* e, const Variant* t, const Variant* n, const std::string& sv1,
               const std::string& sv2, const std::string& label) {
        out_ += somatic_line(b, e, t, n, region_, sv1, sv2, label, conf_, env_.run);
    }
    bool good(const Variant& v, const Vars& vars, const std::string& type) const {
        return v.is_good_var(vars.reference_variant.get(), type, splice_, conf_);
    }
    std::optional<core::RejectReason> reason(const Variant& v, const Vars& vars, const std::string& type) const {
        return v.reject_reason(vars.reference_variant.get(), type, splice_, conf_);
    }
    // -p: note that the allele `description` was rejected at this position for `id`.
    void reject(const std::string& description, std::string id) {
        if (conf_.do_pileup) {
            rejected_.emplace_back(description, std::move(id));
        }
    }
    void reject(const Variant& v, std::optional<core::RejectReason> r, bool tumour) {
        if (r && *r != core::RejectReason::no_ref) {
            reject(v.description_string, reject_id(*r, tumour));
        }
    }
    void flush_rejected(Vars* v1, Vars* v2);

    // -p: rejected alleles of the current position, in order (allele description, FILTER ID).
    std::vector<std::pair<std::string, std::string>> rejected_;

    const PipelineEnv& env_;
    const Config& conf_;
    const Region& region_;
    std::unordered_set<std::string>& splice_;
    std::string& out_;
    int max_read_length_ = 0;
};

void SomaticPostProcess::accept(PipelineResult& bam1, PipelineResult& bam2) {
    max_read_length_ = std::max(bam1.max_read_length, bam2.max_read_length);
    std::set<int> positions;  // HashSet of both key sets, then sorted
    for (int k : bam1.aligned.keys()) {
        positions.insert(k);
    }
    for (int k : bam2.aligned.keys()) {
        positions.insert(k);
    }
    int last_position = 0;
    for (int position : positions) {
        try {
            last_position = position;
            rejected_.clear();
            if (position < region_.start || position > region_.end) {
                continue;
            }
            Vars* v1 = bam1.aligned.get(position);
            Vars* v2 = bam2.aligned.get(position);
            if (v1 == nullptr && v2 == nullptr) {  // both samples have no coverage
                continue;
            }
            if (v1 == nullptr) {  // no coverage for sample 1
                calling_for_one_sample(*v2, true, DELETION);
            } else if (v2 == nullptr) {  // no coverage for sample 2
                calling_for_one_sample(*v1, false, SAMPLE_SPECIFIC);
            } else {  // both samples have coverage
                calling_for_both_samples(position, *v1, *v2);
            }
            flush_rejected(v1, v2);
        } catch (const util::RecoverableError& e) {
            core::print_exception_and_continue(e, "position", std::to_string(last_position), &region_, env_.exceptions);
        }
    }
}

void SomaticPostProcess::calling_for_one_sample(Vars& variants, bool is_first_cover, const std::string& var_label) {
    if (variants.variants.empty()) {
        return;
    }
    for (const VariantPtr& vp : variants.variants) {
        Variant& variant = *vp;
        if (variant.refallele == variant.varallele) {
            continue;
        }
        variant.vartype = variant.var_type();
        if (const auto r = reason(variant, variants, variant.vartype)) {
            reject(variant, r, !is_first_cover);  // -p
            continue;
        }
        if (variant.vartype == "Complex") {
            variant.adj_complex();
        }
        if (is_first_cover) {
            print(&variant, &variant, nullptr, &variant, "", variants.sv, var_label);
        } else {
            print(&variant, &variant, &variant, nullptr, variants.sv, "", var_label);
        }
    }
}

void SomaticPostProcess::calling_for_both_samples(int position, Vars& v1, Vars& v2) {
    if (v1.variants.empty() && v2.variants.empty()) {
        return;
    }
    if (!v1.variants.empty()) {
        print_variations_from_first_sample(position, v1, v2);
    } else if (!v2.variants.empty()) {  // sample 1 has only reference
        print_variations_from_second_sample(position, v1, v2);
    }
}

void SomaticPostProcess::print_variations_from_first_sample(int position, Vars& v1, Vars& v2) {
    std::size_t processed = 0;
    while (processed < v1.variants.size() && good(*v1.variants[processed], v1, v1.variants[processed]->var_type())) {
        // The Variant object stays alive even if combineAnalysis replaces map entries.
        const VariantPtr vref_ptr = v1.variants[processed];
        Variant& vref = *vref_ptr;
        if (vref.refallele == vref.varallele) {
            processed++;
            continue;
        }
        const std::string nt = vref.description_string;
        vref.vartype = vref.var_type();
        if (vref.vartype == "Complex") {
            vref.adj_complex();
        }
        Variant* v2nt = get_varn(v2, nt);
        if (v2nt != nullptr) {
            const std::string type = determinate_type(v2, vref, *v2nt);
            print(&vref, v2nt, &vref, v2nt, v1.sv, v2.sv, type);
        } else {  // sample 1 only, should be strong somatic
            Variant empty_for_print;
            const Variant* var_for_print = &empty_for_print;
            if (!v2.variants.empty()) {
                const Variant* v2r = get_var(v2, 0);
                empty_for_print.total_pos_coverage = v2r != nullptr ? v2r->total_pos_coverage : 0;
                empty_for_print.ref_forward_coverage = v2r != nullptr ? v2r->ref_forward_coverage : 0;
                empty_for_print.ref_reverse_coverage = v2r != nullptr ? v2r->ref_reverse_coverage : 0;
            } else if (v2.reference_variant) {
                var_for_print = v2.reference_variant.get();
            } else {
                var_for_print = nullptr;
            }
            std::string type = STRONG_SOMATIC;
            VariantPtr new_v2nt;
            if (vref.vartype != "SNV" && (nt.size() > 10 || minus_num_num(nt))) {
                new_v2nt = std::make_shared<Variant>();
                // Ensure it's initialized before passing to combineAnalysis
                v2.var_description_string_to_variants.put(nt, new_v2nt);
                if (vref.position_coverage < conf_.minr + 3 && nt.find('<') == std::string::npos) {
                    const CombineAnalysisData tpl = combine_analysis(vref, *new_v2nt, position, nt, max_read_length_);
                    max_read_length_ = tpl.max_read_length;
                    if (tpl.type == FALSE_TYPE) {
                        reject(nt, "REJ_Combine_T");  // -p
                        processed++;
                        continue;
                    }
                    if (!tpl.type.empty()) {
                        type = tpl.type;
                    }
                }
            }
            if (type == STRONG_SOMATIC) {
                print(&vref, &vref, &vref, var_for_print, v1.sv, v2.sv, STRONG_SOMATIC);
            } else {
                print(&vref, &vref, &vref, new_v2nt.get(), v1.sv, v2.sv, type);
            }
        }
        processed++;
    }
    // -p: the loop stops at the first rejected tumour allele; the alleles after it are dropped
    // too (masked), whether or not isGoodVar would reject them. Reported after the LOH branch
    // below, which can still print one of them (Germline / LikelyLOH via the normal).
    std::vector<std::pair<std::string, std::string>> tumour_rejected;
    if (conf_.do_pileup) {
        for (std::size_t k = processed; k < v1.variants.size(); ++k) {
            const Variant& v = *v1.variants[k];
            if (v.refallele == v.varallele) {
                continue;
            }
            const auto r = reason(v, v1, v.var_type());
            if (r == core::RejectReason::no_ref) {
                continue;
            }
            tumour_rejected.emplace_back(v.description_string, r ? reject_id(*r, true) : "REJ_Masked_T");
        }
    }
    std::set<std::string> printed_by_loh;
    const auto report_tumour_rejected = [&] {
        for (auto& [description, id] : tumour_rejected) {
            if (!printed_by_loh.contains(description)) {
                reject(description, std::move(id));
            }
        }
    };
    if (processed == 0) {
        if (v2.variants.empty()) {
            report_tumour_rejected();
            return;
        }
        for (const VariantPtr& v2var_ptr : v2.variants) {
            Variant& v2var = *v2var_ptr;
            v2var.vartype = v2var.var_type();
            if (const auto r = reason(v2var, v2, v2var.vartype)) {
                reject(v2var, r, false);  // -p
                continue;
            }
            // potential LOH
            const std::string nt = v2var.description_string;
            Variant* v1nt = get_varn(v1, nt);
            if (v1nt != nullptr) {
                if (v1nt->refallele == v1nt->varallele) {
                    continue;
                }
                const std::string type = v1nt->frequency < conf_.lofreq ? LIKELY_LOH : GERMLINE;
                if (v2var.vartype == "Complex") {
                    v1nt->adj_complex();
                }
                v1nt->vartype = v1nt->var_type();
                print(v1nt, &v2var, v1nt, &v2var, v1.sv, v2.sv, type);
                printed_by_loh.insert(nt);
            } else {
                if (v2var.refallele == v2var.varallele) {
                    continue;
                }
                const Variant* v1var = get_var(v1, 0);
                const int tcov = v1var != nullptr ? v1var->total_pos_coverage : 0;
                const Variant* v1ref = v1.reference_variant.get();
                const int fwd = v1ref != nullptr ? v1ref->vars_count_on_forward : 0;
                const int rev = v1ref != nullptr ? v1ref->vars_count_on_reverse : 0;
                std::optional<std::string> genotype =
                    v1var != nullptr
                        ? v1var->genotype
                        : std::optional<std::string>(
                              v1ref != nullptr ? v1ref->description_string + "/" + v1ref->description_string : "N/N");
                if (v2var.vartype == "Complex") {
                    v2var.adj_complex();
                }
                Variant var_for_print;
                var_for_print.total_pos_coverage = tcov;
                var_for_print.ref_forward_coverage = fwd;
                var_for_print.ref_reverse_coverage = rev;
                var_for_print.genotype = genotype;
                print(&v2var, &v2var, &var_for_print, &v2var, "", v2.sv, STRONG_LOH);
            }
        }
    }
    report_tumour_rejected();
}

void SomaticPostProcess::print_variations_from_second_sample(int position, Vars& v1, Vars& v2) {
    for (const VariantPtr& v2var_ptr : v2.variants) {
        Variant& v2var = *v2var_ptr;
        if (v2var.refallele == v2var.varallele) {
            continue;
        }
        v2var.vartype = v2var.var_type();
        if (const auto r = reason(v2var, v2, v2var.vartype)) {
            reject(v2var, r, false);  // -p
            continue;
        }
        // potential LOH
        const std::string description_string = v2var.description_string;
        std::string type = STRONG_LOH;
        VariantPtr& v1nt_slot = v1.var_description_string_to_variants.get_or_create(description_string);
        if (!v1nt_slot) {
            v1nt_slot = std::make_shared<Variant>();  // computeIfAbsent
        }
        const VariantPtr v1nt = v1nt_slot;
        v1nt->position_coverage = 0;
        std::string new_type;
        VariantPtr* v2d = v2.var_description_string_to_variants.get(description_string);
        if (v2d == nullptr || !*v2d) {
            util::throw_npe("varDescriptionStringToVariants.get(" + description_string + ")");
        }
        const VariantPtr v2desc = *v2d;
        if (v2desc->position_coverage < conf_.minr + 3 && description_string.find('<') == std::string::npos &&
            (description_string.size() > 10 || minus_num_num(description_string))) {
            const CombineAnalysisData tpl =
                combine_analysis(*v2desc, *v1nt, position, description_string, max_read_length_);
            max_read_length_ = tpl.max_read_length;
            new_type = tpl.type;
            if (new_type == FALSE_TYPE) {
                reject(description_string, "REJ_Combine_N");  // -p
                continue;
            }
        }
        const Variant* var_for_print = nullptr;
        if (!new_type.empty()) {
            type = new_type;
            var_for_print = v1nt.get();
        } else {
            var_for_print = v1.reference_variant.get();
        }
        if (v2var.vartype == "Complex") {
            v2var.adj_complex();
        }
        print(&v2var, &v2var, var_for_print, &v2var, "", v2.sv, type);
    }
}

// -p: one record per rejected allele of the position, FILTER = its reasons (both samples when
// both rejected it). Built read-only, so the calls are unchanged: each sample's columns are that
// sample's variant for the allele (not an empty placeholder), or its reference variant; the alleles are trimmed
// (adjComplex) on a copy; STATUS is "Rejected" (the upstream labelling has side effects and is not run).
void SomaticPostProcess::flush_rejected(Vars* v1, Vars* v2) {
    if (rejected_.empty()) {
        return;
    }
    std::vector<std::pair<std::string, std::string>> merged;  // first-seen order
    for (const auto& [description, id] : rejected_) {
        auto it = std::find_if(merged.begin(), merged.end(), [&](const auto& m) { return m.first == description; });
        if (it == merged.end()) {
            merged.emplace_back(description, id);
        } else if (!(";" + it->second + ";").contains(";" + id + ";")) {
            it->second += ";" + id;
        }
    }
    rejected_.clear();
    for (auto& [description, ids] : merged) {  // tumour reasons first
        std::vector<std::string> parts;
        for (std::size_t start = 0; start <= ids.size();) {
            const auto semi = ids.find(';', start);
            parts.push_back(ids.substr(start, semi == std::string::npos ? std::string::npos : semi - start));
            start = semi == std::string::npos ? ids.size() + 1 : semi + 1;
        }
        std::stable_partition(parts.begin(), parts.end(), [](const std::string& id) { return id.ends_with("_T"); });
        ids.clear();
        for (const std::string& id : parts) {
            ids += (ids.empty() ? "" : ";") + id;
        }
    }
    const auto allele = [](Vars* vars, const std::string& description) -> const Variant* {
        if (vars == nullptr) {
            return nullptr;
        }
        const VariantPtr* p = std::as_const(vars->var_description_string_to_variants).get(description);
        // Not the empty placeholders the somatic module puts in the map before combineAnalysis
        // (computeIfAbsent / `new Variant()`): they have no alleles.
        return p != nullptr && *p && !(*p)->refallele.empty() ? p->get() : nullptr;
    };
    const auto reference = [](Vars* vars) -> const Variant* {
        return vars != nullptr ? vars->reference_variant.get() : nullptr;
    };
    for (const auto& [description, ids] : merged) {
        const Variant* t = allele(v1, description);
        const Variant* n = allele(v2, description);
        const Variant* base = t != nullptr ? t : n;
        if (base == nullptr) {
            continue;
        }
        Variant c = *base;
        if (c.vartype.empty()) {
            c.vartype = c.var_type();
        }
        // idempotent: alleles already trimmed upstream stay as they are. Not with an empty ALT
        // (adjComplex would throw on charAt(0); upstream never trims rejected alleles).
        if (c.vartype == "Complex" && !c.varallele.empty()) {
            c.adj_complex();
        }
        const Variant* tumour = t != nullptr ? &c : reference(v1);
        const Variant* normal = t != nullptr ? (n != nullptr ? n : reference(v2)) : &c;
        out_ += kRejectedLineMarker;
        out_ += ids;
        out_ += '\t';
        out_ += somatic_line(&c, &c, tumour, normal, region_, v1 != nullptr ? v1->sv : "", v2 != nullptr ? v2->sv : "",
                             REJECTED, conf_, env_.run);
    }
}

std::string SomaticPostProcess::determinate_type(Vars& variants, const Variant& standard, Variant& to_compare) {
    std::string type;
    if (good(to_compare, variants, standard.vartype)) {
        if (standard.frequency > (1 - conf_.lofreq) && to_compare.frequency < 0.8 && to_compare.frequency > 0.2) {
            type = LIKELY_LOH;
        } else if (to_compare.frequency < conf_.lofreq || to_compare.position_coverage <= 1) {
            type = LIKELY_SOMATIC;
        } else {
            type = GERMLINE;
        }
    } else {
        if (to_compare.frequency < conf_.lofreq || to_compare.position_coverage <= 1) {
            type = LIKELY_SOMATIC;
        } else {
            type = AF_DIFF;
        }
    }
    if (to_compare.is_noise(conf_) && standard.vartype == "SNV") {
        type = STRONG_SOMATIC;
    }
    return type;
}

CombineAnalysisData SomaticPostProcess::combine_analysis(const Variant& variant1, Variant& variant2, int position,
                                                         const std::string& description_string, int max_read_length) {
    if (conf_.y) {
        std::println(stderr, "Start Combine {} {}", position, description_string);
    }
    // Don't do it for structural variants
    if (variant1.end_position - variant1.start_position > conf_.sv_min_len) {
        return {max_read_length, ""};
    }
    const Region region{region_.chr, variant1.start_position - max_read_length, variant1.end_position + max_read_length,
                        ""};
    auto ref = io::get_reference(env_.provider, region, conf_, env_.run.chr_lengths);
    if (!ref) {
        throw util::RecoverableError(ref.error().message);
    }
    PipelineResult tpl =
        run_pipeline(env_, region, conf_.bam.bam1() + ":" + conf_.bam.bam2(), max_read_length, *ref, splice_);
    max_read_length = tpl.max_read_length;
    Vars* vars = tpl.aligned.get(position);
    Variant* vref = vars ? get_varn(*vars, description_string) : nullptr;
    if (vref == nullptr) {
        return {max_read_length, FALSE_TYPE};
    }
    if (conf_.y) {
        std::println(stderr, "Combine: 1: {} comb: {}", variant1.position_coverage, vref->position_coverage);
    }
    if (vref->position_coverage - variant1.position_coverage >= conf_.minr) {
        variant2.total_pos_coverage = std::max(0, vref->total_pos_coverage - variant1.total_pos_coverage);
        variant2.position_coverage = std::max(0, vref->position_coverage - variant1.position_coverage);
        variant2.ref_forward_coverage = std::max(0, vref->ref_forward_coverage - variant1.ref_forward_coverage);
        variant2.ref_reverse_coverage = std::max(0, vref->ref_reverse_coverage - variant1.ref_reverse_coverage);
        variant2.vars_count_on_forward = std::max(0, vref->vars_count_on_forward - variant1.vars_count_on_forward);
        variant2.vars_count_on_reverse = std::max(0, vref->vars_count_on_reverse - variant1.vars_count_on_reverse);
        if (variant2.position_coverage != 0) {
            const double pc2 = variant2.position_coverage;
            const auto back = [&](double comb, double one) {
                return (comb * vref->position_coverage - one * variant1.position_coverage) / pc2;
            };
            variant2.mean_position = back(vref->mean_position, variant1.mean_position);
            variant2.mean_quality = back(vref->mean_quality, variant1.mean_quality);
            variant2.mean_mapping_quality = back(vref->mean_mapping_quality, variant1.mean_mapping_quality);
            variant2.high_quality_reads_frequency =
                back(vref->high_quality_reads_frequency, variant1.high_quality_reads_frequency);
            variant2.extra_frequency = back(vref->extra_frequency, variant1.extra_frequency);
            variant2.number_of_mismatches = back(vref->number_of_mismatches, variant1.number_of_mismatches);
        } else {
            variant2.mean_position = 0;
            variant2.mean_quality = 0;
            variant2.mean_mapping_quality = 0;
            variant2.high_quality_reads_frequency = 0;
            variant2.extra_frequency = 0;
            variant2.number_of_mismatches = 0;
        }
        variant2.is_at_least_at_2_positions = true;
        variant2.has_at_least_2_diff_qualities = true;
        if (variant2.total_pos_coverage <= 0) {
            return {max_read_length, FALSE_TYPE};
        }
        variant2.frequency = variant2.position_coverage / static_cast<double>(variant2.total_pos_coverage);
        // Can't back calculate and should be inaccurate
        variant2.high_quality_to_low_quality_ratio = variant1.high_quality_to_low_quality_ratio;
        variant2.genotype = vref->genotype;
        variant2.strand_bias_flag =
            std::to_string(core::strand_bias(variant2.ref_forward_coverage, variant2.ref_reverse_coverage, conf_)) +
            ";" +
            std::to_string(core::strand_bias(variant2.vars_count_on_forward, variant2.vars_count_on_reverse, conf_));
        return {max_read_length, GERMLINE};
    }
    if (vref->position_coverage < variant1.position_coverage - 2) {
        if (conf_.y) {
            std::println(stderr, "Combine produce less: {} {} {} {} {}", region_.chr, position, description_string,
                         vref->position_coverage, variant1.position_coverage);
        }
        return {max_read_length, FALSE_TYPE};
    }
    return {max_read_length, ""};
}

}  // namespace

std::string somatic_header(const Config& conf) {
    if (!conf.print_header) {
        return "";
    }
    return "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\t"
           "Depth\tAltDepth\tRefFwdReads\tRefRevReads\tAltFwdReads\tAltRevReads\tGenotype\tAF\t"
           "Bias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\tHiAF\tExtraAF\tNM\t"
           "Depth\tAltDepth\tRefFwdReads\tRefRevReads\tAltFwdReads\tAltRevReads\tGenotype\tAF\t"
           "Bias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\tHiAF\tExtraAF\tNM\t"
           "shift3\tMSI\tMSI_NT\t5pFlankSeq\t3pFlankSeq\tSeg\tVarLabel\tVarType\t"
           "Duprate1\tSV_info1\tDuprate2\tSV_info2\n";
}

void process_somatic_region(const PipelineEnv& env, const Region& region, std::string& out) {
    // AbstractMode.tryToGetReference: a failure stops the run
    auto reference = io::get_reference(env.provider, region, env.conf, env.run.chr_lengths);
    if (!reference) {
        throw core::FatalError(reference.error().message);
    }
    std::unordered_set<std::string> splice;
    PipelineResult bam1 = run_pipeline(env, region, env.conf.bam.bam1(), 0, *reference, splice);
    PipelineResult bam2 = run_pipeline(env, region, env.conf.bam.bam2(), bam1.max_read_length, *reference, splice);
    SomaticPostProcess(env, region, splice, out).accept(bam1, bam2);
}

}  // namespace omnivar::modes
