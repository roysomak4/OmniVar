#include "core/vars_builder.hpp"

#include <algorithm>
#include <cstdio>
#include <print>
#include <set>

#include "core/patterns.hpp"
#include "core/variation_utils.hpp"
#include "util/checked_string.hpp"
#include "util/num_format.hpp"
#include "util/number_parse.hpp"

namespace senri::core {

namespace {

constexpr int kRef20Bases = 20;
constexpr int kRef30Bases = 30;
constexpr int kRef50Bases = 50;
constexpr int kRef70Bases = 70;

using util::round_half_even;

// "((unit)+)$": the longest suffix of s made of whole copies of unit ("" if none).
std::string repeat_suffix(const std::string& s, const std::string& unit) {
    std::size_t len = 0;
    while (len + unit.size() <= s.size() && s.compare(s.size() - len - unit.size(), unit.size(), unit) == 0) {
        len += unit.size();
    }
    return s.substr(s.size() - len);
}

// "^((unit)+)": the longest prefix of s made of whole copies of unit.
std::size_t repeat_prefix_length(const std::string& s, const std::string& unit) {
    std::size_t len = 0;
    while (len + unit.size() <= s.size() && s.compare(len, unit.size(), unit) == 0) {
        len += unit.size();
    }
    return len;
}

int to_int(const std::string& s) {
    auto r = util::parse_int(s);
    if (!r) {
        throw util::RecoverableError("NumberFormatException: " + r.error());
    }
    return *r;
}

// ref.get(pos).toString() if present, else "" (the `ref.containsKey(p) ? ... : ""` idiom)
std::string ref_or_empty(const RefMap& ref, int pos) {
    auto c = ref_get(ref, pos);
    return c ? std::string(1, *c) : std::string();
}

// String.valueOf(Character) / string concatenation of a possibly-null Character
std::string ref_or_null(const RefMap& ref, int pos) {
    auto c = ref_get(ref, pos);
    return c ? std::string(1, *c) : std::string("null");
}

// ToVarsBuilder.validateRefallele: IUPAC ambiguity codes -> first base alphabetically.
std::string validate_refallele(std::string refallele) {
    auto code = [](char c) -> char {
        switch (c) {
            case 'M':
            case 'R':
            case 'W':
            case 'V':
            case 'H':
            case 'D': return 'A';
            case 'S':
            case 'Y':
            case 'B': return 'C';
            case 'K': return 'G';
            default: return 0;
        }
    };
    for (std::size_t i = 0; i < refallele.size(); ++i) {
        const char c = code(refallele[i]);
        if (c != 0) {
            // replaceFirst(refBase, code): the first occurrence of this letter
            const auto p = refallele.find(refallele[i]);
            refallele[p] = c;
        }
    }
    return refallele;
}

class Builder {
public:
    Builder(VariationData& data, const RegionScope& scope)
        : data_(data), scope_(scope), conf_(scope.conf), ref_(scope.reference.reference_sequences) {}

    AlignedVars run();

private:
    int chr_length() const {
        // getOrElse(instance().chrLengths, region.chr, 0) (also inserts the 0 upstream; unused)
        auto it = scope_.run.chr_lengths.find(scope_.region.chr);
        return it != scope_.run.chr_lengths.end() ? it->second : 0;
    }
    bool is_the_same_variation_on_ref(int position, const VariationMap& vars);
    int calc_hicov(const VariationMap& non_insertion);
    void create_variant(double duprate, AlignedVars& aligned, int position, const VariationMap& vars,
                        int total_pos_coverage, std::vector<VariantPtr>& var, std::vector<std::string>& debug_lines,
                        const std::vector<std::string>& keys, int hicov);
    int create_insertion(double duprate, int position, int total_pos_coverage, std::vector<VariantPtr>& var,
                         std::vector<std::string>& debug_lines, int hicov);
    double collect_vars_at_position(AlignedVars& aligned, int position, const std::vector<VariantPtr>& var);
    void collect_reference_variants(int position, int total_pos_coverage, Vars& vars,
                                    std::vector<std::string>& debug_lines);
    void update_ref_variant(int position, int total_pos_coverage, Variant& vref, std::vector<std::string>& debug_lines,
                            int rfc, int rrc);
    MsiResult proceed_vref_is_deletion(int position, int dellen);
    MsiResult proceed_vref_is_insertion(int position, const std::string& vn);
    void construct_debug_lines(const std::vector<std::string>& debug_lines, Variant& vref) const;
    std::string debug_content(const std::string& n, const Variant& v) const;

    VariationData& data_;
    const RegionScope& scope_;
    const Config& conf_;
    const RefMap& ref_;
};

}  // namespace

MsiResult find_msi(const std::string& tseq1, const std::string& tseq2, const std::string* left) {
    int nmsi = 1;
    int shift3 = 0;
    std::string maxmsi;
    double msicnt = 0;
    while (nmsi <= static_cast<int>(tseq1.size()) && nmsi <= 6) {
        const std::string msint = util::substr(tseq1, -nmsi, nmsi);
        std::string msimatch = repeat_suffix(tseq1, msint);
        if (left != nullptr && !left->empty()) {
            // a non-empty suffix of tseq1 is always found again in left + tseq1
            std::string m = repeat_suffix(*left + tseq1, msint);
            if (!m.empty()) {
                msimatch = std::move(m);
            }
        }
        double curmsi = static_cast<double>(msimatch.size()) / nmsi;
        const std::size_t pre = repeat_prefix_length(tseq2, msint);
        if (pre > 0) {
            curmsi += static_cast<double>(pre) / nmsi;
        }
        if (curmsi > msicnt) {
            maxmsi = msint;
            msicnt = curmsi;
        }
        nmsi++;
    }
    const std::string tseq = tseq1 + tseq2;
    while (shift3 < static_cast<int>(tseq2.size()) &&
           tseq[static_cast<std::size_t>(shift3)] == tseq2[static_cast<std::size_t>(shift3)]) {
        shift3++;
    }
    return MsiResult{msicnt, shift3, maxmsi};
}

AlignedVars build_vars(VariationData& data, const RegionScope& scope) { return Builder(data, scope).run(); }

namespace {

AlignedVars Builder::run() {
    if (conf_.y) {
        std::println(stderr, "Current segment: {}:{}-{} ", scope_.region.chr, scope_.region.start, scope_.region.end);
    }
    AlignedVars aligned;
    int last_position = 0;
    data_.non_insertion_variants.for_each([&](const int& position, VariationMap& vars_at_cur_position) {
        try {
            last_position = position;
            if (vars_at_cur_position.empty() && !data_.insertion_variants.contains(position)) {
                return;
            }
            if (!vars_at_cur_position.sv || conf_.delete_duplicate_variants) {
                if (position < scope_.region.start || position > scope_.region.end) {
                    return;
                }
            }
            if (!vars_at_cur_position.sv && !data_.ref_coverage.contains(position)) {
                return;
            }
            if (is_the_same_variation_on_ref(position, vars_at_cur_position)) {
                return;
            }
            const int* cov = data_.ref_coverage.get(position);
            if (cov == nullptr || *cov == 0) {
                if (!vars_at_cur_position.sv) {
                    util::throw_npe("varsAtCurPosition.sv.type");
                }
                std::println(stderr, "Error tcov: {} {} {} {} {}", scope_.region.chr, position, scope_.region.start,
                             scope_.region.end, vars_at_cur_position.sv->type);
                return;
            }
            int total_pos_coverage = *cov;
            const int hicov = calc_hicov(vars_at_cur_position);

            std::vector<VariantPtr> var;
            std::vector<std::string> keys = vars_at_cur_position.keys();
            std::sort(keys.begin(), keys.end());
            std::vector<std::string> debug_lines;

            create_variant(data_.duprate, aligned, position, vars_at_cur_position, total_pos_coverage, var, debug_lines,
                           keys, hicov);
            total_pos_coverage = create_insertion(data_.duprate, position, total_pos_coverage, var, debug_lines, hicov);
            // sortVariants: by meanQuality * positionCoverage descending, then description
            std::stable_sort(var.begin(), var.end(), [](const VariantPtr& a, const VariantPtr& b) {
                const int c = util::double_compare(b->mean_quality * b->position_coverage,
                                                   a->mean_quality * a->position_coverage);
                if (c != 0) {
                    return c < 0;
                }
                return a->description_string < b->description_string;
            });

            const double maxfreq = collect_vars_at_position(aligned, position, var);
            if (maxfreq <= conf_.freq && !scope_.run.amplicon_based_calling) {
                if (!conf_.bam.has_bam2()) {
                    aligned.erase(position);
                    return;
                }
            }
            Vars& variations_at_pos = get_or_put_vars(aligned, position);
            collect_reference_variants(position, total_pos_coverage, variations_at_pos, debug_lines);
        } catch (const util::RecoverableError& e) {
            print_exception_and_continue(e, "position", std::to_string(last_position), &scope_.region,
                                         scope_.exceptions);
        }
    });
    return aligned;
}

bool Builder::is_the_same_variation_on_ref(int position, const VariationMap& vars) {
    std::set<std::string> vk;
    vars.for_each([&](const std::string& k, const Variation&) { vk.insert(k); });
    if (data_.insertion_variants.contains(position)) {
        vk.insert("I");
    }
    if (vk.size() == 1 && ref_.contains(position) && vk.contains(std::string(1, ref_.at(position)))) {
        if (!conf_.bam.has_bam2() && !scope_.run.amplicon_based_calling) {
            return true;
        }
    }
    return false;
}

int Builder::calc_hicov(const VariationMap& non_insertion) {
    int hicov = 0;
    non_insertion.for_each([&](const std::string& key, const Variation& v) {
        if (key == "SV" || key.starts_with("+")) {
            return;
        }
        hicov += v.high_quality_reads_count;
    });
    return hicov;
}

void Builder::create_variant(double duprate, AlignedVars& aligned, int position, const VariationMap& vars,
                             int total_pos_coverage, std::vector<VariantPtr>& var,
                             std::vector<std::string>& debug_lines, const std::vector<std::string>& keys, int hicov) {
    for (const std::string& description_string : keys) {
        if (description_string == "SV") {
            const SvInfo& sv = vars.sv.value();
            get_or_put_vars(aligned, position).sv =
                std::to_string(sv.splits) + "-" + std::to_string(sv.pairs) + "-" + std::to_string(sv.clusters);
            continue;
        }
        const Variation& cnt = *vars.get(description_string);
        if (cnt.vars_count == 0) {
            continue;
        }
        const int fwd = cnt.get_dir(false);
        const int rev = cnt.get_dir(true);
        const int bias = strand_bias(fwd, rev, conf_);
        const double base_quality = round_half_even(cnt.mean_quality / cnt.vars_count, 1);
        const double mapping_quality = round_half_even(cnt.mean_mapping_quality / cnt.vars_count, 1);
        const int hicnt = cnt.high_quality_reads_count;
        const int locnt = cnt.low_quality_reads_count;
        int ttcov = total_pos_coverage;
        if (cnt.vars_count > total_pos_coverage && cnt.extracnt > 0 &&
            cnt.vars_count - total_pos_coverage < cnt.extracnt) {
            ttcov = cnt.vars_count;
        }
        auto tvref = std::make_shared<Variant>();
        tvref->description_string = description_string;
        tvref->position_coverage = cnt.vars_count;
        tvref->vars_count_on_forward = fwd;
        tvref->vars_count_on_reverse = rev;
        tvref->strand_bias_flag = std::to_string(bias);
        tvref->frequency = round_half_even(cnt.vars_count / static_cast<double>(ttcov), 4);
        tvref->mean_position = round_half_even(cnt.mean_position / cnt.vars_count, 1);
        tvref->is_at_least_at_2_positions = cnt.pstd;
        tvref->mean_quality = base_quality;
        tvref->has_at_least_2_diff_qualities = cnt.qstd;
        tvref->mean_mapping_quality = mapping_quality;
        tvref->high_quality_to_low_quality_ratio = hicnt / (locnt != 0 ? static_cast<double>(locnt) : 0.5);
        tvref->high_quality_reads_frequency = hicov > 0 ? hicnt / static_cast<double>(hicov) : 0;
        tvref->extra_frequency = cnt.extracnt != 0 ? cnt.extracnt / static_cast<double>(ttcov) : 0;
        tvref->shift3 = 0;
        tvref->msi = 0;
        tvref->number_of_mismatches = round_half_even(cnt.number_of_mismatches / cnt.vars_count, 1);
        tvref->hicnt = hicnt;
        tvref->hicov = hicov;
        tvref->duprate = duprate;
        var.push_back(tvref);
        if (conf_.debug) {
            debug_lines.push_back(debug_content(description_string, *tvref));
        }
    }
}

int Builder::create_insertion(double duprate, int position, int total_pos_coverage, std::vector<VariantPtr>& var,
                              std::vector<std::string>& debug_lines, int hicov) {
    VariationMap* insertion_variations = data_.insertion_variants.get(position);
    if (insertion_variations == nullptr) {
        return total_pos_coverage;
    }
    std::vector<std::string> descs = insertion_variations->keys();
    std::sort(descs.begin(), descs.end());
    for (const std::string& description_string : descs) {
        if (description_string.contains('&') && data_.ref_coverage.contains(position + 1)) {
            total_pos_coverage = *data_.ref_coverage.get(position + 1);
        }
        const Variation& cnt = *insertion_variations->get(description_string);
        const int fwd = cnt.get_dir(false);
        const int rev = cnt.get_dir(true);
        const int bias = strand_bias(fwd, rev, conf_);
        const double vqual = round_half_even(cnt.mean_quality / cnt.vars_count, 1);
        const double mq = round_half_even(cnt.mean_mapping_quality / static_cast<double>(cnt.vars_count), 1);
        const int hicnt = cnt.high_quality_reads_count;
        const int locnt = cnt.low_quality_reads_count;

        int ttcov = total_pos_coverage;
        if (cnt.vars_count > total_pos_coverage && cnt.extracnt != 0 &&
            cnt.vars_count - total_pos_coverage < cnt.extracnt) {
            ttcov = cnt.vars_count;
        }
        if (ttcov < cnt.vars_count) {
            ttcov = cnt.vars_count;
            const int* next_cov = data_.ref_coverage.get(position + 1);
            if (next_cov != nullptr && ttcov < *next_cov - cnt.vars_count) {
                ttcov = *next_cov;
                if (Variation* v =
                        get_variation_maybe(data_.non_insertion_variants, position + 1, ref_get(ref_, position + 1))) {
                    v->vars_count_on_forward -= fwd;
                    v->vars_count_on_reverse -= rev;
                }
            }
            total_pos_coverage = ttcov;
        }

        auto tvref = std::make_shared<Variant>();
        if (hicov < hicnt) {
            hicov = hicnt;
        }
        tvref->description_string = description_string;
        tvref->position_coverage = cnt.vars_count;
        tvref->vars_count_on_forward = fwd;
        tvref->vars_count_on_reverse = rev;
        tvref->strand_bias_flag = std::to_string(bias);
        tvref->frequency = round_half_even(cnt.vars_count / static_cast<double>(ttcov), 4);
        tvref->mean_position = round_half_even(cnt.mean_position / static_cast<double>(cnt.vars_count), 1);
        tvref->is_at_least_at_2_positions = cnt.pstd;
        tvref->mean_quality = vqual;
        tvref->has_at_least_2_diff_qualities = cnt.qstd;
        tvref->mean_mapping_quality = mq;
        tvref->high_quality_to_low_quality_ratio = hicnt / (locnt != 0 ? static_cast<double>(locnt) : 0.5);
        tvref->high_quality_reads_frequency = hicov > 0 ? hicnt / static_cast<double>(hicov) : 0;
        tvref->extra_frequency = cnt.extracnt != 0 ? cnt.extracnt / static_cast<double>(ttcov) : 0;
        tvref->shift3 = 0;
        tvref->msi = 0;
        tvref->number_of_mismatches =
            round_half_even(cnt.number_of_mismatches / static_cast<double>(cnt.vars_count), 1);
        tvref->hicnt = hicnt;
        tvref->hicov = hicov;
        tvref->duprate = duprate;
        var.push_back(tvref);
        if (conf_.debug) {
            debug_lines.push_back("I" + debug_content(description_string, *tvref));
        }
    }
    return total_pos_coverage;
}

std::string Builder::debug_content(const std::string& n, const Variant& v) const {
    return n + ":" + std::to_string(v.vars_count_on_forward + v.vars_count_on_reverse) + ":F-" +
           std::to_string(v.vars_count_on_forward) + ":R-" + std::to_string(v.vars_count_on_reverse) + ":" +
           util::decimal_format(v.frequency, 4) + ":" + v.strand_bias_flag + ":" +
           util::decimal_format(v.mean_position, 1) + ":" + (v.is_at_least_at_2_positions ? "1" : "0") + ":" +
           util::decimal_format(v.mean_quality, 1) + ":" + (v.has_at_least_2_diff_qualities ? "1" : "0") + ":" +
           util::decimal_format(v.high_quality_reads_frequency, 4) + ":" +
           util::double_to_string(v.mean_mapping_quality) + ":" +
           util::decimal_format(v.high_quality_to_low_quality_ratio, 3);
}

double Builder::collect_vars_at_position(AlignedVars& aligned, int position, const std::vector<VariantPtr>& var) {
    double maxfreq = 0;
    const std::string ref_base = ref_or_null(ref_, position);  // String.valueOf(ref.get(position))
    for (const VariantPtr& tvar : var) {
        if (tvar->description_string == ref_base) {
            get_or_put_vars(aligned, position).reference_variant = tvar;
        } else {
            Vars& v = get_or_put_vars(aligned, position);
            v.variants.push_back(tvar);
            v.var_description_string_to_variants.put(tvar->description_string, tvar);
            if (tvar->frequency > maxfreq) {
                maxfreq = tvar->frequency;
            }
        }
    }
    return maxfreq;
}

MsiResult Builder::proceed_vref_is_deletion(int position, int dellen) {
    const std::string leftseq = join_ref(ref_, std::max(position - kRef70Bases, 1), position - 1);
    const int chr0 = chr_length();
    const std::string tseq = join_ref(ref_, position, std::min(position + dellen + kRef70Bases, chr0));
    MsiResult r = find_msi(util::substr(tseq, 0, dellen), util::substr(tseq, dellen), &leftseq);
    MsiResult without_left = find_msi(leftseq, util::substr(tseq, dellen), nullptr);
    if (r.msi < without_left.msi) {
        r.msi = without_left.msi;
        r.msint = without_left.msint;  // don't change shift3
    }
    if (r.msi <= r.shift3 / static_cast<double>(dellen)) {
        r.msi = r.shift3 / static_cast<double>(dellen);
    }
    return r;
}

MsiResult Builder::proceed_vref_is_insertion(int position, const std::string& vn) {
    const std::string tseq1 = util::substring(vn, 1);
    const std::string leftseq = join_ref(ref_, std::max(position - kRef50Bases, 1), position);
    const int x = chr_length();
    const std::string tseq2 = join_ref(ref_, position + 1, std::min(position + kRef70Bases, x));
    MsiResult r = find_msi(tseq1, tseq2, &leftseq);
    MsiResult without_left = find_msi(leftseq, tseq2, nullptr);
    if (r.msi < without_left.msi) {
        r.msi = without_left.msi;
        r.msint = without_left.msint;
    }
    if (r.msi <= r.shift3 / static_cast<double>(tseq1.size())) {
        r.msi = r.shift3 / static_cast<double>(tseq1.size());
    }
    return r;
}

void Builder::collect_reference_variants(int position, int total_pos_coverage, Vars& vars,
                                         std::vector<std::string>& debug_lines) {
    int reference_forward_coverage = 0;
    int reference_reverse_coverage = 0;
    std::string genotype1;
    if (vars.reference_variant && vars.reference_variant->frequency >= conf_.freq) {
        genotype1 = vars.reference_variant->description_string;
    } else if (!vars.variants.empty()) {
        genotype1 = vars.variants[0]->description_string;
    } else {
        if (!vars.reference_variant) {
            util::throw_npe("variationsAtPos.referenceVariant");
        }
        genotype1 = vars.reference_variant->description_string;
    }
    if (vars.reference_variant) {
        reference_forward_coverage = vars.reference_variant->vars_count_on_forward;
        reference_reverse_coverage = vars.reference_variant->vars_count_on_reverse;
    }
    if (genotype1.starts_with("+")) {
        if (auto d = patterns::dup_num(genotype1)) {
            genotype1 = "+" + std::to_string(Config::kSvFlank + to_int(*d));
        } else {
            genotype1 = "+" + std::to_string(genotype1.size() - 1);
        }
    }
    std::string genotype2;

    const int ref_cov_here = *data_.ref_coverage.get(position);
    if (total_pos_coverage > ref_cov_here && data_.non_insertion_variants.contains(position + 1) &&
        ref_.contains(position + 1) &&
        data_.non_insertion_variants.get(position + 1)->contains(std::string(1, ref_.at(position + 1)))) {
        Variation* tpref = get_variation_maybe(data_.non_insertion_variants, position + 1, ref_get(ref_, position + 1));
        reference_forward_coverage = tpref->vars_count_on_forward;
        reference_reverse_coverage = tpref->vars_count_on_reverse;
    }

    if (!vars.variants.empty()) {
        for (const VariantPtr& vref_ptr : vars.variants) {
            Variant& vref = *vref_ptr;
            std::string genotype1current = genotype1;
            genotype2 = vref.description_string;
            if (genotype2.starts_with("+")) {
                genotype2 = "+" + std::to_string(genotype2.size() - 1);
            }
            const std::string description_string = vref.description_string;
            int deletion_length = 0;
            if (auto d = patterns::begin_minus_number(description_string)) {
                deletion_length = to_int(*d);
            }
            int end_position = position;
            if (description_string.starts_with("-")) {
                end_position = position + deletion_length - 1;
            }
            std::string refallele;
            std::string varallele;
            int shift3 = 0;
            double msi = 0;
            std::string msint;
            int start_position = position;

            if (description_string.starts_with("+")) {
                if (!description_string.contains('&') && !description_string.contains('#') &&
                    !description_string.contains("<dup")) {
                    MsiResult r = proceed_vref_is_insertion(position, description_string);
                    msi = r.msi;
                    shift3 = r.shift3;
                    msint = r.msint;
                }
                if (conf_.move_indels_to_3) {
                    start_position += shift3;
                    end_position += shift3;
                }
                refallele = ref_or_empty(ref_, position);
                varallele = refallele + util::substring(description_string, 1);
                if (static_cast<int>(varallele.size()) > conf_.sv_min_len) {
                    end_position += static_cast<int>(varallele.size());
                    varallele = "<DUP>";
                }
                if (auto d = patterns::dup_num(varallele)) {
                    const int dup_count = to_int(*d);
                    end_position = start_position + (2 * Config::kSvFlank + dup_count) - 1;
                    genotype2 = "+" + std::to_string(2 * Config::kSvFlank + dup_count);
                    varallele = "<DUP>";
                }
            } else if (description_string.starts_with("-")) {
                const bool inv = patterns::inv_num(description_string).has_value();
                const bool start_minus_num_caret = patterns::begin_minus_number_caret(description_string);
                if (deletion_length < conf_.sv_min_len) {
                    varallele = patterns::remove_begin_minus_number(description_string);
                    MsiResult r = proceed_vref_is_deletion(position, deletion_length);
                    msi = r.msi;
                    shift3 = r.shift3;
                    msint = r.msint;
                    if (inv) {
                        varallele = "<INV>";
                        genotype2 = "<INV" + std::to_string(deletion_length) + ">";
                    }
                } else if (start_minus_num_caret) {
                    varallele = "<INV>";
                    genotype2 = "<INV" + std::to_string(deletion_length) + ">";
                } else {
                    varallele = "<DEL>";
                }
                if (!description_string.contains('&') && !description_string.contains('#') &&
                    !description_string.contains('^')) {
                    if (conf_.move_indels_to_3) {
                        start_position += shift3;
                    }
                    if (varallele != "<DEL>") {
                        varallele = ref_or_empty(ref_, position - 1);
                    }
                    refallele = ref_or_empty(ref_, position - 1);
                    start_position--;
                }
                if (patterns::some_sv_numbers(description_string)) {
                    refallele = ref_or_empty(ref_, position);
                } else if (deletion_length < conf_.sv_min_len) {
                    refallele += join_ref(ref_, position, position + deletion_length - 1);
                }
            } else {
                // SNP or MNP
                const std::string tseq1 = join_ref(ref_, std::max(position - kRef30Bases, 1), position + 1);
                const int chr0 = chr_length();
                const std::string tseq2 = join_ref(ref_, position + 2, std::min(position + kRef70Bases, chr0));
                MsiResult r = find_msi(tseq1, tseq2, nullptr);
                msi = r.msi;
                shift3 = r.shift3;
                msint = r.msint;
                refallele = ref_or_empty(ref_, position);
                varallele = description_string;
            }

            if (auto extra = patterns::amp_atgc(description_string)) {
                varallele = util::replace_first(varallele, "&", "");
                std::string tch = join_ref(ref_, end_position + 1, end_position + static_cast<int>(extra->size()));
                refallele += tch;
                genotype1current += tch;
                end_position += static_cast<int>(extra->size());
                if (auto vextra = patterns::amp_atgc(varallele)) {
                    varallele = util::replace_first(varallele, "&", "");
                    tch = join_ref(ref_, end_position + 1, end_position + static_cast<int>(vextra->size()));
                    refallele += tch;
                    genotype1current += tch;
                    end_position += static_cast<int>(vextra->size());
                }
                if (description_string.starts_with("+")) {
                    refallele = util::substring(refallele, 1);
                    varallele = util::substring(varallele, 1);
                    start_position++;
                }
                if (varallele == "<DEL>" && !refallele.empty()) {
                    refallele = ref_or_empty(ref_, start_position);
                    if (const int* c = data_.ref_coverage.get(start_position - 1)) {
                        total_pos_coverage = *c;
                    }
                    if (vref.position_coverage > total_pos_coverage) {
                        total_pos_coverage = vref.position_coverage;
                    }
                    vref.frequency = vref.position_coverage / static_cast<double>(total_pos_coverage);
                }
            }

            if (auto groups = patterns::hash_group_caret_group(description_string)) {
                const std::string& matched_sequence = groups->first;
                const std::string& tail = groups->second;
                end_position += static_cast<int>(matched_sequence.size());
                refallele += join_ref(ref_, end_position - static_cast<int>(matched_sequence.size()) + 1, end_position);
                if (auto d = patterns::begin_digits(tail)) {
                    const int deletion = to_int(*d);
                    refallele += join_ref(ref_, end_position + 1, end_position + deletion);
                    end_position += deletion;
                }
                varallele = patterns::remove_first_caret_digits(util::replace_first(varallele, "#", ""));
                genotype1current = util::replace_first(util::replace_first(genotype1current, "#", "m"), "^", "i");
                genotype2 = util::replace_first(util::replace_first(genotype2, "#", "m"), "^", "i");
            }
            if (patterns::caret_atgnc(description_string)) {
                varallele = util::replace_first(varallele, "^", "");
                genotype1current = util::replace_first(genotype1current, "^", "i");
                genotype2 = util::replace_first(genotype2, "^", "i");
            }

            // Perform adjustment to get as close to CRISPR site as possible
            const int cut_site = conf_.crispr_cutting_site;
            if (cut_site != 0 && refallele.size() > 1 && varallele.size() > 1) {
                long n = 0;
                while (static_cast<long>(refallele.size()) > n + 1 && static_cast<long>(varallele.size()) > n + 1 &&
                       util::substr(refallele, n, 1) == util::substr(varallele, n, 1)) {
                    n++;
                }
                if (n != 0) {
                    start_position += static_cast<int>(n);
                    refallele = util::substr(refallele, n);
                    varallele = util::substr(varallele, n);
                }
            }
            if (cut_site != 0 && (refallele.size() != varallele.size() &&
                                  util::substr(refallele, 0, 1) == util::substr(varallele, 0, 1))) {
                if (!(start_position == cut_site || end_position == cut_site)) {
                    int n = 0;
                    const int dis = std::min(std::abs(cut_site - start_position), std::abs(cut_site - end_position));
                    if (start_position < cut_site) {
                        while (start_position + n < cut_site && n < shift3 && end_position + n != cut_site) {
                            n++;
                        }
                        if (std::abs(start_position + n - cut_site) > dis &&
                            std::abs(end_position + n - cut_site) > dis) {
                            n = 0;
                        }
                    }
                    if (end_position < cut_site && n == 0) {
                        if (std::abs(end_position - cut_site) <= std::abs(start_position - cut_site)) {
                            while (end_position + n < cut_site && n < shift3) {
                                n++;
                            }
                        }
                    }
                    if (n > 0) {
                        start_position += n;
                        end_position += n;
                        refallele.clear();
                        for (int i = start_position; i <= end_position; i++) {
                            refallele += ref_or_null(ref_, i);
                        }
                        std::string tva;
                        if (refallele.size() < varallele.size()) {  // Insertion
                            tva = util::substr(varallele, 1);
                            if (tva.size() > 1) {
                                const int ttn = n % static_cast<int>(tva.size());
                                if (ttn != 0) {
                                    tva = util::substr(tva, ttn) + util::substr(tva, 0, ttn);
                                }
                            }
                        }
                        varallele = ref_or_null(ref_, start_position) + tva;
                    }
                    vref.crispr = n;
                }
            }

            vref.leftseq = join_ref(ref_, std::max(start_position - kRef20Bases, 1), start_position - 1);
            const int chr0 = chr_length();
            vref.rightseq = join_ref(ref_, end_position + 1, std::min(end_position + kRef20Bases, chr0));
            std::string genotype = genotype1current + "/" + genotype2;
            genotype = util::replace_all(util::replace_all(util::replace_all(genotype, "&", ""), "#", ""), "^", "i");
            vref.extra_frequency = round_half_even(vref.extra_frequency, 4);
            vref.frequency = round_half_even(vref.frequency, 4);
            vref.high_quality_reads_frequency = round_half_even(vref.high_quality_reads_frequency, 4);
            vref.msi = round_half_even(msi, 3);
            vref.msint = static_cast<int>(msint.size());
            vref.shift3 = shift3;
            vref.start_position = start_position;
            vref.end_position = end_position;
            vref.refallele = validate_refallele(refallele);
            vref.varallele = varallele;
            vref.genotype = genotype;
            vref.total_pos_coverage = total_pos_coverage;
            vref.ref_forward_coverage = reference_forward_coverage;
            vref.ref_reverse_coverage = reference_reverse_coverage;
            if (vars.reference_variant) {
                vref.strand_bias_flag = vars.reference_variant->strand_bias_flag + ";" + vref.strand_bias_flag;
            } else {
                vref.strand_bias_flag = "0;" + vref.strand_bias_flag;
            }
            // adjustVariantCounts
            const std::string message = "column in variant on position: " + std::to_string(position) + " " +
                                        vref.refallele + "->" + vref.varallele + " was negative, adjusted to zero.";
            if (vref.ref_forward_coverage < 0) {
                vref.ref_forward_coverage = 0;
                std::println(stderr, "Reference forward count {}", message);
            }
            if (vref.ref_reverse_coverage < 0) {
                vref.ref_reverse_coverage = 0;
                std::println(stderr, "Reference reverse count {}", message);
            }
            if (vref.vars_count_on_forward < 0) {
                vref.vars_count_on_forward = 0;
                std::println(stderr, "Variant forward count {}", message);
            }
            if (vref.vars_count_on_reverse < 0) {
                vref.vars_count_on_reverse = 0;
                std::println(stderr, "Variant reverse count {}", message);
            }
            construct_debug_lines(debug_lines, vref);
        }
        if (conf_.disable_sv) {
            std::erase_if(vars.variants,
                          [](const VariantPtr& v) { return patterns::any_sv(v->varallele).has_value(); });
        }
    } else if (vars.reference_variant) {
        update_ref_variant(position, total_pos_coverage, *vars.reference_variant, debug_lines,
                           reference_forward_coverage, reference_reverse_coverage);
    } else {
        vars.reference_variant = std::make_shared<Variant>();
    }
}

void Builder::update_ref_variant(int position, int total_pos_coverage, Variant& vref,
                                 std::vector<std::string>& debug_lines, int rfc, int rrc) {
    vref.total_pos_coverage = total_pos_coverage;
    vref.position_coverage = 0;
    vref.frequency = 0;
    vref.ref_forward_coverage = rfc;
    vref.ref_reverse_coverage = rrc;
    vref.vars_count_on_forward = 0;
    vref.vars_count_on_reverse = 0;
    vref.msi = 0;
    vref.msint = 0;
    if (!vref.strand_bias_flag.contains(';')) {
        vref.strand_bias_flag += ";0";
    }
    vref.shift3 = 0;
    vref.start_position = position;
    vref.end_position = position;
    vref.high_quality_reads_frequency = round_half_even(vref.high_quality_reads_frequency, 4);
    const std::string reference_base = ref_or_empty(ref_, position);
    vref.refallele = validate_refallele(reference_base);
    vref.varallele = validate_refallele(reference_base);
    vref.genotype = reference_base + "/" + reference_base;
    vref.leftseq.clear();
    vref.rightseq.clear();
    vref.duprate = data_.duprate;
    construct_debug_lines(debug_lines, vref);
}

void Builder::construct_debug_lines(const std::vector<std::string>& debug_lines, Variant& vref) const {
    if (conf_.debug) {
        std::string sb;
        for (const auto& s : debug_lines) {
            if (!sb.empty()) {
                sb += " & ";
            }
            sb += s;
        }
        vref.debug = sb;
    }
}

}  // namespace

}  // namespace senri::core
