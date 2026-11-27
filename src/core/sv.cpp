
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <print>

#include "core/realigner.hpp"
#include "io/reference.hpp"
#include "util/checked_string.hpp"
#include "util/sequence.hpp"

namespace senri::core {

namespace {

// SequenceUtil.reverseComplement
std::string reverse_complement(const std::string& s) { return util::complement(util::reverse(s)); }

std::string i2s(int v) { return std::to_string(v); }

// VariationUtils.isHasAndEquals(index, ref, index2): both bases loaded and equal.
bool ref_equal(const RefMap& ref, int i, int j) {
    auto a = ref_get(ref, i);
    auto b = ref_get(ref, j);
    return a && b && *a == *b;
}

// Map.get(k) unboxed: NPE when absent.
int cov_at(util::HashOrderMap<int, int>& cov, int k) {
    const int* v = cov.get(k);
    if (v == nullptr) {
        util::throw_npe("refCoverage.get(" + std::to_string(k) + ")");
    }
    return *v;
}

// LinkedHashMap<Integer, Integer> entries sorted by value, highest first (stable); first key or 0.
int top_soft(const util::LinkedMap<int, int>& soft) {
    std::vector<std::pair<int, int>> v;
    soft.for_each([&](int k, int c) { v.emplace_back(k, c); });
    std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    return v.empty() ? 0 : v[0].first;
}

// StructuralVariantsProcessor.PairsData
struct PairsData {
    int pairs = 0;
    double pmean = 0;
    double qmean = 0;
    double big_qmean = 0;
    double nm = 0;
};

enum class Side { five, three };

class SvProcessor {
public:
    SvProcessor(VariationData& data, RegionScope& scope)
        : d_(data),
          scope_(scope),
          conf_(scope.conf),
          ref_(scope.reference.reference_sequences),
          region_(scope.region),
          max_read_length_(data.max_read_length) {}

    void find_all_svs() {
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants: DEL\n");
        }
        find_del();
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants: INV\n");
        }
        find_inv();
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants\n");
        }
        findsv();
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants: DEL discordant pairs only\n");
        }
        find_del_disc();
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants: INV discordant pairs only\n");
        }
        find_inv_disc();
        if (conf_.y) {
            std::println(stderr, "Start Structural Variants: DUP discordant pairs only\n");
        }
        find_dup_disc();
    }

    void output_clipping();

private:
    void find_del();
    void find_inv();
    void find_inv_sub(std::vector<Sclip>& svref, int dir, Side side);
    void findsv();
    void find_del_disc();
    void find_inv_disc();
    void find_dup_disc();
    PairsData check_pairs(int start, int end, const std::vector<std::vector<Sclip>*>& sv);
    Match find_match_rev(std::string seq, int position, int dir, int seed = Config::kSeed1, int mm = 3) const;
    Match match(const std::string& seq, int position, int dir, int seed = Config::kSeed1, int mm = 3) const {
        return find_match(seq, scope_.reference, position, dir, seed, mm, conf_);
    }
    std::string conseq(Sclip& sc) const { return findconseq(sc, 0, scope_.run, conf_); }
    // isLoaded(region.chr, s, e, reference)
    [[nodiscard]] bool loaded(int s, int e) const { return io::is_loaded(region_.chr, s, e, scope_.reference); }
    void reference(int s, int e, int extension) { extend_reference(scope_, region_.with_bounds(s, e), extension); }
    void partial(int s, int e) { partial_pipeline(d_, scope_, region_.with_bounds(s, e), max_read_length_); }
    void exception(const std::exception& e, int last) const {
        print_exception_and_continue(e, "variant", std::to_string(last), &region_, scope_.exceptions);
    }
    // The pair-support Variation of findDEL / findsv: mcnt pairs, the cluster's means scaled.
    static Variation pair_support(int mcnt, const Sclip& del);

    VariationData& d_;
    RegionScope& scope_;
    const Config& conf_;
    const RefMap& ref_;
    const Region& region_;
    int max_read_length_;
};

Variation SvProcessor::pair_support(int mcnt, const Sclip& del) {
    Variation tv;
    tv.vars_count = mcnt;
    tv.high_quality_reads_count = mcnt;
    tv.vars_count_on_forward = mcnt / 2;
    tv.vars_count_on_reverse = mcnt - mcnt / 2;
    tv.mean_quality = del.mean_quality * mcnt / del.vars_count;
    tv.mean_position = del.mean_position * mcnt / del.vars_count;
    tv.mean_mapping_quality = del.mean_mapping_quality * mcnt / del.vars_count;
    tv.number_of_mismatches = del.number_of_mismatches * mcnt / del.vars_count;
    return tv;
}

void SvProcessor::find_del() {
    PositionMap& niv = d_.non_insertion_variants;
    const RefMap& ref = ref_;
    auto& cov = d_.ref_coverage;
    int last_start = 0;
    for (Sclip& del : d_.sv_structures.svfdel) {
        try {
            last_start = del.start;
            if (del.used) {
                continue;
            }
            if (del.vars_count < conf_.minr) {
                continue;
            }
            int softp = top_soft(del.soft);
            if (conf_.y) {
                std::println(stderr, "\n\nWorking DEL 5' {} mate cluster cnt: {}", softp, del.vars_count);
            }
            if (softp != 0) {
                Sclip* scv = d_.soft_clips_3end.get(softp);
                if (scv == nullptr) {
                    continue;
                }
                if (scv->used) {
                    continue;
                }
                const std::string seq = conseq(*scv);
                if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                    continue;
                }
                if (!loaded(del.mstart, del.mend)) {
                    reference(del.mstart, del.mend, 300);
                    partial(del.mstart - 200, del.mend + 200);
                }
                Match m = match(seq, softp, 1);
                int bp = m.base_position;
                if (bp == 0) {
                    continue;
                }
                if (!(bp - softp > 30 && is_overlap(softp, bp, del.end, del.mstart, max_read_length_))) {
                    continue;
                }
                bp--;
                const int dellen = bp - softp + 1;
                while (ref.contains(bp) && ref.contains(softp - 1) && ref.at(bp) == ref.at(softp - 1)) {
                    bp--;
                    if (bp != 0) {
                        softp--;
                    }
                }
                Variation& variation = get_variation(niv, softp, "-" + i2s(dellen));
                variation.vars_count = 0;

                SvInfo& sv = get_sv(niv, softp);
                sv.type = "DEL";
                sv.pairs += del.vars_count;
                sv.splits += scv->vars_count;
                sv.clusters++;

                if (!(cov.contains(softp) && cov_at(cov, softp) > del.vars_count)) {
                    cov.put(softp, del.vars_count);
                }
                if (cov.contains(bp) && cov_at(cov, softp) < cov_at(cov, bp)) {
                    cov.put(softp, cov_at(cov, bp));
                }
                Variation* lref = nullptr;
                if (VariationMap* at = niv.get(softp); at != nullptr && ref.contains(softp)) {
                    lref = at->get(std::string(1, ref.at(softp)));
                }
                adj_cnt(variation, *scv, lref);

                const int mcnt = del.vars_count;
                adj_cnt(variation, pair_support(mcnt, del));
                del.used = true;
                (void)mark_sv(softp, bp, {&d_.sv_structures.svrdel}, max_read_length_, conf_);
                if (conf_.y) {
                    std::println(stderr, "    Found DEL SV from 5' softclip unhappy reads: {} -{} Cnt: {} AdjCnt: {}",
                                 bp, dellen, del.vars_count, mcnt);
                }
            } else {  // Look within a read length
                if (conf_.y) {
                    std::println(stderr, "\n\nWorking DEL 5' no softp mate cluster cnt: {}", del.vars_count);
                }
                if (!loaded(del.mstart, del.mend)) {
                    reference(del.mstart, del.mend, 300);
                }
                bool done = false;
                std::vector<int> keys = d_.soft_clips_3end.keys();
                for (int i : keys) {
                    if (done) {
                        break;
                    }
                    Sclip& scv = *d_.soft_clips_3end.get(i);
                    if (scv.used) {
                        continue;
                    }
                    if (!(i >= del.end - 3 && i - del.end < 3 * max_read_length_)) {
                        continue;
                    }
                    const std::string seq = conseq(scv);
                    if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                        continue;
                    }
                    softp = i;
                    Match m = match(seq, softp, 1);
                    int bp = m.base_position;
                    if (bp == 0) {
                        m = match(seq, softp, 1, Config::kSeed2, 0);
                        bp = m.base_position;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    if (!(bp - softp > 30 && is_overlap(softp, bp, del.end, del.mstart, max_read_length_))) {
                        continue;
                    }
                    bp--;
                    const int dellen = bp - softp + 1;
                    Variation& variation = get_variation(niv, softp, "-" + i2s(dellen));
                    variation.vars_count = 0;

                    SvInfo& sv = get_sv(niv, softp);
                    sv.type = "DEL";
                    sv.pairs += del.vars_count;
                    sv.splits += scv.vars_count;
                    sv.clusters++;

                    if (!(cov.contains(softp) && cov_at(cov, softp) > del.vars_count)) {
                        cov.put(softp, del.vars_count);
                    }
                    if (cov.contains(bp) && cov_at(cov, softp) < cov_at(cov, bp)) {
                        cov.put(softp, cov_at(cov, bp));
                    }
                    adj_cnt(variation, scv);
                    const int mcnt = del.vars_count;
                    adj_cnt(variation, pair_support(mcnt, del));
                    del.used = true;
                    (void)mark_sv(softp, bp, {&d_.sv_structures.svrdel}, max_read_length_, conf_);
                    if (conf_.y) {
                        std::println(stderr, "    Found DEL SV from 5' softclip happy reads: {} -{} Cnt: {} AdjCnt: {}",
                                     bp, dellen, del.vars_count, mcnt);
                    }
                    done = true;
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
    for (Sclip& del : d_.sv_structures.svrdel) {
        try {
            last_start = del.start;
            if (del.used) {
                continue;
            }
            if (del.vars_count < conf_.minr) {
                continue;
            }
            int softp = top_soft(del.soft);
            if (softp != 0) {
                if (conf_.y) {
                    std::println(stderr, "\n\nWorking DEL 3' {} mate cluster cnt: {}", softp, del.vars_count);
                }
                Sclip* scv = d_.soft_clips_5end.get(softp);
                if (scv == nullptr) {
                    continue;
                }
                if (scv->used) {
                    continue;
                }
                const std::string seq = conseq(*scv);
                if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                    continue;
                }
                if (!loaded(del.mstart, del.mend)) {
                    reference(del.mstart, del.mend, 300);
                    partial(del.mstart - 200, del.mend + 200);
                }
                Match m = match(seq, softp, -1);
                int bp = m.base_position;
                if (bp == 0) {
                    m = match(seq, softp, -1, Config::kSeed2, 0);
                    bp = m.base_position;
                }
                if (bp == 0) {
                    continue;
                }
                if (!(softp - bp > 30 && is_overlap(bp, softp, del.mend, del.start, max_read_length_))) {
                    continue;
                }
                bp++;
                softp--;
                const int dellen = softp - bp + 1;
                Variation& variation = get_variation(niv, bp, "-" + i2s(dellen));
                variation.vars_count = 0;

                SvInfo& sv = get_sv(niv, bp);
                sv.type = "DEL";
                sv.pairs += del.vars_count;
                sv.splits += scv->vars_count;
                sv.clusters++;

                adj_cnt(variation, *scv);
                if (!(cov.contains(bp) && cov_at(cov, bp) > del.vars_count)) {
                    cov.put(bp, del.vars_count);
                }
                if (cov.contains(softp) && cov_at(cov, softp) > cov_at(cov, bp)) {
                    cov.put(bp, cov_at(cov, softp));
                }
                const int mcnt = del.vars_count;
                adj_cnt(variation, pair_support(mcnt, del));
                del.used = true;
                (void)mark_sv(bp, softp, {&d_.sv_structures.svfdel}, max_read_length_, conf_);
                if (conf_.y) {
                    std::println(stderr, "    Found DEL SV from 3' softclip unhappy reads: {} -+{} Cnt: {} AdjCnt: {}",
                                 bp, dellen, del.vars_count, mcnt);
                }
            } else {
                if (conf_.y) {
                    std::println(stderr, "\n\nWorking DEL 3' no softp mate cluster {} {} {} cnt: {}", region_.chr,
                                 del.mstart, del.mend, del.vars_count);
                }
                if (!loaded(del.mstart, del.mend)) {
                    reference(del.mstart, del.mend, 300);
                }
                bool done = false;
                std::vector<int> keys = d_.soft_clips_5end.keys();
                for (int i : keys) {
                    if (done) {
                        break;
                    }
                    Sclip& scv = *d_.soft_clips_5end.get(i);
                    if (scv.used) {
                        continue;
                    }
                    if (!(i <= del.start + 3 && del.start - i < 3 * max_read_length_)) {
                        continue;
                    }
                    const std::string seq = conseq(scv);
                    if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                        continue;
                    }
                    softp = i;
                    Match m = match(seq, softp, -1);
                    int bp = m.base_position;
                    if (bp == 0) {
                        m = match(seq, softp, -1, Config::kSeed2, 0);
                        bp = m.base_position;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    if (!(softp - bp > 30 && is_overlap(bp, softp, del.mend, del.start, max_read_length_))) {
                        continue;
                    }
                    bp++;
                    softp--;
                    const int dellen = softp - bp + 1;
                    Variation& variation = get_variation(niv, bp, "-" + i2s(dellen));
                    variation.vars_count = 0;

                    SvInfo& sv = get_sv(niv, bp);
                    sv.type = "DEL";
                    sv.pairs += del.vars_count;
                    sv.splits += scv.vars_count;
                    sv.clusters++;

                    adj_cnt(variation, scv);
                    if (!cov.contains(bp)) {
                        cov.put(bp, del.vars_count);
                    }
                    if (cov.contains(softp) && cov_at(cov, softp) > cov_at(cov, bp)) {
                        cov.put(bp, cov_at(cov, softp));
                    }
                    inc_cnt(cov, bp, scv.vars_count);
                    const int mcnt = del.vars_count;
                    adj_cnt(variation, pair_support(mcnt, del));
                    del.used = true;
                    (void)mark_sv(bp, softp, {&d_.sv_structures.svfdel}, max_read_length_, conf_);
                    if (conf_.y) {
                        std::println(stderr, "    Found DEL SV from 3' softclip happy reads: {} -{} Cnt: {} AdjCnt: {}",
                                     bp, dellen, del.vars_count, mcnt);
                    }
                    done = true;
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
}

void SvProcessor::find_inv() {
    find_inv_sub(d_.sv_structures.svfinv5, 1, Side::five);
    find_inv_sub(d_.sv_structures.svrinv5, -1, Side::five);
    find_inv_sub(d_.sv_structures.svfinv3, 1, Side::three);
    find_inv_sub(d_.sv_structures.svrinv3, -1, Side::three);
}

void SvProcessor::find_inv_sub(std::vector<Sclip>& svref, int dir, Side side) {
    PositionMap& niv = d_.non_insertion_variants;
    const RefMap& ref = ref_;
    int last_start = 0;
    // dir = 1 means 3' soft clipping
    for (Sclip& inv : svref) {
        try {
            last_start = inv.start;
            if (inv.used) {
                continue;
            }
            if (inv.vars_count < conf_.minr) {
                continue;
            }
            int softp = top_soft(inv.soft);
            auto& sclip = dir == 1 ? d_.soft_clips_3end : d_.soft_clips_5end;
            if (conf_.y) {
                std::println(stderr, "\n\nWorking INV {} {} {} pair_cnt: {}", softp, dir,
                             side == Side::five ? "_5" : "_3", inv.vars_count);
            }
            if (!loaded(inv.mstart, inv.mend)) {
                reference(inv.mstart, inv.mend, 500);
                partial(inv.mstart - 200, inv.mend + 200);
            }
            int bp = 0;
            Sclip empty;
            Sclip* scv = &empty;
            std::string seq;
            std::string extra;
            if (softp != 0) {
                Sclip* found = sclip.get(softp);
                if (found == nullptr) {
                    continue;
                }
                scv = found;
                if (scv->used) {
                    continue;
                }
                seq = conseq(*scv);
                if (seq.empty()) {
                    continue;
                }
                Match m = find_match_rev(seq, softp, dir);
                bp = m.base_position;
                extra = m.matched_sequence;
                if (bp == 0) {
                    m = find_match_rev(seq, softp, dir, Config::kSeed2, 0);
                    bp = m.base_position;
                    extra = m.matched_sequence;
                }
                if (bp == 0) {
                    continue;
                }
            } else {
                // Look within 100bp to see whether a soft cliping can be found but not associated
                // with discordant pairs
                const int sp = dir == 1 ? inv.end : inv.start;  // starting position
                for (int i = 1; i <= 2 * max_read_length_; i++) {
                    const int cp = sp + i * dir;
                    Sclip* found = sclip.get(cp);
                    if (found == nullptr) {
                        continue;
                    }
                    scv = found;
                    if (scv->used) {
                        continue;
                    }
                    seq = conseq(*scv);
                    if (seq.empty()) {
                        continue;
                    }
                    Match m = find_match_rev(seq, cp, dir);
                    bp = m.base_position;
                    extra = m.matched_sequence;
                    if (bp == 0) {
                        m = find_match_rev(seq, cp, dir, Config::kSeed2, 0);
                        bp = m.base_position;
                        extra = m.matched_sequence;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    softp = cp;
                    if ((dir == 1 && std::abs(bp - inv.mend) < Config::kMinSvcDist * max_read_length_) ||
                        (dir == -1 && std::abs(bp - inv.mstart) < Config::kMinSvcDist * max_read_length_)) {
                        break;
                    }
                }
                if (bp == 0) {
                    continue;
                }
            }
            if (conf_.y) {
                std::println(stderr, "    {} {} {} {} {} pair_cnt: {} soft_cnt: {}", softp, bp, dir,
                             side == Side::five ? "_5" : "_3", seq, inv.vars_count, scv->vars_count);
            }
            if (side == Side::five) {
                if (dir == -1) {
                    bp--;
                }
            } else {
                if (dir == 1) {
                    bp++;
                    if (bp != 0) {
                        softp--;
                    }
                } else {
                    softp--;
                }
            }
            if (side == Side::three) {
                std::swap(bp, softp);
            }
            if ((dir == -1 && side == Side::five) || (dir == 1 && side == Side::three)) {
                while (ref.contains(softp) && ref.contains(bp) && ref.at(softp) == util::complement(ref.at(bp))) {
                    softp++;
                    if (softp != 0) {
                        bp--;
                    }
                }
            }
            while (ref.contains(softp - 1) && ref.contains(bp + 1) &&
                   ref.at(softp - 1) == util::complement(ref.at(bp + 1))) {
                softp--;
                if (softp != 0) {
                    bp++;
                }
            }
            if (bp > softp && bp - softp > 150 && (bp - softp) / static_cast<double>(std::abs(inv.mlen)) < 1.5) {
                const int len = bp - softp + 1;
                const std::string ins5 = reverse_complement(join_ref(ref, bp - Config::kSvFlank + 1, bp));
                const std::string ins3 = reverse_complement(join_ref(ref, softp, softp + Config::kSvFlank - 1));
                std::string ins = ins5 + "<inv" + i2s(len - 2 * Config::kSvFlank) + ">" + ins3;
                if (len - 2 * Config::kSvFlank <= 0) {
                    ins = reverse_complement(join_ref(ref, softp, bp));
                }
                if (dir == 1 && !extra.empty()) {
                    extra = reverse_complement(extra);
                    ins = extra + ins;
                } else if (dir == -1 && !extra.empty()) {
                    ins = ins + extra;
                }
                const std::string gt = "-" + i2s(len) + "^" + ins;

                Variation& vref = get_variation(niv, softp, gt);
                inv.used = true;
                vref.pstd = true;
                vref.qstd = true;

                SvInfo& sv = get_sv(niv, softp);
                sv.type = "INV";
                sv.splits += scv->vars_count;
                sv.pairs += inv.vars_count;
                sv.clusters++;

                Variation* vref_softp = nullptr;
                if (dir == -1) {
                    if (VariationMap* at = niv.get(softp); at != nullptr && ref.contains(softp)) {
                        vref_softp = at->get(std::string(1, ref.at(softp)));
                    }
                }
                adj_cnt(vref, *scv, vref_softp);
                CountMap dels5;
                dels5.get_or_create(softp).put(gt, inv.vars_count);
                const int* prev = d_.ref_coverage.get(softp - 1);
                d_.ref_coverage.put(softp, prev ? *prev : inv.vars_count);
                scv->used = true;

                realigndel_with_bams(d_, scope_, dels5);

                if (conf_.y) {
                    const int* c = d_.ref_coverage.get(softp);
                    std::println(stderr,
                                 "  Found INV SV: {} {} {} BP: {} cov: {} Cnt: {} EXTRA: {} {} {} {} cnt: {} {}\t DIR: "
                                 "{} Side: {}",
                                 seq, softp, gt, bp, c ? i2s(*c) : "null", inv.vars_count, extra, inv.mstart, inv.mend,
                                 inv.mlen, scv->vars_count,
                                 inv.mlen != 0 ? i2s((bp - softp) / std::abs(inv.mlen)) : "?", dir,
                                 side == Side::five ? "_5" : "_3");
                }
                return;
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
}

void SvProcessor::findsv() {
    PositionMap& niv = d_.non_insertion_variants;
    const RefMap& ref = ref_;
    auto& cov = d_.ref_coverage;
    struct SortPositionSclip {
        int position;
        Sclip* soft_clip;
        int count;
    };
    // fillAndSortTmpSV: unused clips in CURSEG (the
    // region), count desc, a stable sort of HashMap order.
    const auto fill = [&](util::HashOrderMap<int, Sclip>& clips) {
        std::vector<SortPositionSclip> tmp;
        clips.for_each([&](int position, Sclip& sc) {
            if (sc.used) {
                return;
            }
            if (position < region_.start || position > region_.end) {
                return;
            }
            tmp.push_back({position, &sc, sc.vars_count});
        });
        std::stable_sort(tmp.begin(), tmp.end(),
                         [](const SortPositionSclip& a, const SortPositionSclip& b) { return a.count > b.count; });
        return tmp;
    };
    const auto first_sv_used = [&](int p) {
        const std::vector<Sclip*>* v = d_.softp2sv.get(p);
        if (v == nullptr) {
            return false;
        }
        if (v->empty()) {
            throw util::RecoverableError("IndexOutOfBoundsException: Index 0 out of bounds for length 0");
        }
        return (*v)[0]->used;
    };
    const auto pairs_support = [](const PairsData& pd) {
        Variation tmp;
        tmp.vars_count = pd.pairs;
        tmp.high_quality_reads_count = pd.pairs;
        tmp.vars_count_on_forward = pd.pairs / 2;
        tmp.vars_count_on_reverse = pd.pairs - pd.pairs / 2;
        tmp.mean_position = pd.pmean;
        tmp.mean_quality = pd.qmean;
        tmp.mean_mapping_quality = pd.big_qmean;
        tmp.number_of_mismatches = pd.nm;
        return tmp;
    };
    const std::vector<std::vector<Sclip>*> dels = {&d_.sv_structures.svfdel, &d_.sv_structures.svrdel};

    auto tmp5 = fill(d_.soft_clips_5end);
    auto tmp3 = fill(d_.soft_clips_3end);
    int last_position = 0;
    for (const SortPositionSclip& tuple5 : tmp5) {
        try {
            int p5 = tuple5.position;
            last_position = p5;
            Sclip& sc5v = *tuple5.soft_clip;
            const int cnt5 = tuple5.count;
            if (cnt5 < conf_.minr) {
                break;
            }
            if (sc5v.used) {
                continue;
            }
            if (first_sv_used(p5)) {
                continue;
            }
            const std::string seq = conseq(sc5v);
            if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Finding SV 5': {} {} cnt: {}", seq, p5, cnt5);
            }
            Match m = match(seq, p5, -1);
            int bp = m.base_position;
            std::string EXTRA = m.matched_sequence;
            if (bp != 0) {
                // candidate deletion
                if (bp < p5) {
                    const PairsData pd = check_pairs(bp, p5, dels);
                    if (pd.pairs == 0) {
                        continue;
                    }
                    p5--;
                    bp++;
                    const int dellen = p5 - bp + 1;
                    Variation& vref = get_variation(niv, bp, "-" + i2s(dellen));
                    vref.vars_count = 0;

                    SvInfo& sv = get_sv(niv, bp);
                    sv.type = "DEL";
                    sv.pairs += pd.pairs;
                    sv.splits += cnt5;
                    sv.clusters += pd.pairs != 0 ? 1 : 0;

                    if (!cov.contains(bp)) {
                        cov.put(bp, pd.pairs + sc5v.vars_count);
                    }
                    if (cov.contains(p5 + 1) && cov_at(cov, bp) < cov_at(cov, p5 + 1)) {
                        cov.put(bp, cov_at(cov, p5 + 1));
                    }
                    adj_cnt(vref, sc5v);
                    adj_cnt(vref, pairs_support(pd));
                    if (conf_.y) {
                        std::println(stderr, "    Finding candidate deletion 5'");
                    }
                }  // else: candidate duplication (not implemented upstream)
            } else {  // candidate inversion
                Match mr = find_match_rev(seq, p5, -1);
                bp = mr.base_position;
                EXTRA = mr.matched_sequence;
                if (bp == 0) {
                    continue;
                }
                if (!(std::abs(bp - p5) > Config::kSvFlank)) {
                    continue;
                }
                if (bp <= p5) {  // bp at 5' side
                    std::swap(bp, p5);
                }
                bp--;
                while (ref.contains(bp + 1) && is_has_and_equals(util::complement(ref.at(bp + 1)), ref, p5 - 1)) {
                    p5--;
                    if (p5 != 0) {
                        bp++;
                    }
                }
                const std::string ins5 = reverse_complement(join_ref(ref, bp - Config::kSvFlank + 1, bp));
                const std::string ins3 = reverse_complement(join_ref(ref, p5, p5 + Config::kSvFlank - 1));
                const int mid = bp - p5 - static_cast<int>(ins5.size()) - static_cast<int>(ins3.size()) + 1;
                std::string vn = "-" + i2s(bp - p5 + 1) + "^" + ins5 + "<inv" + i2s(mid) + ">" + ins3 + EXTRA;
                if (mid <= 0) {
                    const std::string tins = reverse_complement(join_ref(ref, p5, bp));
                    vn = "-" + i2s(bp - p5 + 1) + "^" + tins + EXTRA;
                }
                Variation& vref = get_variation(niv, p5, vn);

                SvInfo& sv = get_sv(niv, p5);
                sv.type = "INV";
                sv.splits += cnt5;

                adj_cnt(vref, sc5v);
                inc_cnt(cov, p5, cnt5);
                if (cov.contains(bp) && cov_at(cov, p5) < cov_at(cov, bp)) {
                    cov.put(p5, cov_at(cov, bp));
                }
                if (conf_.y) {
                    std::println(stderr, "    Found INV: {} {} Cnt:{}", p5, vn, cnt5);
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
    for (const SortPositionSclip& tuple3 : tmp3) {
        try {
            int p3 = tuple3.position;
            last_position = p3;
            Sclip& sc3v = *tuple3.soft_clip;
            const int cnt3 = tuple3.count;
            if (cnt3 < conf_.minr) {
                break;
            }
            if (sc3v.used) {
                continue;
            }
            if (first_sv_used(p3)) {
                continue;
            }
            const std::string seq = conseq(sc3v);
            if (seq.empty() || static_cast<int>(seq.size()) < Config::kSeed2) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Finding SV 3': {} {} cnt: {}", seq, p3, cnt3);
            }
            Match m = match(seq, p3, 1);
            int bp = m.base_position;
            std::string EXTRA = m.matched_sequence;
            if (bp != 0) {
                if (bp > p3) {  // candidate deletion
                    const PairsData pd = check_pairs(p3, bp, dels);
                    if (pd.pairs == 0) {
                        continue;
                    }
                    const int dellen = bp - p3;
                    bp--;
                    while (ref_equal(ref, bp, p3 - 1)) {
                        bp--;
                        if (bp != 0) {
                            p3--;
                        }
                    }
                    Variation& vref = get_variation(niv, p3, "-" + i2s(dellen));
                    vref.vars_count = 0;

                    SvInfo& sv = get_sv(niv, p3);
                    sv.type = "DEL";
                    sv.pairs += pd.pairs;
                    sv.splits += cnt3;
                    sv.clusters += pd.pairs != 0 ? 1 : 0;

                    if (!cov.contains(p3)) {
                        cov.put(p3, pd.pairs + sc3v.vars_count);
                    }
                    if (cov.contains(bp) && cov_at(cov, bp) < cov_at(cov, p3)) {
                        cov.put(bp, cov_at(cov, p3));
                    }
                    adj_cnt(vref, sc3v);
                    adj_cnt(vref, pairs_support(pd));
                    if (conf_.y) {
                        std::println(stderr, "    Finding candidate deletion 3'");
                    }
                }  // else: candidate duplication (not implemented upstream)
            } else {  // candidate inversion
                Match mr = find_match_rev(seq, p3, 1);
                bp = mr.base_position;
                EXTRA = mr.matched_sequence;
                if (bp == 0) {
                    continue;
                }
                if (std::abs(bp - p3) <= Config::kSvFlank) {
                    continue;
                }
                if (bp < p3) {  // bp at 5' side
                    std::swap(bp, p3);
                    p3++;
                    bp--;
                }
                while (ref.contains(bp + 1) && is_has_and_equals(util::complement(ref.at(bp + 1)), ref, p3 - 1)) {
                    p3--;
                    if (p3 != 0) {
                        bp++;
                    }
                }
                const std::string ins5 = reverse_complement(join_ref(ref, bp - Config::kSvFlank + 1, bp));
                const std::string ins3 = reverse_complement(join_ref(ref, p3, p3 + Config::kSvFlank - 1));
                const int mid = bp - p3 - 2 * Config::kSvFlank + 1;
                std::string vn = "-" + i2s(bp - p3 + 1) + "^" + EXTRA + ins5 + "<inv" + i2s(mid) + ">" + ins3;
                if (mid <= 0) {
                    const std::string tins = reverse_complement(join_ref(ref, p3, bp));
                    vn = "-" + i2s(bp - p3 + 1) + "^" + EXTRA + tins;
                }
                Variation& vref = get_variation(niv, p3, vn);

                SvInfo& sv = get_sv(niv, p3);
                sv.type = "INV";
                sv.splits += cnt3;

                adj_cnt(vref, sc3v);
                inc_cnt(cov, p3, cnt3);
                if (cov.contains(bp) && cov_at(cov, p3) < cov_at(cov, bp)) {
                    cov.put(p3, cov_at(cov, bp));
                }
                if (conf_.y) {
                    std::println(stderr, "    Found INV: {} BP: {} {} EXTRA: {} Cnt: {}", p3, bp, vn, EXTRA, cnt3);
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
}

void SvProcessor::find_del_disc() {
    PositionMap& niv = d_.non_insertion_variants;
    auto& cov = d_.ref_coverage;
    const int MINDIST = 8 * max_read_length_;  // the minimum distance between two clusters
    int last_start = 0;
    const auto support = [](const Sclip& del) {
        Variation tv;
        tv.vars_count = 2 * del.vars_count;
        tv.high_quality_reads_count = 2 * del.vars_count;
        tv.vars_count_on_forward = del.vars_count;
        tv.vars_count_on_reverse = del.vars_count;
        tv.mean_quality = 2 * del.mean_quality;
        tv.mean_position = 2 * del.mean_position;
        tv.mean_mapping_quality = 2 * del.mean_mapping_quality;
        tv.number_of_mismatches = 2 * del.number_of_mismatches;
        return tv;
    };
    const auto sclip_count = [](util::HashOrderMap<int, Sclip>& m, int p) {
        const Sclip* s = m.get(p);
        return s ? s->vars_count : 0;
    };
    for (Sclip& del : d_.sv_structures.svfdel) {
        try {
            last_start = del.start;
            if (del.used) {
                continue;
            }
            if (!scope_.splice.empty() && std::abs(del.mlen) < 250000) {
                continue;  // more stringent for RNA-Seq
            }
            if (del.vars_count < conf_.minr + 5) {
                continue;
            }
            if (del.mstart <= del.end + MINDIST) {
                continue;
            }
            if (del.mean_mapping_quality / del.vars_count <= Config::kDiscPairQual) {
                continue;
            }
            const int mlen = del.mstart - del.end - max_read_length_ / (del.vars_count + 1);
            if (!(mlen > 0 && mlen > MINDIST)) {
                continue;
            }
            int bp = del.end + (max_read_length_ / (del.vars_count + 1)) / 2;
            if (del.softp != 0) {
                bp = del.softp;
            }
            if (!ref_.contains(bp)) {
                reference(bp - 150, bp + 150, mlen < 1000 ? mlen : 1000);
            }
            Variation& vref = get_variation(niv, bp, "-" + i2s(mlen));
            vref.vars_count = 0;
            SvInfo& sv = get_sv(niv, bp);
            sv.type = "DEL";
            sv.splits += sclip_count(d_.soft_clips_3end, del.end + 1);
            sv.splits += sclip_count(d_.soft_clips_5end, del.mstart);
            sv.pairs += del.vars_count;
            sv.clusters++;
            if (conf_.y) {
                std::println(stderr, "  Found DEL with discordant pairs only: cnt: {} BP: {} Len: {} {}-{}<->{}-{}",
                             del.vars_count, bp, mlen, del.start, del.end, del.mstart, del.mend);
            }
            adj_cnt(vref, support(del));
            if (!cov.contains(bp)) {
                cov.put(bp, 2 * del.vars_count);
            }
            del.used = true;
            (void)mark_sv(del.end, del.mstart, {&d_.sv_structures.svrdel}, max_read_length_, conf_);
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
    for (Sclip& del : d_.sv_structures.svrdel) {
        try {
            last_start = del.start;
            if (del.used) {
                continue;
            }
            if (!scope_.splice.empty() && std::abs(del.mlen) < 250000) {
                continue;  // more stringent for RNA-Seq
            }
            if (del.vars_count < conf_.minr + 5) {
                continue;
            }
            if (del.start <= del.mend + MINDIST) {
                continue;
            }
            if (del.mean_mapping_quality / del.vars_count <= Config::kDiscPairQual) {
                continue;
            }
            const int mlen = del.start - del.mend - max_read_length_ / (del.vars_count + 1);
            if (!(mlen > 0 && mlen > MINDIST)) {
                continue;
            }
            const int bp = del.mend + ((max_read_length_ / (del.vars_count + 1)) / 2);
            if (!ref_.contains(bp)) {
                reference(bp - 150, bp + 150, mlen < 1000 ? mlen : 1000);
            }
            Variation& vref = get_variation(niv, bp, "-" + i2s(mlen));
            vref.vars_count = 0;

            SvInfo& sv = get_sv(niv, bp);
            sv.type = "DEL";
            sv.splits += sclip_count(d_.soft_clips_3end, del.mend + 1);
            sv.splits += sclip_count(d_.soft_clips_5end, del.start);
            sv.pairs += del.vars_count;
            sv.clusters += 1;
            if (conf_.y) {
                std::println(stderr,
                             "  Found DEL with discordant pairs only (reverse): cnt: {} BP: {} Len: {} {}-{}<->{}-{}",
                             del.vars_count, bp, mlen, del.start, del.end, del.mstart, del.mend);
            }
            if (del.softp != 0) {
                if (Sclip* s = d_.soft_clips_5end.get(del.softp)) {
                    s->used = true;
                }
            }
            adj_cnt(vref, support(del));
            if (!cov.contains(bp)) {
                cov.put(bp, 2 * del.vars_count);
            }
            if (cov.contains(del.start) && cov_at(cov, bp) < cov_at(cov, del.start)) {
                cov.put(bp, cov_at(cov, del.start));
            }
            del.used = true;
            reference(del.mstart - 100, del.mend + 100, 200);
            (void)mark_sv(del.mend, del.start, {&d_.sv_structures.svfdel}, max_read_length_, conf_);
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
}

void SvProcessor::find_inv_disc() {
    PositionMap& niv = d_.non_insertion_variants;
    const RefMap& ref = ref_;
    int last_start = 0;
    const auto found = [&](int bp, int pe, int cnt, int rcnt, const Sclip& f, const Sclip& r) -> Variation& {
        const int len = pe - bp + 1;
        const std::string ins5 = reverse_complement(join_ref(ref, bp, bp + Config::kSvFlank - 1));
        const std::string ins3 = reverse_complement(join_ref(ref, pe - Config::kSvFlank + 1, pe));
        std::string ins = ins3 + "<inv" + i2s(len - 2 * Config::kSvFlank) + ">" + ins5;
        if (len - 2 * Config::kSvFlank <= 0) {
            ins = reverse_complement(join_ref(ref, bp, pe));
        }
        if (conf_.y) {
            std::println(stderr, "  Found INV with discordant pairs only: cnt: {} Len: {} {}", cnt, len, ins);
        }
        Variation& vref = get_variation(niv, bp, "-" + i2s(len) + "^" + ins);
        vref.pstd = true;
        vref.qstd = true;
        Variation tmp;
        tmp.vars_count = cnt + rcnt;
        tmp.high_quality_reads_count = cnt + rcnt;
        tmp.vars_count_on_forward = cnt;
        tmp.vars_count_on_reverse = rcnt;
        tmp.mean_quality = f.mean_quality + r.mean_quality;
        tmp.mean_position = f.mean_position + r.mean_position;
        tmp.mean_mapping_quality = f.mean_mapping_quality + r.mean_mapping_quality;
        tmp.number_of_mismatches = f.number_of_mismatches + r.number_of_mismatches;
        adj_cnt(vref, tmp);
        return vref;
    };
    const auto sclip_count = [](util::HashOrderMap<int, Sclip>& m, int p) {
        const Sclip* s = m.get(p);
        return s ? s->vars_count : 0;
    };
    for (Sclip& invf5 : d_.sv_structures.svfinv5) {
        try {
            if (invf5.used) {
                continue;
            }
            // Values read once, before the inner loop (locals upstream).
            const Sclip f = invf5;
            const int cnt = f.vars_count;
            const int me = f.mend;
            const int ms = f.mstart;
            const int end = f.end;
            const int start = f.start;
            last_start = start;
            if (!(f.mean_mapping_quality / static_cast<double>(cnt) > Config::kDiscPairQual)) {
                continue;
            }
            for (Sclip& invr5 : d_.sv_structures.svrinv5) {
                try {
                    if (invr5.used) {
                        continue;
                    }
                    const int rcnt = invr5.vars_count;
                    const int rstart = invr5.start;
                    last_start = rstart;
                    const int rms = invr5.mstart;
                    if (!(invr5.mean_mapping_quality / static_cast<double>(rcnt) > Config::kDiscPairQual)) {
                        continue;
                    }
                    if (!(cnt + rcnt > conf_.minr + 5)) {
                        continue;
                    }
                    if (is_overlap(end, me, rstart, rms, max_read_length_)) {
                        const int bp = std::abs((end + rstart) / 2);
                        const int pe = std::abs((me + rms) / 2);
                        if (!ref.contains(pe)) {
                            reference(pe - 150, pe + 150, 300);
                        }
                        found(bp, pe, cnt, rcnt, f, invr5);
                        invf5.used = true;
                        invr5.used = true;

                        SvInfo& sv = get_sv(niv, bp);
                        sv.type = "INV";
                        sv.pairs += cnt;
                        sv.splits += sclip_count(d_.soft_clips_5end, start);
                        sv.splits += sclip_count(d_.soft_clips_5end, ms);
                        sv.clusters++;
                        if (!d_.ref_coverage.contains(bp)) {
                            d_.ref_coverage.put(bp, 2 * cnt);
                        }
                        (void)mark_sv(bp, pe, {&d_.sv_structures.svfinv3, &d_.sv_structures.svrinv3}, max_read_length_,
                                      conf_);
                    }
                } catch (const util::RecoverableError& e) {
                    exception(e, last_start);
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
    for (Sclip& invf3 : d_.sv_structures.svfinv3) {
        try {
            if (invf3.used) {
                continue;
            }
            const Sclip f = invf3;
            const int cnt = f.vars_count;
            const int me = f.mend;
            const int end = f.end;
            for (Sclip& invr3 : d_.sv_structures.svrinv3) {
                try {
                    if (invr3.used) {
                        continue;
                    }
                    const int rcnt = invr3.vars_count;
                    const int rstart = invr3.start;
                    last_start = rstart;
                    const int rms = invr3.mstart;
                    if (!(invr3.mean_mapping_quality / static_cast<double>(rcnt) > Config::kDiscPairQual)) {
                        continue;
                    }
                    if (!(cnt + rcnt > conf_.minr + 5)) {
                        continue;
                    }
                    if (is_overlap(me, end, rms, rstart, max_read_length_)) {
                        const int pe = std::abs((end + rstart) / 2);
                        const int bp = std::abs((me + rms) / 2);
                        if (!ref.contains(bp)) {
                            reference(bp - 150, bp + 150, 300);
                        }
                        found(bp, pe, cnt, rcnt, f, invr3);
                        invf3.used = true;
                        invr3.used = true;

                        SvInfo& sv = get_sv(niv, bp);
                        sv.type = "INV";
                        sv.pairs += cnt;
                        sv.splits += sclip_count(d_.soft_clips_3end, end + 1);
                        sv.splits += sclip_count(d_.soft_clips_3end, me + 1);
                        sv.clusters++;
                        if (!d_.ref_coverage.contains(bp)) {
                            d_.ref_coverage.put(bp, 2 * cnt);
                        }
                        (void)mark_sv(bp, pe, {&d_.sv_structures.svfinv5, &d_.sv_structures.svrinv5}, max_read_length_,
                                      conf_);
                    }
                } catch (const util::RecoverableError& e) {
                    exception(e, last_start);
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
}

void SvProcessor::find_dup_disc() {
    PositionMap& niv = d_.non_insertion_variants;
    PositionMap& iv = d_.insertion_variants;
    const RefMap& ref = ref_;
    auto& cov = d_.ref_coverage;
    int last_start = 0;
    struct Means {
        int cnt;
        double qmean, big_qmean, pmean, nm;
    };
    const auto means_of = [](const Variation& v) {
        return Means{v.vars_count, v.mean_quality, v.mean_mapping_quality, v.mean_position, v.number_of_mismatches};
    };
    const auto add_support = [](Variation& vref, const Means& f, const Means& r) {
        const int tcnt = r.cnt + f.cnt;
        Variation tmp;
        tmp.vars_count = tcnt;
        tmp.extracnt = tcnt;
        tmp.high_quality_reads_count = tcnt;
        tmp.vars_count_on_forward = f.cnt;
        tmp.vars_count_on_reverse = r.cnt;
        tmp.mean_quality = f.qmean + r.qmean;
        tmp.mean_position = f.pmean + r.pmean;
        tmp.mean_mapping_quality = f.big_qmean + r.big_qmean;
        tmp.number_of_mismatches = f.nm + r.nm;
        adj_cnt(vref, tmp);
        return tcnt;
    };
    for (Sclip& dup : d_.sv_structures.svfdup) {
        try {
            if (dup.used) {
                continue;
            }
            const int ms = dup.mstart;
            const int me = dup.mend;
            const int cnt = dup.vars_count;
            int end = dup.end;
            const int start = dup.start;
            last_start = start;
            if (!(cnt >= conf_.minr + 5)) {
                continue;
            }
            if (!(dup.mean_mapping_quality / cnt > Config::kDiscPairQual)) {
                continue;
            }
            int mlen = end - ms + max_read_length_ / cnt;
            int bp = ms - (max_read_length_ / cnt) / 2;
            int pe = end;
            if (!loaded(ms, me)) {
                if (!ref.contains(bp)) {
                    reference(bp - 150, bp + 150, 300);
                }
                partial(ms - 200, me + 200);
            }
            Means f = means_of(dup);
            Means r = f;
            if (!dup.soft.empty()) {
                pe = top_soft(dup.soft);
                Sclip* sc3 = d_.soft_clips_3end.get(pe);
                if (sc3 == nullptr) {
                    continue;
                }
                if (sc3->used) {
                    continue;
                }
                f = means_of(*sc3);
                const std::string seq = conseq(*sc3);
                const Match m = match(seq, bp, 1);
                int tbp = m.base_position;
                if (tbp != 0 && tbp < pe) {
                    sc3->used = true;
                    while (ref_equal(ref, pe - 1, tbp - 1)) {
                        tbp--;
                        if (tbp != 0) {
                            pe--;
                        }
                    }
                    mlen = pe - tbp;
                    bp = tbp;
                    pe--;
                    end = pe;
                    if (const Sclip* sc5 = d_.soft_clips_5end.get(bp)) {
                        r = means_of(*sc5);
                    }
                }
            }
            const std::string ins5 = join_ref(ref, bp, bp + Config::kSvFlank - 1);
            const std::string ins3 = join_ref(ref, pe - Config::kSvFlank + 1, pe);
            const std::string ins = ins5 + "<dup" + i2s(mlen - 2 * Config::kSvFlank) + ">" + ins3;

            Variation& vref = get_variation(iv, bp, "+" + ins);
            vref.vars_count = 0;

            SvInfo& sv = get_sv(niv, bp);
            sv.type = "DUP";
            sv.pairs += cnt;
            const Sclip* soft_sc = dup.softp != 0 ? d_.soft_clips_3end.get(dup.softp) : nullptr;
            sv.splits += soft_sc ? soft_sc->vars_count : 0;
            sv.clusters++;
            if (conf_.y) {
                std::println(stderr,
                             "  Found DUP with discordant pairs only (forward): cnt: {} BP: {} END: {} {} Len: {} "
                             "{}-{}<->{}-{}",
                             cnt, bp, pe, ins, mlen, start, end, ms, me);
            }
            const int tcnt = add_support(vref, f, r);
            dup.used = true;
            if (!cov.contains(bp)) {
                cov.put(bp, tcnt);
            }
            if (cov.contains(end) && cov_at(cov, bp) < cov_at(cov, end)) {
                cov.put(bp, cov_at(cov, end));
            }
            const SvMark mk = mark_dup_sv(bp, pe, {&d_.sv_structures.svrdup}, max_read_length_, conf_);
            sv.clusters += mk.count;
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
    for (Sclip& dup : d_.sv_structures.svrdup) {
        try {
            if (dup.used) {
                continue;
            }
            const int ms = dup.mstart;
            const int me = dup.mend;
            const int cnt = dup.vars_count;
            const int end = dup.end;
            const int start = dup.start;
            last_start = start;
            if (cnt < conf_.minr + 5) {
                continue;
            }
            if (!(dup.mean_mapping_quality / cnt > Config::kDiscPairQual)) {
                continue;
            }
            int mlen = me - start + max_read_length_ / cnt;
            int bp = start - (max_read_length_ / cnt) / 2;
            int pe = mlen + bp - 1;
            int tpe = pe;
            if (!loaded(ms, me)) {
                if (!ref.contains(pe)) {
                    reference(pe - 150, pe + 150, 300);
                }
                partial(ms - 200, me + 200);
            }
            Means f = means_of(dup);
            Means r = f;
            if (!dup.soft.empty()) {
                bp = top_soft(dup.soft);
                Sclip* sc5 = d_.soft_clips_5end.get(bp);
                if (sc5 == nullptr) {
                    continue;
                }
                if (sc5->used) {
                    continue;
                }
                r = means_of(*sc5);
                const std::string seq = conseq(*sc5);
                const Match m = match(seq, pe, -1);
                const int tbp = m.base_position;
                if (tbp != 0 && tbp > bp) {
                    sc5->used = true;
                    pe = tbp;
                    mlen = pe - bp + 1;
                    tpe = pe + 1;
                    while (ref_equal(ref, tpe, bp + (tpe - pe - 1))) {
                        tpe++;
                    }
                    if (const Sclip* sc3 = d_.soft_clips_3end.get(tpe)) {
                        f = means_of(*sc3);
                    }
                }
            }
            const std::string ins5 = join_ref(ref, bp, bp + Config::kSvFlank - 1);
            const std::string ins3 = join_ref(ref, pe - Config::kSvFlank + 1, pe);
            const std::string ins = ins5 + "<dup" + i2s(mlen - 2 * Config::kSvFlank) + ">" + ins3;

            Variation& vref = get_variation(iv, bp, "+" + ins);
            vref.vars_count = 0;

            SvInfo& sv = get_sv(niv, bp);
            sv.type = "DUP";
            sv.pairs += cnt;
            const Sclip* s5 = d_.soft_clips_5end.get(bp);
            sv.splits += s5 ? s5->vars_count : 0;
            const Sclip* s3 = d_.soft_clips_3end.get(tpe);
            sv.splits += s3 ? s3->vars_count : 0;
            sv.clusters++;
            if (conf_.y) {
                std::println(stderr,
                             "  Found DUP with discordant pairs only (reverse): cnt: {} BP: {} Len: {} {}-{}<->{}-{}",
                             cnt, bp, mlen, start, end, ms, me);
            }
            const int tcnt = add_support(vref, f, r);
            dup.used = true;
            if (!cov.contains(bp)) {
                cov.put(bp, tcnt);
            }
            if (cov.contains(me) && cov_at(cov, bp) < cov_at(cov, me)) {
                cov.put(bp, cov_at(cov, me));
            }
            const SvMark mk = mark_dup_sv(bp, pe, {&d_.sv_structures.svfdup}, max_read_length_, conf_);
            sv.clusters += mk.count;
        } catch (const util::RecoverableError& e) {
            exception(e, last_start);
        }
    }
}

PairsData SvProcessor::check_pairs(int start, int end, const std::vector<std::vector<Sclip>*>& sv) {
    int pairs = 0;
    PairsData pd;
    for (std::vector<Sclip>* cluster : sv) {
        for (Sclip& svr : *cluster) {
            if (svr.used) {
                continue;
            }
            int s = (svr.start + svr.end) / 2;
            int e = (svr.mstart + svr.mend) / 2;
            if (s > e) {
                std::swap(s, e);
            }
            if (!is_overlap(start, end, s, e, max_read_length_)) {
                continue;
            }
            if (svr.vars_count > pairs) {
                pd = PairsData{svr.vars_count, svr.mean_position, svr.mean_quality, svr.mean_mapping_quality,
                               svr.number_of_mismatches};
                pairs = svr.vars_count;
            }
            svr.used = true;
            if (conf_.y) {
                std::println(stderr, "      Pair [{}:{}-{}] overlapping [{}:{}-{}] found and marked.", region_.chr, s,
                             e, region_.chr, start, end);
            }
        }
    }
    return pd;
}

Match SvProcessor::find_match_rev(std::string seq, int position, int dir, int seed, int mm) const {
    const RefMap& ref = ref_;
    // dir = 1 means from 3' soft clip
    if (dir == 1) {
        seq = util::reverse(seq);
    }
    seq = util::complement(seq);
    if (conf_.y) {
        std::println(stderr, "    Working MatchRev {} {} {}", position, seq, dir);
    }
    std::string extra;
    const int seq_len = static_cast<int>(seq.size());
    const auto not_equal = [&](int c, int pos) {  // isHasAndNotEquals(Utils.charAt(...), ref, pos)
        auto r = ref_get(ref, pos);
        return r && static_cast<unsigned char>(*r) != c;
    };
    for (int i = seq_len - seed; i >= 0; i--) {
        const std::string s = util::substr(seq, i, seed);
        auto it = scope_.reference.seed.find(s);
        if (it == scope_.reference.seed.end() || it->second.size() != 1) {
            continue;
        }
        const int first_seed = it->second[0];
        int bp = dir == 1 ? first_seed + seq_len - i - 1 : first_seed - i;
        if (ismatchref(seq, ref, bp, -1 * dir, mm)) {
            if (conf_.y) {
                std::println(stderr, "      Found SV BP (reverse): {} BP: {} SEEDpos: {} {} {} {} {} ", dir, bp,
                             first_seed, position, s, i, seq);
            }
            return {bp, extra};
        }
        // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
        std::string sseq = seq;
        int eqcnt = 0;
        for (int j = 1; j <= 15; j++) {
            bp -= dir;
            sseq = dir == -1 ? util::substr(sseq, 1) : util::substr(sseq, 0, -1);
            if (dir == -1) {
                if (not_equal(util::utils_char_at(sseq, 0), bp)) {
                    continue;
                }
                eqcnt++;
                if (not_equal(util::utils_char_at(sseq, 1), bp + 1)) {
                    continue;
                }
                extra = util::substr(seq, 0, j);
            } else {
                if (not_equal(util::utils_char_at(sseq, -1), bp)) {
                    continue;
                }
                eqcnt++;
                if (not_equal(util::utils_char_at(sseq, -2), bp - 1)) {
                    continue;
                }
                extra = util::substr(seq, -j);
            }
            if (eqcnt >= 3 && eqcnt / static_cast<double>(j) > 0.5) {
                break;
            }
            if (conf_.y) {
                std::println(stderr, "      FoundSEED SV BP (reverse): {} BP: {} SEEDpos{} {} {} {} {} EXTRA: {}", dir,
                             bp, first_seed, position, s, i, seq, extra);
            }
            if (ismatchref(sseq, ref, bp, -1 * dir, 1)) {
                return {bp, extra};
            }
        }
    }
    return {0, ""};
}

void SvProcessor::output_clipping() {
    std::println(stderr, "5' Remaining clipping reads");
    d_.soft_clips_5end.for_each([&](int p, Sclip& sc) {
        if (sc.used || sc.vars_count < conf_.minr) {
            return;
        }
        std::string seq = conseq(sc);
        if (!seq.empty() && static_cast<int>(seq.size()) > Config::kSeed2) {
            std::println(stderr, "  P: {} Cnt: {} Seq: {}", p, sc.vars_count, util::reverse(seq));
        }
    });
    std::println(stderr, "3' Remaining clipping reads");
    d_.soft_clips_3end.for_each([&](int p, Sclip& sc) {
        if (sc.used || sc.vars_count < conf_.minr) {
            return;
        }
        const std::string seq = conseq(sc);
        if (!seq.empty() && static_cast<int>(seq.size()) > Config::kSeed2) {
            std::println(stderr, "  P: {} Cnt: {} Seq: {}", p, sc.vars_count, seq);
        }
    });
}

}  // namespace

void sv_process(VariationData& data, RegionScope& scope) {
    SvProcessor sv(data, scope);
    if (!scope.conf.disable_sv) {
        sv.find_all_svs();
    }
    adj_snv(data, scope);
    if (scope.conf.y) {
        sv.output_clipping();
        std::println(stderr, "TIME: Finish realign");
    }
}

}  // namespace senri::core
