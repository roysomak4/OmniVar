#include "core/realigner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <print>
#include <set>

#include "core/cigar.hpp"
#include "core/cigar_parser.hpp"
#include "io/bam.hpp"
#include "io/record_preprocessor.hpp"
#include "io/reference.hpp"
#include "util/checked_string.hpp"
#include "util/num_format.hpp"
#include "util/number_parse.hpp"
#include "util/regex_engine.hpp"
#include "util/sequence.hpp"

namespace omnivar::core {

namespace {

struct SortPositionDescription {
    int position;
    std::string description_string;
    int count;
};

// VariationRealigner.fillAndSortTmp: count desc, position asc,
// description desc. A total order, so the HashMap order of `changes` does not matter.
std::vector<SortPositionDescription> fill_and_sort_tmp(const CountMap& changes) {
    std::vector<SortPositionDescription> tmp;
    changes.for_each([&](int position, const util::HashOrderMap<std::string, int>& v) {
        v.for_each([&](const std::string& desc, int cnt) { tmp.push_back({position, desc, cnt}); });
    });
    std::stable_sort(tmp.begin(), tmp.end(), [](const SortPositionDescription& a, const SortPositionDescription& b) {
        if (a.count != b.count) {
            return a.count > b.count;
        }
        if (a.position != b.position) {
            return a.position < b.position;
        }
        return a.description_string > b.description_string;
    });
    return tmp;
}

// StringBuilder(s).insert(1, "&")
std::string insert_amp(const std::string& s) { return s.substr(0, 1) + "&" + s.substr(1); }

// 16-bit char from Utils.charAt: (char) -1 when a negative index runs past the start.
constexpr int kNoChar = 0xFFFF;
// The same char appended to a String (only ever compared, never a base).
char as_char(int c) { return c == kNoChar ? '\xff' : static_cast<char>(c); }

// VariationUtils.isHasAndEquals / isHasAndNotEquals (Character, ref, index) for a 16-bit char.
bool has_and_equals(int ch, const RefMap& ref, int index) {
    auto r = ref_get(ref, index);
    return r && static_cast<unsigned char>(*r) == ch;
}
bool has_and_not_equals(int ch, const RefMap& ref, int index) {
    auto r = ref_get(ref, index);
    return r && static_cast<unsigned char>(*r) != ch;
}
// VariationUtils.isEquals(Character, Character) with a char on one side.
bool equals_char(std::optional<char> r, char c) { return r && *r == c; }

// String.replaceAll("#|\\^", "")
std::string strip_hash_caret(std::string_view s) {
    std::string out;
    for (char c : s) {
        if (c != '#' && c != '^') {
            out += c;
        }
    }
    return out;
}

std::string i2s(int v) { return std::to_string(v); }

// Patterns regexes used on description strings (a few per indel candidate, so the
// general engine is fast enough).
const util::Regex& re(std::string_view pattern) {
    thread_local std::vector<std::pair<std::string, util::Regex>> cache;  // one per worker thread
    for (const auto& [p, r] : cache) {
        if (p == pattern) {
            return r;
        }
    }
    cache.emplace_back(std::string(pattern), util::Regex::compile(pattern).value());
    return cache.back().second;
}
// Matcher.find(); group(i) or nullopt (no match, or the group did not take part).
struct Find {
    std::optional<util::Regex::Match> m;
    std::string text;
    explicit operator bool() const { return m.has_value(); }
    [[nodiscard]] std::string group(std::size_t i) const {
        auto g = m->group(text, i);
        return g ? std::string(*g) : std::string();
    }
    [[nodiscard]] int num(std::size_t i) const {
        auto v = util::parse_int(group(i));
        if (!v) {
            throw util::RecoverableError("NumberFormatException: For input string: \"" + group(i) + "\"");
        }
        return *v;
    }
};
Find find(std::string_view pattern, const std::string& text) { return Find{re(pattern).find(text), text}; }
// String.replaceFirst(pattern, "")
std::string remove_first(std::string_view pattern, const std::string& text) {
    auto m = re(pattern).find(text);
    if (!m) {
        return text;
    }
    return text.substr(0, m->begin) + text.substr(m->end);
}

CountMap single_count(int position, const std::string& description, int count) {
    CountMap m;
    m.get_or_create(position).put(description, count);
    return m;
}

int chr_length_or_npe(const RunContext& run, const std::string& chr) {
    auto it = run.chr_lengths.find(chr);
    if (it == run.chr_lengths.end()) {
        util::throw_npe("chrLengths.get(" + chr + ")");
    }
    return it->second;
}

struct Mismatch {
    std::string mismatch_sequence;  // $mm
    int mismatch_position;          // $mp
    int end;                        // 3 or 5
};

struct MismatchResult {
    std::vector<Mismatch> mismatches;
    std::vector<int> scp;
    int nm = 0;
    int misp = 0;
    std::string misnt;
};

struct SortPositionSclip {
    int position;
    Sclip* soft_clip;
    int count;
};

// VariationRealigner COMP2: softClip.varsCount desc, then position asc (a total order).
void sort_comp2(std::vector<SortPositionSclip>& v) {
    std::stable_sort(v.begin(), v.end(), [](const SortPositionSclip& a, const SortPositionSclip& b) {
        if (a.soft_clip->vars_count != b.soft_clip->vars_count) {
            return a.soft_clip->vars_count > b.soft_clip->vars_count;
        }
        return a.position < b.position;
    });
}
// VariationRealigner COMP3: count desc, then position asc.
void sort_comp3(std::vector<SortPositionSclip>& v) {
    std::stable_sort(v.begin(), v.end(), [](const SortPositionSclip& a, const SortPositionSclip& b) {
        if (a.count != b.count) {
            return a.count > b.count;
        }
        return a.position < b.position;
    });
}

class Realigner {
public:
    Realigner(VariationData& data, RegionScope& scope)
        : d_(data),
          scope_(scope),
          conf_(scope.conf),
          ref_(scope.reference.reference_sequences),
          region_(scope.region),
          chr_(chr_name(scope.conf, scope.region)),
          max_read_length_(data.max_read_length) {
        if (!scope.bam.empty()) {
            static const util::Regex colon = util::Regex::compile(":").value();
            for (auto p : util::regex_split(scope.bam, colon)) {
                bams_.emplace_back(p);
            }
        }
    }

    void realign_indels();
    void realigndel(bool use_bams, const CountMap& position_to_deletions_count);

private:
    std::string realignins(const CountMap& position_to_insertion_count);
    void realignlgdel();
    void realignlgins30();
    void realignlgins();

    MismatchResult find_mm5(int position, const std::string& wupseq);
    MismatchResult find_mm3(int p, const std::string& sanpseq);
    [[nodiscard]] bool no_passing_reads(int start, int end) const;
    [[nodiscard]] int findbp(const std::string& sequence, int start_position, int direction) const;
    [[nodiscard]] std::pair<int, std::string> findbi(const std::string& seq, int position, int dir) const;
    void adj_ref_cnt(const Variation& tv, Variation* ref, int len) const;
    void adj_ref_factor(Variation* ref, double factor_f) const;
    static void add_var_factor(Variation* vref, double factor_f);
    static void rm_cnt(Variation& vref, const Variation& tv);
    void partial_pipeline(int tts, int tte);
    std::vector<SortPositionSclip> sclips_in_extended_region(util::HashOrderMap<int, Sclip>& sclips, bool with_count);
    std::string conseq(Sclip& sc, int dir) const { return findconseq(sc, dir, scope_.run, conf_); }
    void exception(const std::exception& e, int last_position) const {
        print_exception_and_continue(e, "variant", std::to_string(last_position), &region_, scope_.exceptions);
    }

    VariationData& d_;
    RegionScope& scope_;
    const Config& conf_;
    const RefMap& ref_;
    const Region& region_;
    std::string chr_;
    int max_read_length_;
    std::vector<std::string> bams_;
};

void Realigner::realign_indels() {
    if (conf_.y) {
        std::println(stderr, "Start Realigndel");
    }
    realigndel(true, d_.position_to_deletion_count);
    if (conf_.y) {
        std::println(stderr, "Start Realignins");
    }
    realignins(d_.position_to_insertion_count);
    if (conf_.y) {
        std::println(stderr, "Start Realignlgdel");
    }
    realignlgdel();
    if (conf_.y) {
        std::println(stderr, "Start Realignlgins30");
    }
    realignlgins30();
    if (conf_.y) {
        std::println(stderr, "Start Realignlgins");
    }
    realignlgins();
}

void Realigner::realigndel(bool use_bams, const CountMap& position_to_deletions_count) {
    const RefMap& ref = ref_;
    const bool have_bams = use_bams && !bams_.empty();
    PositionMap& niv = d_.non_insertion_variants;
    auto tmp = fill_and_sort_tmp(position_to_deletions_count);
    int last_position = 0;
    for (const auto& tpl : tmp) {
        try {
            const int p = tpl.position;
            last_position = p;
            const std::string& vn = tpl.description_string;
            const int dcnt = tpl.count;
            if (conf_.y) {
                const int* cov = d_.ref_coverage.get(p);
                std::println(stderr, "  Realigndel for: {} {} {} cov: {}", p, vn, dcnt, cov ? i2s(*cov) : "null");
            }
            Variation& vref = get_variation(niv, p, vn);
            int dellen = 0;
            if (auto m = find("^-(\\d+)", vn)) {
                dellen = m.num(1);
            }
            if (auto m = find("\\^(\\d+)$", vn)) {
                dellen += m.num(1);
            }
            std::string extrains;
            std::string extra;
            std::string inv5;
            std::string inv3;
            if (auto m = find("^-\\d+\\^([ATGNC]+)<...\\d+>([ATGNC]+)$", vn)) {
                inv5 = m.group(1);
                inv3 = m.group(2);
            } else if (auto m2 = find("^-\\d+(.*)", vn)) {
                for (char c : m2.group(1)) {  // replaceAll("\\^|&|#", "")
                    if (c != '^' && c != '&' && c != '#') {
                        extra += c;
                    }
                }
                if (auto m3 = find("\\^([ATGNC]+)", vn)) {
                    extrains = m3.group(1);
                }
            }

            const int wustart = (p - 200) > 1 ? (p - 200) : 1;
            std::string wupseq = join_ref(ref, wustart, p - 1) + extra;  // 5' flanking seq
            if (!inv3.empty()) {
                wupseq = inv3;
            }
            auto len_it = scope_.run.chr_lengths.find(chr_);
            const int sanend =
                (len_it != scope_.run.chr_lengths.end() && p + 200 > len_it->second) ? len_it->second : p + 200;

            // 3' flanking seq
            const int extra_len = static_cast<int>(extra.size());
            const int extrains_len = static_cast<int>(extrains.size());
            std::string sanpseq = extra + join_ref(ref, p + dellen + extra_len - extrains_len, sanend);
            if (!inv5.empty()) {
                sanpseq = inv5;
            }

            // mismatches, mismatch positions, 5 or 3 ends
            MismatchResult r3 = find_mm3(p, sanpseq);
            MismatchResult r5 = find_mm5(p + dellen + extra_len - extrains_len - 1, wupseq);
            if (conf_.y) {
                std::println(stderr, "  Mismatches: misp3: {}-{} misp5: {}-{} sclip3: {} sclip5: {}", r3.misp, r3.misnt,
                             r5.misp, r5.misnt, r3.scp.size(), r5.scp.size());
            }

            std::vector<Mismatch> mmm = r3.mismatches;
            mmm.insert(mmm.end(), r5.mismatches.begin(), r5.mismatches.end());
            for (const Mismatch& mismatch : mmm) {
                std::string mm = mismatch.mismatch_sequence;
                const int mp = mismatch.mismatch_position;
                const int me = mismatch.end;
                if (mm.size() > 1) {
                    mm = insert_amp(mm);
                }
                VariationMap* at_mp = niv.get(mp);
                if (at_mp == nullptr) {
                    continue;
                }
                Variation* tv = at_mp->get(mm);
                if (tv == nullptr) {
                    continue;
                }
                if (tv->vars_count == 0) {
                    continue;
                }
                if (tv->mean_quality / tv->vars_count < conf_.goodq) {
                    continue;
                }
                if (tv->mean_position / tv->vars_count > (me == 3 ? r3.nm + 4 : r5.nm + 4)) {
                    continue;
                }
                if (tv->vars_count >= dcnt + dellen || tv->vars_count / dcnt >= 8) {
                    continue;
                }
                if (conf_.y) {
                    std::println(stderr, "  Realigndel Adj: {} {} {} {} {} {} {} {}", mm, mp, me, r3.nm, r5.nm, p,
                                 tv->vars_count, util::double_to_string(tv->mean_quality));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    double f = tv->mean_position != 0 ? (mp - p) / (tv->mean_position / tv->vars_count) : 1;
                    if (f > 1) {
                        f = 1;
                    }
                    inc_cnt(d_.ref_coverage, p, util::d2i_cast(tv->vars_count * f));
                    adj_ref_cnt(*tv, get_variation_maybe(niv, p, ref_get(ref, p)), dellen);
                }
                Variation* lref = nullptr;
                if (mp > p && me == 3) {
                    if (VariationMap* at_p = niv.get(p); at_p != nullptr) {
                        lref = at_p->get(std::string(1, ref_at(ref, p)));
                    }
                }
                adj_cnt(vref, *tv, lref);
                at_mp->erase(mm);
                if (at_mp->empty()) {
                    niv.erase(mp);
                }
            }
            if (r3.misp != 0 && r3.mismatches.size() == 1) {
                if (VariationMap* m = niv.get(r3.misp); m != nullptr) {
                    if (Variation* v = m->get(r3.misnt); v != nullptr && v->vars_count < dcnt) {
                        m->erase(r3.misnt);
                    }
                }
            }
            if (r5.misp != 0 && r5.mismatches.size() == 1) {
                if (VariationMap* m = niv.get(r5.misp); m != nullptr) {
                    if (Variation* v = m->get(r5.misnt); v != nullptr && v->vars_count < dcnt) {
                        m->erase(r5.misnt);
                    }
                }
            }

            for (int sc5pp : r5.scp) {
                Sclip* tv = d_.soft_clips_5end.get(sc5pp);
                if (tv == nullptr || tv->used) {
                    continue;
                }
                const std::string seq = conseq(*tv, 0);
                // Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
                if (dcnt <= 2 && tv->vars_count / dcnt > 5) {
                    continue;
                }
                if (conf_.y) {
                    std::println(stderr, "  Realigndel 5: {} {} seq: '{}' Wuseq: {} cnt: {} {} {} {}", p, sc5pp, seq,
                                 util::reverse(wupseq), tv->vars_count, dcnt, vn, p);
                }
                if (!seq.empty() && ismatch(seq, wupseq, -1, conf_)) {
                    if (sc5pp > p) {
                        inc_cnt(d_.ref_coverage, p, tv->vars_count);
                    }
                    adj_cnt(vref, *tv);
                    tv->used = true;
                }
            }

            for (int sc3pp : r3.scp) {
                Sclip* tv = d_.soft_clips_3end.get(sc3pp);
                if (tv == nullptr || tv->used) {
                    continue;
                }
                const std::string seq = conseq(*tv, 0);
                if (dcnt <= 2 && tv->vars_count / dcnt > 5) {
                    continue;
                }
                if (conf_.y) {
                    std::println(stderr, "  Realigndel 3: {} {} seq '{}' Sanseq: {} cnt: {} {} {} {} {}", p, sc3pp, seq,
                                 sanpseq, tv->vars_count, dcnt, vn, p, dellen);
                }
                if (!seq.empty() && ismatch(seq, util::substr(sanpseq, sc3pp - p), 1, conf_)) {
                    if (sc3pp <= p) {
                        inc_cnt(d_.ref_coverage, p, tv->vars_count);
                    }
                    Variation* lref = sc3pp <= p ? nullptr : get_variation_maybe(niv, p, ref_get(ref, p));
                    adj_cnt(vref, *tv, lref);
                    tv->used = true;
                }
            }
            const int pe = p + dellen + extra_len - extrains_len;
            Variation* h = get_variation_maybe(niv, p, ref_get(ref, p));
            // taking the size of gap into account
            if (have_bams && pe - p >= 5 && pe - p < max_read_length_ - 10 && h != nullptr && h->vars_count != 0 &&
                no_passing_reads(p, pe) &&
                vref.vars_count > 2 * h->vars_count * (1 - (pe - p) / static_cast<double>(max_read_length_))) {
                adj_cnt(vref, *h, h);
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }

    for (int i = static_cast<int>(tmp.size()) - 1; i > 0; i--) {
        try {
            const SortPositionDescription& tpl = tmp[static_cast<std::size_t>(i)];
            const int p = tpl.position;
            last_position = p;
            const std::string& vn = tpl.description_string;
            VariationMap* at_p = niv.get(p);
            if (at_p == nullptr) {
                continue;
            }
            Variation* vref = at_p->get(vn);
            if (vref == nullptr) {
                continue;
            }
            if (auto m = find("(-\\d+)&[ATGC]+$", vn)) {
                const std::string tn = m.group(1);
                if (Variation* tref = at_p->get(tn); tref != nullptr) {
                    if (vref->vars_count < tref->vars_count) {
                        adj_cnt(*tref, *vref);
                        at_p->erase(vn);
                    }
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
}

std::string Realigner::realignins(const CountMap& position_to_insertion_count) {
    const RefMap& ref = ref_;
    PositionMap& niv = d_.non_insertion_variants;
    PositionMap& iv = d_.insertion_variants;
    auto tmp = fill_and_sort_tmp(position_to_insertion_count);
    std::string NEWINS;
    int last_position = 0;
    for (const auto& tpl : tmp) {
        try {
            const int position = tpl.position;
            last_position = position;
            const std::string& vn = tpl.description_string;
            const int insertion_count = tpl.count;
            if (conf_.y) {
                std::println(stderr, "  Realign Ins: {} {} {}", position, vn, insertion_count);
            }
            std::string insert;
            if (auto m = find("^\\+([ATGC]+)", vn)) {
                insert = m.group(1);
            } else {
                continue;
            }
            std::string ins3;
            int inslen = static_cast<int>(insert.size());
            if (auto m = find("<dup(\\d+)>([ATGC]+)$", vn)) {
                ins3 = m.group(2);
                inslen += m.num(1) + static_cast<int>(ins3.size());
            }
            std::string extra;
            if (auto m = find("&([ATGC]+)", vn)) {
                extra = m.group(1);
            }
            std::string compm;  // the match part for a complex variant
            if (auto m = find("#([ATGC]+)", vn)) {
                compm = m.group(1);
            }
            int newdel = 0;  // the adjacent deletion
            if (auto m = find("\\^(\\d+)$", vn)) {
                newdel = m.num(1);
            }
            std::string tn = remove_first("^\\+", vn);
            tn = remove_first("&", tn);
            tn = remove_first("#", tn);
            tn = remove_first("\\^\\d+$", tn);
            tn = remove_first("\\^", tn);

            const int wustart = position - 150 > 1 ? (position - 150) : 1;
            const std::string wupseq = join_ref(ref, wustart, position) + tn;  // 5prime flanking seq
            auto len_it = scope_.run.chr_lengths.find(chr_);
            int sanend = position + static_cast<int>(vn.size()) + 100;
            if (len_it != scope_.run.chr_lengths.end() && len_it->second < sanend) {
                sanend = len_it->second;
            }
            std::string sanpseq;  // 3prime flanking seq
            MismatchResult r3;
            if (!ins3.empty()) {
                const int p3 = position + inslen - static_cast<int>(ins3.size()) + Config::kSvFlank;
                if (static_cast<int>(ins3.size()) > Config::kSvFlank) {
                    sanpseq = util::substr(ins3, Config::kSvFlank - static_cast<long>(ins3.size()));
                }
                sanpseq += join_ref(ref, position + 1, position + 101);
                r3 = find_mm3(p3 + 1, sanpseq);
            } else {
                sanpseq = tn + join_ref(ref,
                                        position + static_cast<int>(extra.size()) + 1 + static_cast<int>(compm.size()) +
                                            newdel,
                                        sanend);
                r3 = find_mm3(position + 1, sanpseq);
            }
            MismatchResult r5 =
                find_mm5(position + static_cast<int>(extra.size()) + static_cast<int>(compm.size()) + newdel, wupseq);

            std::vector<Mismatch> mmm = r3.mismatches;
            mmm.insert(mmm.end(), r5.mismatches.begin(), r5.mismatches.end());
            Variation* vref = &get_variation(iv, position, vn);
            for (const Mismatch& mismatch : mmm) {
                std::string mismatch_bases = mismatch.mismatch_sequence;
                const int mismatch_position = mismatch.mismatch_position;
                const int mismatch_end = mismatch.end;
                if (mismatch_bases.size() > 1) {
                    mismatch_bases = insert_amp(mismatch_bases);
                }
                VariationMap* at_mp = niv.get(mismatch_position);
                if (at_mp == nullptr) {
                    continue;
                }
                Variation* variation = at_mp->get(mismatch_bases);
                if (variation == nullptr) {
                    continue;
                }
                if (variation->vars_count == 0) {
                    continue;
                }
                if (variation->mean_quality / variation->vars_count < conf_.goodq) {
                    continue;
                }
                if (variation->mean_position / variation->vars_count > (mismatch_end == 3 ? r3.nm + 4 : r5.nm + 4)) {
                    continue;
                }
                if (variation->vars_count >= insertion_count + static_cast<int>(insert.size()) ||
                    variation->vars_count / insertion_count >= 8) {
                    continue;
                }
                if (conf_.y) {
                    std::println(stderr, "    insMM: {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", mismatch_bases, mismatch_position,
                                 mismatch_end, r3.nm, r5.nm, vn, insertion_count, variation->vars_count);
                }
                // Adjust ref cnt so that AF won't > 1
                if (mismatch_position > position && mismatch_end == 5) {
                    inc_cnt(d_.ref_coverage, position, variation->vars_count);
                }
                Variation* lref = nullptr;
                if (mismatch_position > position && mismatch_end == 3) {
                    if (VariationMap* at_p = niv.get(position); at_p != nullptr && ref.contains(position)) {
                        lref = at_p->get(std::string(1, ref.at(position)));
                    }
                }
                adj_cnt(*vref, *variation, lref);
                at_mp->erase(mismatch_bases);
                if (at_mp->empty()) {
                    niv.erase(mismatch_position);
                }
            }
            if (r3.misp != 0 && r3.mismatches.size() == 1) {
                if (VariationMap* m = niv.get(r3.misp); m != nullptr) {
                    if (Variation* v = m->get(r3.misnt); v != nullptr && v->vars_count < insertion_count) {
                        m->erase(r3.misnt);
                    }
                }
            }
            if (r5.misp != 0 && r5.mismatches.size() == 1) {
                if (VariationMap* m = niv.get(r5.misp); m != nullptr) {
                    if (Variation* v = m->get(r5.misnt); v != nullptr && v->vars_count < insertion_count) {
                        m->erase(r5.misnt);
                    }
                }
            }
            for (int sc5pp : r5.scp) {
                Sclip* tv = d_.soft_clips_5end.get(sc5pp);
                if (conf_.y) {
                    std::println(stderr, "    55: {} {} VN: '{}'  5' seq: ^{}^", position, sc5pp, vn, wupseq);
                }
                if (tv != nullptr && !tv->used) {
                    const std::string seq = conseq(*tv, 0);
                    if (conf_.y) {
                        std::println(stderr, "    ins5: {} {} {} {} VN: {} iCnt: {} vCnt: {}", position, sc5pp, seq,
                                     wupseq, vn, insertion_count, tv->vars_count);
                    }
                    if (!seq.empty() && ismatch(seq, wupseq, -1, conf_)) {
                        if (sc5pp > position) {
                            inc_cnt(d_.ref_coverage, position, tv->vars_count);
                        }
                        adj_cnt(*vref, *tv);
                        tv->used = true;
                    }
                }
            }
            for (int sc3pp : r3.scp) {
                Sclip* tv = d_.soft_clips_3end.get(sc3pp);
                if (conf_.y) {
                    std::println(stderr, "    33: {} {} VN: '{}'  3' seq: ^{}^", position, sc3pp, vn, sanpseq);
                }
                if (tv != nullptr && !tv->used) {
                    const std::string seq = conseq(*tv, 0);
                    if (conf_.y) {
                        std::println(stderr, "    ins3: {} {} {} {} VN: {} iCnt: {} vCnt: {}", position, sc3pp, seq,
                                     sanpseq, vn, insertion_count, tv->vars_count);
                    }
                    const std::string mseq = !ins3.empty() ? sanpseq : util::substr(sanpseq, sc3pp - position - 1);
                    if (!seq.empty() && ismatch(seq, mseq, 1, conf_)) {
                        const double insert_len = static_cast<double>(insert.size());
                        if (sc3pp <= position || insert_len > tv->mean_position / tv->vars_count) {
                            inc_cnt(d_.ref_coverage, position, tv->vars_count);
                        }
                        Variation* lref = nullptr;
                        if (sc3pp > position) {
                            if (VariationMap* at_p = niv.get(position); at_p != nullptr && ref.contains(position)) {
                                lref = at_p->get(std::string(1, ref.at(position)));
                            }
                        }
                        if (insert_len > tv->mean_position / tv->vars_count) {
                            lref = nullptr;
                        }
                        adj_cnt(*vref, *tv, lref);
                        tv->used = true;
                        const int ilen = static_cast<int>(insert.size());
                        if (ilen + 1 == static_cast<int>(vn.size()) && ilen > max_read_length_ &&
                            sc3pp >= position + 1 + ilen) {
                            int flag = 0;
                            const int offset = (sc3pp - position - 1) % ilen;
                            std::string tvn = vn;
                            for (int seqi = 0; seqi < static_cast<int>(seq.size()) && seqi + offset < ilen; seqi++) {
                                if (util::substr(seq, seqi, 1) != util::substr(insert, seqi + offset, 1)) {
                                    flag++;
                                    const int shift = seqi + offset + 1;
                                    tvn = util::substring(tvn, 0, shift) + util::substr(seq, seqi, 1) +
                                          util::substring(tvn, shift + 1);
                                }
                            }
                            if (flag > 0) {
                                // insertionVariants.get(position): put(tvn, get(vn)); remove(vn)
                                VariationMap* at_pos = iv.get(position);
                                if (at_pos == nullptr) {
                                    util::throw_npe("insertionVariants.get(position)");
                                }
                                if (!at_pos->rename(vn, tvn)) {
                                    at_pos->put(tvn, Variation{});
                                }
                                NEWINS = tvn;
                            }
                        }
                    }
                }
            }
            const int first3 = r3.scp.at(0);
            const int first5 = r5.scp.at(0);
            if (!r3.scp.empty() && !r5.scp.empty() && first3 > first5 + 3 &&
                first3 - first5 < max_read_length_ * 0.75) {
                if (ref.contains(position)) {
                    if (VariationMap* at_p = niv.get(position); at_p != nullptr) {
                        if (Variation* rv = at_p->get(std::string(1, ref.at(position)))) {
                            adj_ref_factor(rv, (first3 - first5 - 1) / static_cast<double>(max_read_length_));
                        }
                    }
                }
                adj_ref_factor(vref, -(first3 - first5 - 1) / static_cast<double>(max_read_length_));
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        } catch (const std::out_of_range& e) {
            exception(e, last_position);
        }
    }

    for (int i = static_cast<int>(tmp.size()) - 1; i > 0; i--) {
        try {
            const SortPositionDescription& tpl = tmp[static_cast<std::size_t>(i)];
            const int p = tpl.position;
            last_position = p;
            const std::string& vn = tpl.description_string;
            VariationMap* at_p = iv.get(p);
            if (at_p == nullptr) {
                continue;
            }
            Variation* vref = at_p->get(vn);
            if (vref == nullptr) {
                continue;
            }
            if (auto m = find("(\\+[ATGC]+)&[ATGC]+$", vn)) {
                const std::string tn = m.group(1);
                if (Variation* tref = at_p->get(tn); tref != nullptr) {
                    if (vref->vars_count < tref->vars_count) {
                        adj_cnt(*tref, *vref, get_variation_maybe(niv, p, ref_get(ref, p)));
                        at_p->erase(vn);
                    }
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
    return NEWINS;
}

std::vector<SortPositionSclip> Realigner::sclips_in_extended_region(util::HashOrderMap<int, Sclip>& sclips,
                                                                    bool with_count) {
    std::vector<SortPositionSclip> tmp;
    sclips.for_each([&](int p, Sclip& sc) {
        if (p < region_.start - Config::kExtension || p > region_.end + Config::kExtension) {
            return;
        }
        tmp.push_back({p, &sc, with_count ? sc.vars_count : 0});
    });
    return tmp;
}

// The reference extension and partialPipeline call of realignlgdel / realignlgins for a
// breakpoint outside the region.
void Realigner::partial_pipeline(int tts, int tte) {
    const Region modified_region = region_.with_bounds(tts, tte);
    if (!io::is_loaded(chr_, tts, tte, scope_.reference)) {
        extend_reference(scope_, modified_region, max_read_length_);
    }
    core::partial_pipeline(d_, scope_, modified_region, max_read_length_);
}

void Realigner::realignlgdel() {
    const int longmm = 3;
    const RefMap& ref = ref_;
    PositionMap& niv = d_.non_insertion_variants;
    const std::vector<std::vector<Sclip>*> dels = {&d_.sv_structures.svfdel, &d_.sv_structures.svrdel};

    auto tmp = sclips_in_extended_region(d_.soft_clips_5end, false);
    sort_comp2(tmp);
    int svcov = 0;
    int clusters = 0;
    int pairs = 0;
    int last_position = 0;
    for (const SortPositionSclip& t : tmp) {
        try {
            const int p = t.position;
            last_position = p;
            Sclip& sc5v = *t.soft_clip;
            const int cnt = sc5v.vars_count;
            if (cnt < conf_.minr) {
                break;
            }
            if (sc5v.used) {  // already been used
                continue;
            }
            const std::string seq = conseq(sc5v, 5);
            if (seq.empty()) {
                continue;
            }
            if (seq.size() < 7) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Working Realignlgdel: 5' {} '{}' {}", p, seq, cnt);
            }
            int bp = findbp(seq, p - 5, -1);
            std::string extra;
            std::string EXTRA;
            if (bp == 0) {
                if (islowcomplexseq(seq)) {
                    continue;
                }
                Match match = find_match(seq, scope_.reference, p, -1, Config::kSeed1, 1, conf_);
                bp = match.base_position;
                EXTRA = match.matched_sequence;
                if (!(bp != 0 && p - bp > 15 && p - bp < Config::kSvMaxLen)) {
                    continue;
                }
                bp++;
                const SvMark sv_mark = mark_sv(bp, p, dels, max_read_length_, conf_);
                svcov = sv_mark.cov;
                clusters = sv_mark.count;
                pairs = sv_mark.pairs;
                if (svcov == 0) {
                    if (cnt <= conf_.minr) {
                        continue;
                    }
                }
                SvInfo& sv = get_sv(niv, bp);
                sv.type = "DEL";
                sv.pairs += pairs;
                sv.splits += cnt;
                sv.clusters += clusters;

                if (bp < region_.start) {
                    const int tts = bp - max_read_length_;
                    int tte = bp + max_read_length_;
                    if (bp + max_read_length_ >= region_.start) {
                        tte = region_.start - 1;
                    }
                    partial_pipeline(tts, tte);
                }
            }
            int dellen = p - bp;
            int en = 0;
            std::string gt = "-" + i2s(dellen);
            if (EXTRA.empty()) {
                while (en < static_cast<int>(seq.size()) && !equals_char(ref_get(ref, bp - en - 1), seq[en])) {
                    extra += util::substr(seq, en, 1);
                    en++;
                }
                if (!extra.empty()) {
                    extra = util::reverse(extra);
                    gt = "-" + i2s(dellen) + "&" + extra;
                    bp -= static_cast<int>(extra.size());
                }
            } else {
                dellen -= static_cast<int>(EXTRA.size());
                gt = dellen == 0 ? "-" + i2s(static_cast<int>(EXTRA.size())) + "^" + EXTRA
                                 : "-" + i2s(dellen) + "&" + EXTRA;
            }
            if (conf_.y) {
                std::println(stderr, "  Found Realignlgdel: {} {} 5' {} {} {}", bp, gt, p, seq, cnt);
            }

            // Work on the softclipped read at 3'
            int n = 0;
            if (extra.empty() && EXTRA.empty()) {
                while (ref.contains(bp + n) && ref.contains(bp + dellen + n) &&
                       ref.at(bp + n) == ref.at(bp + dellen + n)) {
                    n++;
                }
            }
            int sc3p = bp + n;
            std::string str;
            int mcnt = 0;
            while (mcnt <= longmm && ref.contains(bp + n) && ref.contains(bp + dellen + n) &&
                   ref.at(bp + n) != ref.at(bp + dellen + n)) {
                str += ref.at(bp + dellen + n);
                n++;
                mcnt++;
            }
            if (str.size() == 1) {
                int nm = 0;
                while (ref.contains(bp + n) && ref.contains(bp + dellen + n) &&
                       ref.at(bp + n) == ref.at(bp + dellen + n)) {
                    n++;
                    if (n != 0) {
                        nm++;
                    }
                }
                if (nm >= 3 && !d_.soft_clips_3end.contains(sc3p)) {
                    sc3p = bp + n;
                }
            }

            // likely a false positive
            if (VariationMap* at_bp = niv.get(bp);
                at_bp != nullptr && at_bp->sv && !d_.soft_clips_3end.contains(sc3p)) {
                if (svcov == 0 && cnt <= conf_.minr) {
                    remove_sv(niv, bp);
                    continue;
                }
            }
            Variation& tv = get_variation(niv, bp, gt);
            tv.qstd = true;  // more accurate implementation later
            tv.pstd = true;  // more accurate implementation later
            adj_cnt(tv, sc5v);
            sc5v.used = bp != 0;
            if (!d_.ref_coverage.contains(bp) && d_.ref_coverage.contains(p)) {
                d_.ref_coverage.put(bp, *d_.ref_coverage.get(p));
            }
            if (dellen < conf_.indelsize) {
                for (int tp = bp; tp < bp + dellen; tp++) {
                    inc_cnt(d_.ref_coverage, tp, sc5v.vars_count);
                }
            }

            if (Sclip* sclip = d_.soft_clips_3end.get(sc3p); sclip != nullptr && !sclip->used) {
                if (sc3p > bp) {
                    adj_cnt(tv, *sclip, get_variation_maybe(niv, bp, ref_get(ref, bp)));
                } else {
                    adj_cnt(tv, *sclip);
                }
                if (sc3p == bp) {
                    if (dellen < conf_.indelsize) {
                        for (int tp = bp; tp < bp + dellen; tp++) {
                            inc_cnt(d_.ref_coverage, tp, sclip->vars_count);
                        }
                    }
                }
                for (int ip = bp + 1; ip < sc3p; ip++) {
                    const std::string base(1, ref_at(ref, dellen + ip));
                    Variation& vv = get_variation(niv, ip, base);
                    rm_cnt(vv, *sclip);
                    if (vv.vars_count == 0) {
                        niv.get(ip)->erase(base);
                    }
                    if (niv.get(ip)->empty()) {
                        niv.erase(ip);
                    }
                }
                sclip->used = bp != 0;
            }
            const int tv_count = tv.vars_count;
            const CountMap dels5 = single_count(bp, gt, tv_count);
            realigndel(true, dels5);

            if (VariationMap* at_bp = niv.get(bp); at_bp != nullptr && at_bp->sv) {
                at_bp->sv->splits += tv.vars_count - tv_count;
            }
            if (svcov > tv.vars_count) {
                add_var_factor(&tv, (svcov - tv.vars_count) / static_cast<double>(tv.vars_count));
            }
            if (conf_.y) {
                std::println(stderr, "  Found lgdel done: {} {} {} 5' {} {}\n", bp, gt, p, seq, tv.vars_count);
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }

    // Work on 3' clipped reads
    tmp = sclips_in_extended_region(d_.soft_clips_3end, false);
    sort_comp2(tmp);
    svcov = 0;
    clusters = 0;
    pairs = 0;
    for (const SortPositionSclip& t : tmp) {
        try {
            const int p = t.position;
            last_position = p;
            Sclip& sc3v = *t.soft_clip;
            const int cnt = sc3v.vars_count;
            if (cnt < conf_.minr) {
                break;
            }
            if (sc3v.used) {
                continue;
            }
            const std::string seq = conseq(sc3v, 3);
            if (seq.empty()) {
                continue;
            }
            if (seq.size() < 7) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Working Realignlgdel: 3' {} '{}' {}", p, seq, cnt);
            }
            int bp = findbp(seq, p + 5, 1);
            std::string extra;
            std::string EXTRA;
            if (bp == 0) {
                if (islowcomplexseq(seq)) {
                    continue;
                }
                Match match = find_match(seq, scope_.reference, p, 1, Config::kSeed1, 1, conf_);
                bp = match.base_position;
                EXTRA = match.matched_sequence;
                if (!(bp != 0 && bp - p > 15 && p - bp < Config::kSvMaxLen)) {
                    continue;
                }
                const SvMark sv_mark = mark_sv(p, bp, dels, max_read_length_, conf_);
                svcov = sv_mark.cov;
                clusters = sv_mark.count;
                pairs = sv_mark.pairs;
                if (svcov == 0) {
                    if (cnt <= conf_.minr) {  // a little more stringent
                        continue;
                    }
                }
                SvInfo& sv = get_sv(niv, p);
                sv.type = "DEL";
                sv.pairs += pairs;
                sv.splits += cnt;
                sv.clusters += clusters;

                if (bp > region_.end) {
                    int tts = bp - max_read_length_;
                    const int tte = bp + max_read_length_;
                    if (bp - max_read_length_ <= region_.end) {
                        tts = region_.end + 1;
                    }
                    partial_pipeline(tts, tte);
                }
            }

            int dellen = bp - p;
            int en = 0;
            if (!EXTRA.empty()) {
                dellen -= static_cast<int>(EXTRA.size());
            } else {
                while (en < static_cast<int>(seq.size()) && !equals_char(ref_get(ref, bp + en), seq[en])) {
                    extra += seq[en];
                    en++;
                }
            }
            std::string gt = "-" + i2s(dellen);
            int sc5p = bp;
            bp = p;  // Set it to 5'
            if (!extra.empty()) {
                gt = "-" + i2s(dellen) + "&" + extra;
                sc5p += static_cast<int>(extra.size());
            } else if (!EXTRA.empty()) {
                gt = "-" + i2s(dellen) + "&" + EXTRA;
            } else {
                // 5' adjustment
                while (ref.contains(bp - 1) && ref.contains(bp + dellen - 1) &&
                       ref.at(bp - 1) == ref.at(bp + dellen - 1)) {
                    bp--;
                    if (bp != 0) {
                        sc5p--;
                    }
                }
                if (bp != p) {
                    if (VariationMap* at_p = niv.get(p); at_p != nullptr && at_p->sv) {
                        const SvInfo old = *at_p->sv;
                        SvInfo& sv = get_sv(niv, bp);
                        sv.clusters = old.clusters;
                        sv.pairs = old.pairs;
                        sv.splits = old.splits;
                        sv.type = old.type;
                        remove_sv(niv, p);
                    }
                }
            }

            if (VariationMap* at_bp = niv.get(bp);
                at_bp != nullptr && at_bp->sv && !d_.soft_clips_5end.contains(sc5p)) {
                if (svcov == 0 && cnt <= conf_.minr) {
                    remove_sv(niv, bp);
                    continue;
                }
            }
            if (conf_.y) {
                std::println(stderr, "  Found Realignlgdel: bp: {} {} 3' {} 5'clip: {} '{}' {}", bp, gt, p, sc5p, seq,
                             cnt);
            }

            Variation& tv = get_variation(niv, bp, gt);
            tv.qstd = true;  // more accurate implementation later
            tv.pstd = true;  // more accurate implementation later
            if (dellen < conf_.indelsize) {
                for (int tp = bp; tp < bp + dellen + static_cast<int>(extra.size()) + static_cast<int>(EXTRA.size());
                     tp++) {
                    inc_cnt(d_.ref_coverage, tp, sc3v.vars_count);
                }
            }
            if (!d_.ref_coverage.contains(bp)) {
                if (const int* c = d_.ref_coverage.get(p - 1)) {
                    d_.ref_coverage.put(bp, *c);
                } else {
                    d_.ref_coverage.put(bp, sc3v.vars_count);
                }
            }
            sc3v.mean_position += dellen * sc3v.vars_count;
            adj_cnt(tv, sc3v);
            sc3v.used = p + dellen != 0;

            const int tv_count = tv.vars_count;
            const CountMap dels5 = single_count(bp, gt, tv_count);
            realigndel(true, dels5);

            if (VariationMap* at_bp = niv.get(bp); at_bp != nullptr && at_bp->sv) {
                at_bp->sv->splits += tv.vars_count - tv_count;
            }
            if (conf_.y) {
                std::println(stderr, "  Found lgdel: {} {} {} 3' '{}' {}\n", bp, gt, p, seq, tv.vars_count);
            }
            if (svcov > tv.vars_count) {
                add_var_factor(&tv, (svcov - tv.vars_count) / static_cast<double>(tv.vars_count));
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
    if (conf_.y) {
        std::println(stderr, "  Done: Realignlgdel\n");
    }
}

void Realigner::realignlgins30() {
    const RefMap& ref = ref_;
    PositionMap& niv = d_.non_insertion_variants;
    PositionMap& iv = d_.insertion_variants;

    auto tmp5 = sclips_in_extended_region(d_.soft_clips_5end, true);
    sort_comp3(tmp5);
    auto tmp3 = sclips_in_extended_region(d_.soft_clips_3end, true);
    sort_comp3(tmp3);
    int last_position = 0;
    for (const SortPositionSclip& t5 : tmp5) {
        try {
            const int p5 = t5.position;
            last_position = p5;
            Sclip& sc5v = *t5.soft_clip;
            const int cnt5 = t5.count;
            if (sc5v.used) {
                continue;
            }
            for (const SortPositionSclip& t3 : tmp3) {
                try {
                    const int p3 = t3.position;
                    last_position = p3;
                    Sclip& sc3v = *t3.soft_clip;
                    const int cnt3 = t3.count;
                    if (sc5v.used) {
                        break;
                    }
                    if (sc3v.used) {
                        continue;
                    }
                    if (p5 - p3 > max_read_length_ * 2.5) {
                        continue;
                    }
                    if (p3 - p5 > max_read_length_ - 10) {  // if they're too far away, don't even try
                        continue;
                    }
                    const std::string seq5 = conseq(sc5v, 5);
                    const std::string seq3 = conseq(sc3v, 3);
                    // next until at least one of consensus sequences has length > 10
                    if (seq5.size() <= 10 || seq3.size() <= 10) {
                        continue;
                    }
                    if (conf_.y) {
                        std::println(stderr, "  Working lgins30: {} {} 3: {} {} 5: {} {}", p3, p5, seq3, cnt3,
                                     util::reverse(seq5), cnt5);
                    }
                    // Don't even try if there're extreme bias
                    if (!(cnt5 / static_cast<double>(cnt3) >= 0.08 && cnt5 / static_cast<double>(cnt3) <= 12)) {
                        continue;
                    }
                    const Match35 match35 = find35match(seq5, seq3, conf_);
                    const int bp5 = match35.matched5end;
                    const int bp3 = match35.matched3end;
                    const int score = match35.max_matched_length;  // length of match
                    if (score == 0) {
                        continue;
                    }
                    // to ensure higher quality of bases are used as read ends are usually poor quality
                    const int smscore = score / 2;
                    // add matched part of seq3 (left clip)
                    std::string ins = bp3 + smscore > 1 ? util::substr(seq3, 0, -(bp3 + smscore) + 1) : seq3;
                    if (bp5 + smscore > 0) {  // add not matched part of seq5 (right clip)
                        ins += util::reverse(util::substr(seq5, 0, bp5 + smscore));
                    }
                    if (islowcomplexseq(ins)) {
                        if (conf_.y) {
                            std::println(stderr, "  Discard low complexity insertion found {}.", ins);
                        }
                        continue;
                    }
                    int bi = 0;
                    Variation* vref = nullptr;
                    if (conf_.y) {
                        std::println(stderr, "  Found candidate lgins30: {} {} {}", p3, p5, ins);
                    }
                    const int seq3_len = static_cast<int>(seq3.size());
                    const int seq5_len = static_cast<int>(seq5.size());
                    if (p5 > p3) {
                        int ins_len = static_cast<int>(ins.size());
                        if (seq3_len > ins_len && !ismatch(util::substr(seq3, ins_len),
                                                           join_ref(ref, p5, p5 + seq3_len - ins_len + 2), 1, conf_)) {
                            continue;
                        }
                        if (seq5_len > ins_len &&
                            !ismatch(util::substr(seq5, ins_len), join_ref(ref, p3 - (seq5_len - ins_len) - 2, p3 - 1),
                                     -1, conf_)) {
                            continue;
                        }
                        if (conf_.y) {
                            std::println(stderr, "  Found lgins30 complex: {} {} {} {}", p3, p5, ins_len, ins);
                        }
                        const std::string tmp = join_ref(ref, p3, p5 - 1);
                        if (tmp.size() > ins.size()) {  // deletion is longer
                            ins = i2s(p3 - p5) + "^" + ins;
                            bi = p3;
                            vref = &get_variation(niv, p3, ins);
                        } else if (tmp.size() < ins.size()) {
                            ins = util::substr(ins, 0, static_cast<long>(ins.size() - tmp.size())) + "&" +
                                  util::substr(ins, p3 - p5);
                            ins = "+" + ins;
                            bi = p3 - 1;
                            vref = &get_variation(iv, bi, ins);
                        } else {  // long MNP
                            ins = "-" + i2s(static_cast<int>(ins.size())) + "^" + ins;
                            bi = p3;
                            vref = &get_variation(niv, p3, ins);
                        }
                    } else {
                        const int ins_len = static_cast<int>(ins.size());
                        if (seq3_len > ins_len && !ismatch(util::substr(seq3, ins_len),
                                                           join_ref(ref, p5, p5 + seq3_len - ins_len + 2), 1, conf_)) {
                            continue;
                        }
                        if (seq5_len > ins_len &&
                            !ismatch(util::substr(seq5, ins_len), join_ref(ref, p3 - (seq5_len - ins_len) - 2, p3 - 1),
                                     -1, conf_)) {
                            continue;
                        }
                        std::string tmp;
                        if (ins_len <= p3 - p5) {  // Tandem duplication
                            int rpt = 2;
                            int tnr = 3;
                            while (((p3 - p5 + ins_len) / static_cast<double>(tnr)) / static_cast<double>(ins_len) >
                                   1) {
                                if ((p3 - p5 + ins_len) % tnr == 0) {
                                    rpt++;
                                }
                                tnr++;
                            }
                            tmp += join_ref_double(ref, p5,
                                                   (p5 + (p3 - p5 + ins_len) / static_cast<double>(rpt) - ins_len));
                            ins = "+" + tmp + ins;
                        } else {
                            tmp += join_ref(ref, p5, p3 - 1);
                            if ((ins_len - static_cast<int>(tmp.size())) % 2 == 0) {
                                const int tex = (ins_len - static_cast<int>(tmp.size())) / 2;
                                ins = (tmp + util::substr(ins, 0, tex)) == util::substr(ins, tex)
                                          ? ("+" + util::substr(ins, tex))
                                          : "+" + tmp + ins;
                            } else {
                                ins = "+" + tmp + ins;
                            }
                        }
                        if (conf_.y) {
                            std::println(stderr, "Found lgins30: {} {} {} {} + {}", p3, p5, ins.size(), tmp, ins);
                        }
                        bi = p5 - 1;
                        vref = &get_variation(iv, bi, ins);
                    }
                    sc3v.used = true;
                    sc5v.used = true;
                    vref->pstd = true;
                    vref->qstd = true;
                    inc_cnt(d_.ref_coverage, bi, sc5v.vars_count);
                    if (conf_.y) {
                        std::println(stderr, " lgins30 Found: '{}' {} {} {}", ins, bi, bp3, bp5);
                    }

                    if (ins.starts_with("+")) {
                        Variation* mvref = get_variation_maybe(niv, bi, ref_get(ref, bi));
                        adj_cnt(*vref, sc3v, mvref);
                        adj_cnt(*vref, sc5v);
                        if (!bams_.empty() && p3 - p5 >= 5 && p3 - p5 < max_read_length_ - 10 && mvref != nullptr &&
                            mvref->vars_count != 0 && no_passing_reads(p5, p3) &&
                            vref->vars_count > 2 * mvref->vars_count) {
                            adj_cnt(*vref, *mvref, mvref);
                        }
                        realignins(single_count(bi, ins, vref->vars_count));
                    } else if (ins.starts_with("-")) {
                        adj_cnt(*vref, sc3v, get_variation_maybe(niv, bi, ref_get(ref, bi)));
                        adj_cnt(*vref, sc5v);
                        realigndel(false, single_count(bi, ins, vref->vars_count));
                    } else {
                        adj_cnt(*vref, sc3v);
                        adj_cnt(*vref, sc5v);
                    }
                    break;
                } catch (const util::RecoverableError& e) {
                    exception(e, last_position);
                }
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
    if (conf_.y) {
        std::println(stderr, "Done: lgins30\n");
    }
}

void Realigner::realignlgins() {
    const RefMap& ref = ref_;
    PositionMap& niv = d_.non_insertion_variants;
    PositionMap& iv = d_.insertion_variants;
    const std::vector<std::vector<Sclip>*> dups = {&d_.sv_structures.svfdup, &d_.sv_structures.svrdup};

    auto tmp = sclips_in_extended_region(d_.soft_clips_5end, false);
    sort_comp2(tmp);
    int last_position = 0;
    for (const SortPositionSclip& t : tmp) {
        try {
            const int p = t.position;
            last_position = p;
            Sclip& sc5v = *t.soft_clip;
            const int cnt = sc5v.vars_count;
            if (cnt < conf_.minr) {
                break;
            }
            if (sc5v.used) {  // already been used
                continue;
            }
            const std::string seq = conseq(sc5v, 0);
            if (seq.empty()) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Working lgins: 5: {} {} cnt: {}", p, seq, cnt);
            }
            if (seq.size() < 12) {
                continue;
            }
            auto [bi, ins] = findbi(seq, p, -1);
            std::string EXTRA;
            if (bi == 0) {
                if (islowcomplexseq(seq)) {
                    continue;
                }
                Match match = find_match(seq, scope_.reference, p, -1, Config::kSeed1, 1, conf_);
                bi = match.base_position;
                EXTRA = match.matched_sequence;
                if (!(bi != 0 && bi - p > 15 && bi - p < Config::kSvMaxLen)) {
                    continue;
                }
                // For large insertions
                if (bi > region_.end) {
                    int tts = bi - max_read_length_;
                    const int tte = bi + max_read_length_;
                    if (bi - max_read_length_ <= region_.end) {
                        tts = region_.end + 1;
                    }
                    partial_pipeline(tts, tte);
                }
                if (bi - p > conf_.sv_min_len + 2 * Config::kSvFlank) {
                    ins = join_ref(ref, p, p + Config::kSvFlank - 1);
                    ins += "<dup" + i2s(bi - p - 2 * Config::kSvFlank + 1) + ">";
                    ins += join_ref_for_5_lgins(ref, bi - Config::kSvFlank + 1, bi, seq, EXTRA);
                } else {
                    ins = join_ref_for_5_lgins(ref, p, bi, seq, EXTRA);
                }
                ins += EXTRA;

                const SvMark tp2 = mark_dup_sv(p, bi, dups, max_read_length_, conf_);
                const int clusters = tp2.count;
                const int pairs = tp2.pairs;
                const int* cov_p1 = d_.ref_coverage.get(p - 1);
                const int* cov_bi = d_.ref_coverage.get(bi);
                if (cov_p1 == nullptr || (cov_bi != nullptr && *cov_p1 < *cov_bi)) {
                    if (cov_bi != nullptr) {
                        d_.ref_coverage.put(p - 1, *cov_bi);
                    } else {
                        d_.ref_coverage.put(p - 1, sc5v.vars_count);
                    }
                } else {
                    if (sc5v.vars_count > *cov_p1) {
                        inc_cnt(d_.ref_coverage, p - 1, sc5v.vars_count);
                    }
                }

                bi = p - 1;
                SvInfo& sv = get_sv(niv, bi);
                sv.type = "DUP";
                sv.pairs += pairs;
                sv.splits += cnt;
                sv.clusters += clusters;
            }

            if (conf_.y) {
                std::println(stderr, "  Found candidate lgins from 5: {} +{} {} {}", bi, ins, p, seq);
            }
            Variation& iref = get_variation(iv, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            adj_cnt(iref, sc5v);
            bool rpflag = true;  // A flag to indicate whether an insertion is a repeat
            for (int i = 0; i < static_cast<int>(ins.size()); i++) {
                if (!equals_char(ref_get(ref, bi + 1 + i), ins[static_cast<std::size_t>(i)])) {
                    rpflag = false;
                    break;
                }
            }
            if (VariationMap* at_bi = niv.get(bi); at_bi != nullptr && !at_bi->sv) {
                inc_cnt(d_.ref_coverage, bi, sc5v.vars_count);
            }
            int len = static_cast<int>(ins.size());
            if (ins.find('&') != std::string::npos) {
                len--;
            }
            if (sc5v.seq.empty()) {
                throw util::RecoverableError("NoSuchElementException: TreeMap.lastKey");
            }
            const int seq_len = sc5v.seq.rbegin()->first + 1;
            for (int ii = len + 1; ii < seq_len; ii++) {
                const int pii = bi - ii + len;
                auto sit = sc5v.seq.find(ii);
                if (sit == sc5v.seq.end()) {
                    continue;
                }
                sit->second.for_each([&](char tnt, Variation& tv) {
                    Variation& tvr = get_variation(niv, pii, std::string(1, tnt));
                    adj_cnt(tvr, tv);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    inc_cnt(d_.ref_coverage, pii, tv.vars_count);
                });
            }
            sc5v.used = bi + len != 0;

            const int iref_count = iref.vars_count;
            std::string newins = realignins(single_count(bi, "+" + ins, iref_count));
            if (newins.empty()) {
                newins = "+" + ins;
            }
            Variation* mref = get_variation_maybe(niv, bi, ref_get(ref, bi));
            Variation& kref = get_variation(iv, bi, newins);

            if (VariationMap* at_bi = niv.get(bi); at_bi != nullptr && at_bi->sv) {
                at_bi->sv->splits += kref.vars_count - iref_count;
            }
            const int ins_len = static_cast<int>(ins.size());
            if (rpflag && !bams_.empty() && ins_len >= 5 && ins_len < max_read_length_ - 10 && mref != nullptr &&
                mref->vars_count != 0 && no_passing_reads(bi, bi + ins_len) && kref.vars_count > 2 * mref->vars_count) {
                adj_cnt(kref, *mref, mref);
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }

    tmp = sclips_in_extended_region(d_.soft_clips_3end, false);
    sort_comp2(tmp);
    for (const SortPositionSclip& t : tmp) {
        try {
            int p = t.position;
            last_position = p;
            Sclip& sc3v = *t.soft_clip;
            const int cnt = sc3v.vars_count;
            if (cnt < conf_.minr) {
                break;
            }
            if (sc3v.used) {
                continue;
            }
            const std::string seq = conseq(sc3v, 0);
            if (seq.empty()) {
                continue;
            }
            if (conf_.y) {
                std::println(stderr, "  Working lgins 3: {} {} cnt: {}", p, seq, cnt);
            }
            if (seq.size() < 12) {
                continue;
            }
            auto [bi, ins] = findbi(seq, p, 1);
            std::string EXTRA;
            if (bi == 0) {
                if (islowcomplexseq(seq)) {
                    continue;
                }
                Match match = find_match(seq, scope_.reference, p, 1, Config::kSeed1, 1, conf_);
                bi = match.base_position;
                EXTRA = match.matched_sequence;
                if (!(bi != 0 && p - bi > 15 && p - bi < Config::kSvMaxLen)) {
                    continue;
                }
                // For large insertions
                if (bi < region_.start) {
                    const int tts = bi - max_read_length_;
                    int tte = bi + max_read_length_;
                    if (bi + max_read_length_ >= region_.start) {
                        tte = region_.start - 1;
                    }
                    partial_pipeline(tts, tte);
                }
                int shift5 = 0;
                while (ref.contains(p - 1) && ref.contains(bi - 1) && ref.at(p - 1) == ref.at(bi - 1)) {
                    p--;
                    bi--;
                    shift5++;
                }
                if (p - bi > conf_.sv_min_len + 2 * Config::kSvFlank) {
                    ins = join_ref_for_3_lgins(ref, bi, bi + Config::kSvFlank - 1, shift5, seq, EXTRA);
                    ins += "<dup" + i2s(p - bi - 2 * Config::kSvFlank) + ">";
                    ins += join_ref(ref, p - Config::kSvFlank, p - 1);
                } else {
                    ins = join_ref_for_3_lgins(ref, bi, p - 1, shift5, seq, EXTRA);
                }
                ins += EXTRA;
                const SvMark tp2 = mark_dup_sv(bi, p - 1, dups, max_read_length_, conf_);
                const int clusters = tp2.count;
                const int pairs = tp2.pairs;
                bi = bi - 1;

                SvInfo& sv = get_sv(niv, bi);
                sv.type = "DUP";
                sv.pairs += pairs;
                sv.splits += cnt;
                sv.clusters += clusters;
                const int* cov_bi = d_.ref_coverage.get(bi);
                const int* cov_p = d_.ref_coverage.get(p);
                if (cov_bi == nullptr || (cov_p != nullptr && *cov_bi < *cov_p)) {
                    if (cov_p != nullptr) {
                        d_.ref_coverage.put(bi, *cov_p);
                    } else {
                        d_.ref_coverage.put(bi, sc3v.vars_count);
                    }
                } else {
                    if (sc3v.vars_count > *cov_bi) {
                        inc_cnt(d_.ref_coverage, bi, sc3v.vars_count);
                    }
                }
            }

            if (conf_.y) {
                std::println(stderr, "  Found candidate lgins from 3: {} +{} {} {}", bi, ins, p, seq);
            }
            Variation& iref = get_variation(iv, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            Variation* lref = get_variation_maybe(niv, bi, ref_get(ref, bi));
            if (p - bi > sc3v.mean_position / cnt) {
                lref = nullptr;
            }
            adj_cnt(iref, sc3v, lref);
            bool rpflag = true;
            for (int i = 0; i < static_cast<int>(ins.size()); i++) {
                if (!equals_char(ref_get(ref, bi + 1 + i), ins[static_cast<std::size_t>(i)])) {
                    rpflag = false;
                    break;
                }
            }
            int len = static_cast<int>(ins.size());
            if (ins.find('&') != std::string::npos) {
                len--;
            }
            if (sc3v.seq.empty()) {
                throw util::RecoverableError("NoSuchElementException: TreeMap.lastKey");
            }
            const int len_seq = sc3v.seq.rbegin()->first + 1;
            for (int ii = len; ii < len_seq; ii++) {
                const int pii = p + ii - len;
                auto sit = sc3v.seq.find(ii);
                if (sit == sc3v.seq.end()) {
                    continue;
                }
                sit->second.for_each([&](char tnt, Variation& tv) {
                    Variation& vref = get_variation(niv, pii, std::string(1, tnt));
                    adj_cnt(vref, tv);
                    vref.pstd = true;
                    vref.qstd = true;
                    inc_cnt(d_.ref_coverage, pii, tv.vars_count);
                });
            }
            sc3v.used = true;
            const int iref_count = iref.vars_count;
            realignins(single_count(bi, "+" + ins, iref_count));
            if (VariationMap* at_bi = niv.get(bi); at_bi != nullptr && at_bi->sv) {
                at_bi->sv->splits += iref.vars_count - iref_count;
            }
            Variation* mref = get_variation_maybe(niv, bi, ref_get(ref, bi));
            const int ins_len = static_cast<int>(ins.size());
            if (rpflag && !bams_.empty() && ins_len >= 5 && ins_len < max_read_length_ - 10 && mref != nullptr &&
                mref->vars_count != 0 && no_passing_reads(bi, bi + ins_len) && iref.vars_count > 2 * mref->vars_count) {
                adj_cnt(iref, *mref, mref);
            }
        } catch (const util::RecoverableError& e) {
            exception(e, last_position);
        }
    }
}

MismatchResult Realigner::find_mm5(int position, const std::string& wupseq) {
    const RefMap& ref = ref_;
    const std::string seq = strip_hash_caret(wupseq);
    const int longmm = 3;
    MismatchResult r;
    int n = 0;
    int mn = 0;
    int mcnt = 0;
    std::string str;
    while (has_and_not_equals(util::utils_char_at(seq, -1 - n), ref, position - n) && mcnt < longmm) {
        str.insert(str.begin(), as_char(util::utils_char_at(seq, -1 - n)));
        r.mismatches.push_back({str, position - n, 5});
        n++;
        mcnt++;
    }
    r.scp.push_back(position + 1);
    // Adjust clipping position if only one mismatch
    if (str.size() == 1) {
        while (has_and_equals(util::utils_char_at(seq, -1 - n), ref, position - n)) {
            n++;
            if (n != 0) {
                mn++;
            }
        }
        if (mn > 1) {
            int n2 = 0;
            while (-1 - n - 1 - n2 >= 0 &&
                   has_and_equals(util::utils_char_at(seq, -1 - n - 1 - n2), ref, position - n - 1 - n2)) {
                n2++;
            }
            if (n2 > 2) {
                r.scp.push_back(position - n - n2);
                r.misp = position - n;
                r.misnt = std::string(1, as_char(util::utils_char_at(seq, -1 - n)));
                if (Sclip* sc = d_.soft_clips_5end.get(position - n - n2)) {
                    sc->used = true;
                }
                mn += n2;
            } else {
                r.scp.push_back(position - n);
                if (Sclip* sc = d_.soft_clips_5end.get(position - n)) {
                    sc->used = true;
                }
            }
        }
    }
    r.nm = mn;
    return r;
}

MismatchResult Realigner::find_mm3(int p, const std::string& sanpseq) {
    const RefMap& ref = ref_;
    const std::string seq = strip_hash_caret(sanpseq);
    const int longmm = 3;
    const int len = static_cast<int>(seq.size());
    MismatchResult r;
    int n = 0;
    int mn = 0;
    int mcnt = 0;
    std::string str;
    while (n < len && ref.contains(p + n) && ref.at(p + n) == seq[static_cast<std::size_t>(n)]) {
        n++;
    }
    r.scp.push_back(p + n);
    const int Tbp = p + n;
    while (mcnt <= longmm && n < len && !equals_char(ref_get(ref, p + n), seq[static_cast<std::size_t>(n)])) {
        str += seq[static_cast<std::size_t>(n)];
        r.mismatches.push_back({str, Tbp, 3});
        n++;
        mcnt++;
    }
    // Adjust clipping position if only one mismatch
    if (str.size() == 1) {
        while (n < len && is_has_and_equals(seq[static_cast<std::size_t>(n)], ref, p + n)) {
            n++;
            if (n != 0) {
                mn++;
            }
        }
        if (mn > 1) {
            int n2 = 0;
            while (n + n2 + 1 < len &&
                   is_has_and_equals(seq[static_cast<std::size_t>(n + n2 + 1)], ref, p + n + 1 + n2)) {
                n2++;
            }
            if (n2 > 2 && n + n2 + 1 < len) {
                r.scp.push_back(p + n + n2);
                r.misp = p + n;
                r.misnt = std::string(1, util::char_at(seq, n));
                if (Sclip* sc = d_.soft_clips_3end.get(p + n + n2)) {
                    sc->used = true;
                }
                mn += n2;
            } else {
                r.scp.push_back(p + n);
                if (Sclip* sc = d_.soft_clips_3end.get(p + n)) {
                    sc->used = true;
                }
            }
        }
    }
    r.nm = mn;
    return r;
}

bool Realigner::no_passing_reads(int start, int end) const {
    int cnt = 0;
    int midcnt = 0;  // Reads end in the middle
    const int dlen = end - start;
    const std::string dlenqr = i2s(dlen) + "D";
    const Region region{chr_, start, end, ""};
    if (scope_.bam_cache == nullptr) {
        util::throw_npe("SamView");
    }
    for (const std::string& bam : bams_) {
        auto file = scope_.bam_cache->get(bam);
        if (!file) {
            throw util::RecoverableError(file.error().message);
        }
        auto view = io::SamView::query(**file, "0", region);
        if (!view) {
            throw util::RecoverableError(view.error().message);
        }
        io::BamRecord record;
        while (view->read(record)) {
            if (record.cigar_string().find(dlenqr) != std::string::npos) {
                continue;
            }
            const int read_start = record.alignment_start();
            // The total aligned length, excluding soft-clipped bases and insertions
            const int read_end = read_start + aligned_length(cigar_of(record));
            if (read_end > end + 2 && read_start < start - 2) {
                cnt++;
            }
            if (read_start < start - 2 && read_end > start && read_end < end) {
                midcnt++;
            }
        }
        if (view->error()) {
            throw util::RecoverableError(*view->error());
        }
    }
    if (conf_.y) {
        std::println(stderr, "    Passing Read CNT: {} {} {} {} {}", cnt, chr_, start, end, midcnt);
    }
    return cnt <= 0 && midcnt + 1 > 0;
}

int Realigner::findbp(const std::string& sequence, int start_position, int direction) const {
    const RefMap& ref = ref_;
    const int maxmm = 3;  // maximum mismatches allowed
    int bp = 0;
    int score = 0;
    auto len_it = scope_.run.chr_lengths.find(chr_);
    const int idx = len_it != scope_.run.chr_lengths.end() ? len_it->second : 0;
    const int seq_len = static_cast<int>(sequence.size());
    for (int n = 0; n < conf_.indelsize; n++) {
        int mm = 0;
        int i = 0;
        std::set<char> m;
        for (i = 0; i < seq_len; i++) {
            const int pos = start_position + direction * n + direction * i;
            if (pos < 1) {
                break;
            }
            if (pos > idx) {
                break;
            }
            if (equals_char(ref_get(ref, pos), sequence[static_cast<std::size_t>(i)])) {
                m.insert(sequence[static_cast<std::size_t>(i)]);
            } else {
                mm++;
            }
            if (mm > maxmm - n / 100) {
                break;
            }
        }
        if (m.size() < 3) {
            continue;
        }
        if (mm <= maxmm - n / 100 && i >= seq_len - 2 && i >= 8 + n / 10 && mm / static_cast<double>(i) < 0.12) {
            const int lbp = start_position + direction * n - (direction < 0 ? direction : 0);
            if (mm == 0 && i == seq_len) {
                if (conf_.y) {
                    std::println(stderr, "  Findbp: {} {} {} {} {}", sequence, start_position, lbp, mm, i);
                }
                return lbp;
            }
            if (i - mm > score) {
                bp = lbp;
                score = i - mm;
            }
        }
    }
    if (conf_.y && bp != 0) {
        std::println(stderr, "  Findbp with mismatches: {} {} {} {} {}", sequence, start_position, bp, direction,
                     score);
    }
    return bp;
}

std::pair<int, std::string> Realigner::findbi(const std::string& seq, int position, int dir) const {
    const RefMap& ref = ref_;
    const int maxmm = 3;  // maximum mismatches allowed
    const int dir_ext = dir == -1 ? 1 : 0;
    int score = 0;
    int bi = 0;
    std::string ins;
    int bi2 = 0;
    const int seq_len = static_cast<int>(seq.size());
    for (int n = 6; n < seq_len; n++) {
        if (position + 6 >= chr_length_or_npe(scope_.run, chr_)) {
            break;
        }
        int mm = 0;
        int i = 0;
        std::set<char> m;
        for (i = 0; i + n < seq_len; i++) {
            if (position + dir * i - dir_ext < 1) {
                break;
            }
            if (position + dir * i - dir_ext > chr_length_or_npe(scope_.run, chr_)) {
                break;
            }
            const char c = seq[static_cast<std::size_t>(i + n)];
            if (!equals_char(ref_get(ref, position + dir * i - dir_ext), c)) {
                mm++;
            } else {
                m.insert(c);
            }
            if (mm > maxmm) {
                break;
            }
        }
        const int mnt = static_cast<int>(m.size());
        if (mnt < 2) {  // at least three different NT for overhang sequences, weeding out low complexity seq
            continue;
        }
        if ((mnt >= 3 && i + n >= seq_len - 1 && i >= 8 && mm / static_cast<double>(i) < 0.15) ||
            (mnt >= 2 && mm == 0 && i + n == seq_len && n >= 20 && i >= 8)) {
            std::string insert = util::substr(seq, 0, n);
            std::string extra;
            int ept = 0;
            while (n + ept + 1 < seq_len && (!equals_char(ref_get(ref, position + ept * dir - dir_ext),
                                                          seq[static_cast<std::size_t>(n + ept)]) ||
                                             !equals_char(ref_get(ref, position + (ept + 1) * dir - dir_ext),
                                                          seq[static_cast<std::size_t>(n + ept + 1)]))) {
                extra += seq[static_cast<std::size_t>(n + ept)];
                ept++;
            }
            const int extra_len = static_cast<int>(extra.size());
            if (dir == -1) {
                insert += extra;
                std::reverse(insert.begin(), insert.end());
                if (extra_len > 0) {
                    insert.insert(insert.size() - extra.size(), "&");
                }
                if (mm == 0 && i + n == seq_len) {
                    bi = position - 1 - extra_len;
                    ins = insert;
                    bi2 = position - 1;
                    if (extra_len == 0) {
                        auto [b, s] = adj_ins_pos(bi, ins, ref);
                        bi = b;
                        ins = s;
                        bi2 = b;
                    }
                    return {bi, ins};
                }
                if (i - mm > score) {
                    bi = position - 1 - extra_len;
                    ins = insert;
                    bi2 = position - 1;
                    score = i - mm;
                }
            } else {
                int s = -1;
                if (extra_len > 0) {
                    insert += "&" + extra;
                } else {
                    while (s >= -n &&
                           equals_char(ref_get(ref, position + s), as_char(util::utils_char_at(insert, s)))) {
                        s--;
                    }
                    if (s < -1) {
                        const std::string tins = util::substr(insert, s + 1, 1 - s);
                        insert.erase(insert.size() + s + 1);
                        insert.insert(0, tins);
                    }
                }
                if (mm == 0 && i + n == seq_len) {
                    bi = position + s;
                    ins = insert;
                    bi2 = position + s + extra_len;
                    if (extra_len == 0) {
                        auto [b, t] = adj_ins_pos(bi, ins, ref);
                        bi = b;
                        ins = t;
                        bi2 = b;
                    }
                    return {bi, ins};
                }
                if (i - mm > score) {
                    bi = position + s;
                    ins = insert;
                    bi2 = position + s + extra_len;
                    score = i - mm;
                }
            }
        }
    }
    if (bi2 == bi && !ins.empty() && bi != 0) {
        auto [b, s] = adj_ins_pos(bi, ins, ref);
        bi = b;
        ins = s;
    }
    return {bi, ins};
}

void Realigner::adj_ref_cnt(const Variation& tv, Variation* ref, int len) const {
    if (ref == nullptr) {
        return;
    }
    if (conf_.y) {
        const std::string ref_cnt = ref->vars_count != 0 ? i2s(ref->vars_count) : "NA";
        std::println(stderr, "    AdjRefCnt: '+' {} {} {} {} Ref: {}", ref->vars_count, tv.vars_count,
                     ref->get_dir(false), tv.get_dir(false), ref_cnt);
        std::println(stderr, "    AdjRefCnt: '-' {} {} {} {} Ref: {}", ref->vars_count, tv.vars_count,
                     ref->get_dir(true), tv.get_dir(true), ref_cnt);
    }
    // the adjustment factor
    double f =
        tv.mean_position != 0 ? (tv.mean_position / tv.vars_count - len + 1) / (tv.mean_position / tv.vars_count) : 0;
    if (f < 0) {
        return;
    }
    if (f > 1) {
        f = 1;
    }
    ref->vars_count -= util::d2i_cast(f * tv.vars_count);
    ref->high_quality_reads_count -= util::d2i_cast(f * tv.high_quality_reads_count);
    ref->low_quality_reads_count -= util::d2i_cast(f * tv.low_quality_reads_count);
    ref->mean_position -= f * tv.mean_position;
    ref->mean_quality -= f * tv.mean_quality;
    ref->mean_mapping_quality -= f * tv.mean_mapping_quality;
    ref->number_of_mismatches -= f * tv.number_of_mismatches;
    ref->sub_dir(true, util::d2i_cast(f * tv.get_dir(true)));
    ref->sub_dir(false, util::d2i_cast(f * tv.get_dir(false)));
    correct_cnt(*ref);
}

void Realigner::adj_ref_factor(Variation* ref, double factor_f) const {
    if (ref == nullptr) {
        return;
    }
    if (factor_f > 1) {
        factor_f = 1;
    }
    if (factor_f < -1) {
        return;
    }
    if (conf_.y) {
        std::println(stderr, "    AdjRefFactor: {} {}", ref->vars_count, util::double_to_string(factor_f));
    }
    const int old_vars_count = ref->vars_count;
    ref->vars_count -= util::d2i_cast(factor_f * ref->vars_count);
    ref->high_quality_reads_count -= util::d2i_cast(factor_f * ref->high_quality_reads_count);
    ref->low_quality_reads_count -= util::d2i_cast(factor_f * ref->low_quality_reads_count);

    // Adjust mean mapping quality, mean quality and mean position only on number of changed counts
    double factor_cnt =
        old_vars_count != 0 ? std::abs(ref->vars_count - old_vars_count) / static_cast<double>(old_vars_count) : 1;
    // Factors must be the same sign
    if ((factor_f < 0 && factor_cnt > 0) || (factor_f > 0 && factor_cnt < 0)) {
        factor_cnt = -factor_cnt;
    }
    ref->mean_position -= ref->mean_position * factor_cnt;
    ref->mean_quality -= ref->mean_quality * factor_cnt;
    ref->mean_mapping_quality -= ref->mean_mapping_quality * factor_cnt;

    ref->number_of_mismatches -= factor_f * ref->number_of_mismatches;
    ref->vars_count_on_forward -= util::d2i_cast(factor_f * ref->vars_count_on_forward);
    ref->vars_count_on_reverse -= util::d2i_cast(factor_f * ref->vars_count_on_reverse);
    correct_cnt(*ref);
}

void Realigner::add_var_factor(Variation* vref, double factor_f) {
    if (vref == nullptr) {
        return;
    }
    if (factor_f < -1) {
        return;
    }
    vref->vars_count += util::d2i_cast(factor_f * vref->vars_count);
    vref->high_quality_reads_count += util::d2i_cast(factor_f * vref->high_quality_reads_count);
    vref->low_quality_reads_count += util::d2i_cast(factor_f * vref->low_quality_reads_count);
    vref->mean_position += factor_f * vref->mean_position;
    vref->mean_quality += factor_f * vref->mean_quality;
    vref->mean_mapping_quality += factor_f * vref->mean_mapping_quality;
    vref->number_of_mismatches += factor_f * vref->number_of_mismatches;
    vref->vars_count_on_forward += util::d2i_cast(factor_f * vref->vars_count_on_forward);
    vref->vars_count_on_reverse += util::d2i_cast(factor_f * vref->vars_count_on_reverse);
}

void Realigner::rm_cnt(Variation& vref, const Variation& tv) {
    vref.vars_count -= tv.vars_count;
    vref.high_quality_reads_count -= tv.high_quality_reads_count;
    vref.low_quality_reads_count -= tv.low_quality_reads_count;
    vref.mean_position -= tv.mean_position;
    vref.mean_quality -= tv.mean_quality;
    vref.mean_mapping_quality -= tv.mean_mapping_quality;
    vref.sub_dir(true, tv.get_dir(true));
    vref.sub_dir(false, tv.get_dir(false));
    correct_cnt(vref);
}

}  // namespace

void realigner_process(VariationData& data, RegionScope& scope) {
    if (!scope.conf.disable_sv) {
        filter_all_sv_structures(data, scope);
    }
    adjust_mnp(data, scope);
    if (scope.conf.y) {
        std::println(stderr, "TIME: Start realign");
    }
    if (scope.conf.perform_local_realignment) {
        realign_indels(data, scope);
    }
}

void realign_indels(VariationData& data, RegionScope& scope) { Realigner(data, scope).realign_indels(); }

void realigndel_with_bams(VariationData& data, RegionScope& scope, const CountMap& dels) {
    Realigner(data, scope).realigndel(true, dels);
}

void extend_reference(RegionScope& scope, const Region& region, int extension) {
    if (scope.provider == nullptr) {
        util::throw_npe("referenceResource");
    }
    if (auto r =
            io::get_reference(*scope.provider, region, extension, scope.reference, scope.conf, scope.run.chr_lengths);
        !r) {
        throw util::RecoverableError(r.error().message);
    }
}

void partial_pipeline(VariationData& data, RegionScope& scope, const Region& region, int max_read_length) {
    if (scope.bam_cache == nullptr || scope.bam.empty()) {
        return;
    }
    static const util::Regex colon = util::Regex::compile(":").value();
    std::vector<std::string> bams;
    for (auto p : util::regex_split(scope.bam, colon)) {
        bams.emplace_back(p);
    }
    RegionScope sub{scope.conf,       scope.run,      region,          scope.reference, scope.splice,
                    scope.exceptions, scope.provider, scope.bam_cache, scope.bam};
    try {
        auto pre = io::RecordPreprocessor::create(bams, region, scope.conf, *scope.bam_cache);
        if (!pre) {
            return;
        }
        CigarParser parser(sub, true, max_read_length);
        parser.process_into(*pre, data);
    } catch (const util::RecoverableError&) {
    } catch (const FatalError&) {}
}

namespace {

// Cluster (Mate + cnt)
struct Cluster {
    int cnt = 0;
    Mate m;
};

Cluster check_cluster(std::vector<Mate>& mates, int rlen, const Config& conf) {
    std::stable_sort(mates.begin(), mates.end(),
                     [](const Mate& a, const Mate& b) { return a.mate_start_ms < b.mate_start_ms; });
    if (mates.empty()) {
        throw util::RecoverableError("IndexOutOfBoundsException: Index 0 out of bounds for length 0");
    }
    std::vector<Cluster> clusters;
    const auto start_cluster = [](const Mate& m) {
        Cluster c;
        c.m.mate_start_ms = m.mate_start_ms;
        c.m.mate_end_me = m.mate_end_me;
        c.m.start_s = m.start_s;
        c.m.end_e = m.end_e;
        return c;
    };
    clusters.push_back(start_cluster(mates[0]));
    std::size_t cur = 0;
    for (const Mate& mate_m : mates) {
        if (mate_m.mate_start_ms - clusters[cur].m.mate_end_me > Config::kMinSvcDist * rlen) {
            cur++;
            clusters.insert(clusters.begin() + static_cast<long>(cur), start_cluster(mate_m));
        }
        Cluster& c = clusters[cur];
        c.cnt++;
        c.m.mate_length_mlen += mate_m.mate_length_mlen;
        if (mate_m.mate_end_me > c.m.mate_end_me) {
            c.m.mate_end_me = mate_m.mate_end_me;
        }
        if (mate_m.start_s < c.m.start_s) {
            c.m.start_s = mate_m.start_s;
        }
        if (mate_m.end_e > c.m.end_e) {
            c.m.end_e = mate_m.end_e;
        }
        c.m.pmean_rp += mate_m.pmean_rp;
        c.m.qmean_q += mate_m.qmean_q;
        c.m.qmean_big_q += mate_m.qmean_big_q;
        c.m.nm += mate_m.nm;
    }
    std::stable_sort(clusters.begin(), clusters.end(),
                     [](const Cluster& a, const Cluster& b) { return a.cnt > b.cnt; });
    if (conf.y) {
        std::string line = "Clusters; ";
        for (const Cluster& c : clusters) {
            line += std::to_string(c.cnt) + "; " + std::to_string(c.m.start_s) + "; " + std::to_string(c.m.end_e) +
                    "; " + std::to_string(c.m.mate_start_ms) + "; " + std::to_string(c.m.mate_end_me) + "; ";
        }
        std::println(stderr, "{}; ; out of; {}", line, mates.size());
    }
    const Cluster& first = clusters[0];
    Cluster r;
    if (first.cnt / static_cast<double>(mates.size()) >= 0.60) {
        r = first;
        r.m.mate_length_mlen = first.m.mate_length_mlen / first.cnt;
    }
    return r;
}

// LinkedHashMap entries sorted by value, highest first (a stable sort); the first key or 0.
int top_soft(const util::LinkedMap<int, int>& soft) {
    std::vector<std::pair<int, int>> v;
    soft.for_each([&](int k, int c) { v.emplace_back(k, c); });
    std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    return v.empty() ? 0 : v[0].first;
}

void filter_sv(std::vector<Sclip>& sv_list, VariationData& data, const RegionScope& scope) {
    for (Sclip& sv : sv_list) {
        try {
            const Cluster cluster = check_cluster(sv.mates, data.max_read_length, scope.conf);
            if (cluster.m.mate_start_ms != 0) {
                sv.mstart = cluster.m.mate_start_ms;
                sv.mend = cluster.m.mate_end_me;
                sv.vars_count = cluster.cnt;
                sv.mlen = cluster.m.mate_length_mlen;
                sv.start = cluster.m.start_s;
                sv.end = cluster.m.end_e;
                sv.mean_position = cluster.m.pmean_rp;
                sv.mean_quality = cluster.m.qmean_q;
                sv.mean_mapping_quality = cluster.m.qmean_big_q;
                sv.number_of_mismatches = cluster.m.nm;
            } else {
                sv.used = true;
            }
            // Too many unhappy mates are false positive
            if (sv.disc != 0 && sv.vars_count / static_cast<double>(sv.disc) < 0.5) {
                if (!(sv.vars_count / static_cast<double>(sv.disc) >= 0.35 && sv.vars_count >= 5)) {
                    sv.used = true;
                }
            }
            sv.softp = top_soft(sv.soft);
            if (sv.softp != 0) {
                data.softp2sv.get_or_create(sv.softp).push_back(&sv);
            }
            if (scope.conf.y) {
                std::println(stderr, "SV cluster: {} {} {} {} Cnt: {} Discordant Cnt: {} Softp: {} Used: {}",
                             cluster.m.start_s, cluster.m.end_e, cluster.m.mate_start_ms, cluster.m.mate_end_me,
                             sv.vars_count, sv.disc, sv.softp, sv.used);
            }
        } catch (const util::RecoverableError& e) {
            print_exception_and_continue(e, "cluster", std::to_string(sv.start), &scope.region, scope.exceptions);
        }
    }
}

}  // namespace

void filter_all_sv_structures(VariationData& data, const RegionScope& scope) {
    SVStructures& s = data.sv_structures;
    filter_sv(s.svfinv3, data, scope);
    filter_sv(s.svrinv3, data, scope);
    filter_sv(s.svfinv5, data, scope);
    filter_sv(s.svrinv5, data, scope);
    filter_sv(s.svfdel, data, scope);
    filter_sv(s.svrdel, data, scope);
    filter_sv(s.svfdup, data, scope);
    filter_sv(s.svrdup, data, scope);
    s.svffus.for_each([&](const std::string&, std::vector<Sclip>& v) { filter_sv(v, data, scope); });
    s.svrfus.for_each([&](const std::string&, std::vector<Sclip>& v) { filter_sv(v, data, scope); });
    data.softp2sv.for_each([](int, std::vector<Sclip*>& v) {
        std::stable_sort(v.begin(), v.end(),
                         [](const Sclip* a, const Sclip* b) { return a->vars_count > b->vars_count; });
    });
}

void adjust_mnp(VariationData& data, const RegionScope& scope) {
    const RefMap& reference = scope.reference.reference_sequences;
    auto tmp = fill_and_sort_tmp(data.mnp);
    for (const auto& tpl : tmp) {
        int last_position = 0;
        try {
            const int position = tpl.position;
            last_position = position;
            const std::string& vn = tpl.description_string;
            VariationMap* vars_on_position = data.non_insertion_variants.get(position);
            if (vars_on_position == nullptr) {
                continue;
            }
            Variation* vref = vars_on_position->get(vn);
            if (vref == nullptr) {  // The variant is likely already been used by indel realignment
                continue;
            }
            if (scope.conf.y) {
                std::println(stderr, "  AdjMnt: {} {} {}", position, vn, vref->vars_count);
            }
            const std::string mnt = util::replace_first(vn, "&", "");
            for (int i = 0; i < static_cast<int>(mnt.size()) - 1; i++) {
                std::string left = util::substr(mnt, 0, i + 1);
                if (left.size() > 1) {
                    left = insert_amp(left);
                }
                std::string right = util::substr(mnt, -(static_cast<long>(mnt.size()) - i - 1));
                if (right.size() > 1) {
                    right = insert_amp(right);
                }
                if (Variation* tref = vars_on_position->get(left)) {
                    if (tref->vars_count <= 0) {
                        continue;
                    }
                    if (tref->vars_count < vref->vars_count && tref->mean_position / tref->vars_count <= i + 1) {
                        if (scope.conf.y) {
                            std::println(stderr, "    AdjMnt Left: {} {} Left: {} Cnt: {}", position, vn, left,
                                         tref->vars_count);
                        }
                        adj_cnt(*vref, *tref);
                        vars_on_position->erase(left);
                    }
                }
                if (VariationMap* next = data.non_insertion_variants.get(position + i + 1)) {
                    if (Variation* tref = next->get(right)) {
                        if (tref->vars_count < 0) {
                            continue;
                        }
                        if (tref->vars_count < vref->vars_count) {
                            if (scope.conf.y) {
                                std::println(stderr, "    AdjMnt Right: {} {} Right: {} Cnt: {}", position, vn, right,
                                             tref->vars_count);
                            }
                            adj_cnt(*vref, *tref);
                            inc_cnt(data.ref_coverage, position, tref->vars_count);
                            next->erase(right);
                        }
                    }
                }
            }
            if (Sclip* sc3v = data.soft_clips_3end.get(position)) {
                if (!sc3v->used) {
                    const std::string seq = findconseq(*sc3v, 0, scope.run, scope.conf);
                    if (seq.starts_with(mnt)) {
                        if (seq.size() == mnt.size() ||
                            ismatchref(util::substring(seq, static_cast<long>(mnt.size())), reference,
                                       position + static_cast<int>(mnt.size()), 1)) {
                            // nonInsertionVariants.get(position).get(vn): NPE if vn was removed above
                            Variation* target = data.non_insertion_variants.get(position)->get(vn);
                            if (target == nullptr) {
                                util::throw_npe("adjustMNP: " + vn + " removed");
                            }
                            adj_cnt(*target, *sc3v);
                            inc_cnt(data.ref_coverage, position, sc3v->vars_count);
                            sc3v->used = true;
                        }
                    }
                }
            }
            if (Sclip* sc5v = data.soft_clips_5end.get(position + static_cast<int>(mnt.size()))) {
                if (!sc5v->used) {
                    std::string seq = findconseq(*sc5v, 0, scope.run, scope.conf);
                    if (!seq.empty() && seq.size() >= mnt.size()) {
                        std::reverse(seq.begin(), seq.end());
                        if (seq.ends_with(mnt)) {
                            if (seq.size() == mnt.size() ||
                                ismatchref(util::substring(seq, 0, static_cast<long>(seq.size() - mnt.size())),
                                           reference, position - 1, -1)) {
                                Variation* target = data.non_insertion_variants.get(position)->get(vn);
                                if (target == nullptr) {
                                    util::throw_npe("adjustMNP: " + vn + " removed");
                                }
                                adj_cnt(*target, *sc5v);
                                inc_cnt(data.ref_coverage, position, sc5v->vars_count);
                                sc5v->used = true;
                            }
                        }
                    }
                }
            }
        } catch (const util::RecoverableError& e) {
            print_exception_and_continue(e, "MNP", std::to_string(last_position), &scope.region, scope.exceptions);
        }
    }
}

void adj_snv(VariationData& data, const RegionScope& scope) {
    const RefMap& reference = scope.reference.reference_sequences;
    data.soft_clips_5end.for_each([&](int position, Sclip& sclip) {
        if (sclip.used) {
            return;
        }
        const std::string seq = findconseq(sclip, 0, scope.run, scope.conf);
        if (seq.size() > 5) {
            return;
        }
        const std::string bp = util::substr(seq, 0, 1);
        const int previous_position = position - 1;
        VariationMap* prev = data.non_insertion_variants.get(previous_position);
        if (prev != nullptr && prev->contains(bp)) {
            if (seq.size() > 1) {
                if (is_not_equals(ref_get(reference, position - 2), seq[1])) {
                    return;
                }
            }
            adj_cnt(*prev->get(bp), sclip);
            inc_cnt(data.ref_coverage, previous_position, sclip.vars_count);
        }
    });
    data.soft_clips_3end.for_each([&](int position, Sclip& sclip) {
        if (sclip.used) {
            return;
        }
        const std::string seq = findconseq(sclip, 0, scope.run, scope.conf);
        if (seq.size() > 5) {
            return;
        }
        const std::string bp = util::substr(seq, 0, 1);
        VariationMap* here = data.non_insertion_variants.get(position);
        if (here != nullptr && here->contains(bp)) {
            if (seq.size() > 1) {
                if (is_not_equals(ref_get(reference, position + 1), seq[1])) {
                    return;
                }
            }
            adj_cnt(*here->get(bp), sclip);
            inc_cnt(data.ref_coverage, position, sclip.vars_count);
        }
    });
}

Match find_match(std::string seq, const io::Reference& reference, int position, int dir, int seed, int mm,
                 const Config& conf) {
    const RefMap& ref = reference.reference_sequences;
    if (dir == -1) {
        seq = util::reverse(seq);  // dir==-1 means 5' clip
    }
    if (conf.y) {
        std::println(stderr, "    Working Match {} {} {} SEED: {}", position, seq, dir, seed);
    }
    std::string extra;
    const int seq_len = static_cast<int>(seq.size());
    for (int i = seq_len - seed; i >= 0; i--) {
        const std::string s = util::substr(seq, i, seed);
        auto it = reference.seed.find(s);
        if (it == reference.seed.end()) {
            continue;
        }
        const std::vector<int>& seeds = it->second;
        if (seeds.size() != 1) {
            continue;
        }
        const int first_seed = seeds[0];
        int bp = dir == 1 ? first_seed - i : first_seed + seq_len - i - 1;
        if (ismatchref(seq, ref, bp, dir, mm)) {
            int m = dir == -1 ? -1 : 0;
            while (has_and_not_equals(util::utils_char_at(seq, m), ref, bp)) {
                extra += util::substr(seq, m, 1);
                bp += dir;
                m += dir;
            }
            if (!extra.empty() && dir == -1) {
                extra = util::reverse(extra);
            }
            if (conf.y) {
                std::println(stderr, "      Found SV BP: {} BP: {} SEEDpos {} {} {} {} {} extra: {}", dir, bp,
                             first_seed, position, s, i, seq, extra);
            }
            return {bp, extra};
        }
        // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
        std::string sseq = seq;
        int eqcnt = 0;
        for (int ii = 1; ii <= 15; ii++) {
            bp += dir;
            sseq = dir == 1 ? util::substr(sseq, 1) : util::substr(sseq, 0, -1);
            if (dir == 1) {
                if (has_and_not_equals(util::utils_char_at(sseq, 0), ref, bp)) {
                    continue;
                }
                eqcnt++;
                if (has_and_not_equals(util::utils_char_at(sseq, 1), ref, bp + 1)) {
                    continue;
                }
                extra = util::substr(seq, 0, ii);
            } else {
                if (has_and_not_equals(util::utils_char_at(sseq, -1), ref, bp)) {
                    continue;
                }
                eqcnt++;
                if (has_and_not_equals(util::utils_char_at(sseq, -2), ref, bp - 1)) {
                    continue;
                }
                extra = util::substr(seq, -ii);
            }
            if (eqcnt >= 3 && eqcnt / static_cast<double>(ii) > 0.5) {
                break;
            }
            if (conf.y) {
                std::println(stderr, "      FoundSEED SV BP: {} BP: {} SEEDpos{} {} {} {} {} EXTRA: {}", dir, bp,
                             first_seed, position, s, i, seq, extra);
            }
            if (ismatchref(sseq, ref, bp, dir, 1)) {
                return {bp, extra};
            }
        }
    }
    return {0, ""};
}

bool is_overlap(int start1, int end1, int start2, int end2, int rlen) {
    if (start1 >= end2 || start2 >= end1) {
        return false;
    }
    int positions[] = {start1, end1, start2, end2};
    std::sort(std::begin(positions), std::end(positions));
    const int ins = positions[2] - positions[1];
    if (end1 != start1 && end2 != start2 && ins / static_cast<double>(end1 - start1) > 0.75 &&
        ins / static_cast<double>(end2 - start2) > 0.75) {
        return true;
    }
    return positions[1] - positions[0] + positions[3] - positions[2] < 3 * rlen;
}

namespace {

SvMark mark(int start, int end, const std::vector<std::vector<Sclip>*>& svs, int rlen, const Config& conf, bool dup) {
    SvMark r;
    for (std::vector<Sclip>* list : svs) {
        for (Sclip& sv_r : *list) {
            int start2 = 0;
            int end2 = 0;
            if (dup) {
                if (sv_r.start < sv_r.mstart) {
                    start2 = sv_r.start;
                    end2 = sv_r.mend;
                } else {
                    start2 = sv_r.mstart;
                    end2 = sv_r.end;
                }
            } else {
                if (sv_r.start < sv_r.mstart) {
                    start2 = sv_r.end;
                    end2 = sv_r.mstart;
                } else {
                    start2 = sv_r.mend;
                    end2 = sv_r.start;
                }
            }
            if (conf.y) {
                std::println(stderr, "   Marking {}SV {} {} {} {} cnt: {}", dup ? "DUP " : "", start, end, start2, end2,
                             sv_r.vars_count);
            }
            if (is_overlap(start, end, start2, end2, rlen)) {
                if (conf.y) {
                    std::println(stderr, "       {}SV {} {} {} {} cnt: {} marked", dup ? "DUP " : "", start, end,
                                 start2, end2, sv_r.vars_count);
                }
                sv_r.used = true;
                r.count++;
                r.pairs += sv_r.vars_count;
                if (sv_r.end != sv_r.start) {
                    r.cov += (sv_r.vars_count * rlen) / (sv_r.end - sv_r.start) + 1;
                }
            }
        }
    }
    return r;
}

}  // namespace

SvMark mark_sv(int start, int end, const std::vector<std::vector<Sclip>*>& svs, int rlen, const Config& conf) {
    return mark(start, end, svs, rlen, conf, false);
}

SvMark mark_dup_sv(int start, int end, const std::vector<std::vector<Sclip>*>& svs, int rlen, const Config& conf) {
    return mark(start, end, svs, rlen, conf, true);
}

SvInfo& get_sv(PositionMap& hash, int start) {
    VariationMap& map = hash.get_or_create(start);
    if (!map.sv) {
        map.sv = SvInfo{};
        map.put("SV", Variation{});
    }
    if (!map.contains("SV")) {
        map.put("SV", Variation{});
    }
    return *map.sv;
}

void remove_sv(PositionMap& hash, int start) {
    VariationMap* map = hash.get(start);
    if (map == nullptr) {
        util::throw_npe("removeSV: hash.get(start)");
    }
    map->sv.reset();
    map->erase("SV");
}

std::string join_ref_double(const RefMap& ref, int from, double to) {
    std::string sb;
    for (int i = from; i < to; i++) {
        if (auto c = ref_get(ref, i)) {
            sb += *c;
        }
    }
    return sb;
}

std::string join_ref_for_5_lgins(const RefMap& ref, int from, int to, const std::string& seq,
                                 const std::string& extra) {
    std::string sb;
    const int n = static_cast<int>(seq.size()) - static_cast<int>(extra.size());
    for (int i = from; i <= to; i++) {
        if (to - i < n) {
            sb += as_char(util::utils_char_at(seq, to - i + static_cast<int>(extra.size())));
        } else if (auto c = ref_get(ref, i)) {
            sb += *c;
        }
    }
    return sb;
}

std::string join_ref_for_3_lgins(const RefMap& ref, int from, int to, int shift5, const std::string& seq,
                                 const std::string& extra) {
    std::string sb;
    const int n = static_cast<int>(seq.size()) - static_cast<int>(extra.size());
    for (int i = from; i <= to; i++) {
        if (i - from >= shift5 && i - from - shift5 < n) {
            sb += as_char(util::utils_char_at(seq, i - from - shift5 + static_cast<int>(extra.size())));
        } else if (auto c = ref_get(ref, i)) {
            sb += *c;
        }
    }
    return sb;
}

bool ismatch(const std::string& seq1, const std::string& seq2_in, int dir, const Config& conf, int mm_limit) {
    if (conf.y) {
        std::println(stderr, "    Matching two seqs {} {} {} {}", seq1, seq2_in, dir, mm_limit);
    }
    const std::string seq2 = strip_hash_caret(seq2_in);
    int mm = 0;
    for (int n = 0; n < static_cast<int>(seq1.size()) && n < static_cast<int>(seq2.size()); n++) {
        const std::string c = util::substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1);
        if (seq1[static_cast<std::size_t>(n)] != util::char_at(c, 0)) {
            mm++;
        }
    }
    return mm <= mm_limit && mm / static_cast<double>(seq1.size()) < 0.15;
}

std::pair<int, std::string> adj_ins_pos(int bi, std::string ins, const RefMap& ref) {
    int n = 1;
    const int len = static_cast<int>(ins.size());
    while (equals_char(ref_get(ref, bi), util::char_at(ins, len - n))) {
        n++;
        if (n > len) {
            n = 1;
        }
        bi--;
    }
    if (n > 1) {
        ins = util::substr(ins, 1 - n) + util::substr(ins, 0, 1 - n);
    }
    return {bi, ins};
}

Match35 find35match(const std::string& seq5, const std::string& seq3, const Config& conf) {
    const int long_mismatch = 2;
    int max_matched_length = 0;
    int b3 = 0;
    int b5 = 0;
    const int len5 = static_cast<int>(seq5.size());
    const int len3 = static_cast<int>(seq3.size());
    for (int i = 0; i < len5 - 8; i++) {
        for (int j = 1; j < len3 - 8; j++) {
            int number_of_mismatch = 0;
            int total_length = 0;
            while (total_length + j <= len3 && i + total_length <= len5) {
                if (util::substr(seq3, -j - total_length, 1) != util::substr(seq5, i + total_length, 1)) {
                    number_of_mismatch++;
                }
                if (number_of_mismatch > long_mismatch) {
                    break;
                }
                total_length++;
            }
            if (total_length - number_of_mismatch > max_matched_length && total_length - number_of_mismatch > 8 &&
                number_of_mismatch / static_cast<double>(total_length) < 0.1 &&
                (total_length + j >= len3 || i + total_length >= len5)) {
                max_matched_length = total_length - number_of_mismatch;
                b3 = j;
                b5 = i;
                if (conf.y) {
                    std::println(stderr, "      Found 35 Match, {} {} {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", seq5, seq3,
                                 total_length + j, len3, i + total_length, len5, total_length, number_of_mismatch, i,
                                 j);
                }
                return {b5, b3, max_matched_length};
            }
        }
    }
    return {b5, b3, max_matched_length};
}

}  // namespace omnivar::core
