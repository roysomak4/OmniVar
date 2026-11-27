#include "core/cigar_modifier.hpp"

#include <cstdio>
#include <cstdlib>
#include <print>
#include <set>
#include <string_view>

#include "core/variation_utils.hpp"
#include "util/checked_string.hpp"
#include "util/number_parse.hpp"
#include "util/regex_engine.hpp"
#include "util/sequence.hpp"

namespace senri::core {

namespace {

using util::Regex;

// A compiled Patterns.* regex. the upstream regex libraries agree on every pattern used here
// (greedy/lazy quantifiers, \d, \D, classes, anchors), so one engine serves both.
struct Rx {
    explicit Rx(std::string_view pattern) : re(Regex::compile(pattern).value()) {}
    Regex re;
};

// Patterns (CIGAR patterns)
const Rx& rx(int id) {
    static const Rx all[] = {
        Rx("^(\\d+)D"),                                                 // 0 BEGIN_NUMBER_D
        Rx("(\\d+)D$"),                                                 // 1 END_NUMBER_D
        Rx("^(\\d+)I"),                                                 // 2 BEGIN_NUMBER_I
        Rx("(\\d+)I$"),                                                 // 3 END_NUMBER_I
        Rx("^(\\d\\d+)S"),                                              // 4 SA_CIGAR_D_S_5clip_GROUP
        Rx("(\\d\\d+)S$"),                                              // 5 SA_CIGAR_D_S_3clip_GROUP
        Rx("^\\d+S"),                                                   // 6 SA_CIGAR_D_S_5clip_GROUP_Repl
        Rx("\\d\\d+S$"),                                                // 7 SA_CIGAR_D_S_3clip_GROUP_Repl
        Rx("^(\\d+)S(\\d+)([ID])"),                                     // 8 BEGIN_NUMBER_S_NUMBER_IorD
        Rx("(\\d+)([ID])(\\d+)S$"),                                     // 9 NUMBER_IorD_NUMBER_S_END
        Rx("^(\\d+)S(\\d+)M(\\d+)([ID])"),                              // 10 BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD
        Rx("(\\d+)([ID])(\\d+)M(\\d+)S$"),                              // 11 NUMBER_IorD_NUMBER_M_NUMBER_S_END
        Rx("^(\\d)M(\\d+)([ID])(\\d+)M"),                               // 12 BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M
        Rx("^\\dM\\d+[ID]\\d+M"),                                       // 13 BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_
        Rx("(\\d+)([ID])(\\d)M$"),                                      // 14 NUMBER_IorD_DIGIT_M_END
        Rx("(\\d+)([ID])(\\d+)M$"),                                     // 15 NUMBER_IorD_NUMBER_M_END
        Rx("^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"),  // 16 D_M_D_DD_M_D_I_D_M_D_DD
        Rx("(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"),        // 17 D_M_D_DD_M_D_I_D_M_D_DD_prim
        Rx("^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M"),  // 18 threeDeletionsPattern
        Rx("^(.*?)(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M"),  // 19 threeIndelsPattern
        Rx("\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M"),  // 20 DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM
        Rx("\\d+M\\d+D\\d+M\\d+D\\d+M\\d+D\\d+M"),           // 21 DM_DD_DM_DD_DM_DD_DM
        Rx("(\\d+)D(\\d+)M(\\d+)([DI])(\\d+I)?"),            // 22 DIG_D_DIG_M_DIG_DI_DIGI
        Rx("(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"),            // 23 DIG_I_DIG_M_DIG_DI_DIGI
        Rx("(\\D)(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"),       // 24 NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI
        Rx("(\\d+)D(\\d+)D"),                                // 25 DIG_D_DIG_D
        Rx("(\\d+)I(\\d+)I"),                                // 26 DIG_I_DIG_I
        Rx("^(.*?)(\\d+)M$"),                                // 27 BEGIN_ANY_DIG_M_END
        Rx("\\d+M$"),                                        // 28 DIG_M_END
        Rx("^(\\d+)M"),                                      // 29 BEGIN_DIG_M
        Rx("^(\\d+)S(\\d+)M"),                               // 30 DIG_S_DIG_M
        Rx("\\d+M\\d+S$"),                                   // 31 DIG_M_DIG_S_END
        Rx("^(.*?)(\\d+)M(\\d+)S$"),                         // 32 ANY_NUMBER_M_NUMBER_S_END
        Rx("(\\d+)[MND]"),                                   // 33 ALIGNED_LENGTH_MND
        Rx("(\\d+)[MIS]"),                                   // 34 SOFT_CLIPPED
    };
    return all[id];
}

enum : int {
    BEGIN_NUMBER_D,
    END_NUMBER_D,
    BEGIN_NUMBER_I,
    END_NUMBER_I,
    SA_5CLIP_GROUP,
    SA_3CLIP_GROUP,
    SA_5CLIP_REPL,
    SA_3CLIP_REPL,
    BEGIN_NUMBER_S_NUMBER_IorD,
    NUMBER_IorD_NUMBER_S_END,
    BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD,
    NUMBER_IorD_NUMBER_M_NUMBER_S_END,
    BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M,
    BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_,
    NUMBER_IorD_DIGIT_M_END,
    NUMBER_IorD_NUMBER_M_END,
    D_M_D_DD_M_D_I_D_M_D_DD,
    D_M_D_DD_M_D_I_D_M_D_DD_prim,
    THREE_DELETIONS,
    THREE_INDELS,
    DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM,
    DM_DD_DM_DD_DM_DD_DM,
    DIG_D_DIG_M_DIG_DI_DIGI,
    DIG_I_DIG_M_DIG_DI_DIGI,
    NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI,
    DIG_D_DIG_D,
    DIG_I_DIG_I,
    BEGIN_ANY_DIG_M_END,
    DIG_M_END,
    BEGIN_DIG_M,
    DIG_S_DIG_M,
    DIG_M_DIG_S_END,
    ANY_NUMBER_M_NUMBER_S_END,
    ALIGNED_LENGTH_MND,
    SOFT_CLIPPED,
};

// A successful find on a fixed text.
class Found {
public:
    Found(std::string text, Regex::Match m) : text_(std::move(text)), m_(std::move(m)) {}

    // Matcher.group(i); nullopt for a group that did not take part.
    [[nodiscard]] std::optional<std::string> group(std::size_t i) const {
        auto g = m_.group(text_, i);
        if (!g) {
            return std::nullopt;
        }
        return std::string(*g);
    }
    // toInt(matcher.group(i)): Integer.parseInt (NPE on a missing group).
    [[nodiscard]] int num(std::size_t i) const {
        auto g = group(i);
        if (!g) {
            util::throw_npe("regex group " + std::to_string(i));
        }
        auto v = util::parse_int(*g);
        if (!v) {
            throw util::RecoverableError("NumberFormatException: For input string: \"" + *g + "\"");
        }
        return *v;
    }
    [[nodiscard]] std::string str(std::size_t i) const { return group(i).value_or(""); }
    [[nodiscard]] bool is(std::size_t i, std::string_view s) const {
        auto g = group(i);
        return g && *g == s;
    }
    // Matcher.replaceFirst on the text the matcher was created with.
    [[nodiscard]] std::string replace(std::string_view with) const {
        return text_.substr(0, m_.begin) + std::string(with) + text_.substr(m_.end);
    }

private:
    std::string text_;
    Regex::Match m_;
};

std::optional<Found> find(int id, const std::string& text) {
    auto m = rx(id).re.find(text);
    if (!m) {
        return std::nullopt;
    }
    return Found(text, std::move(*m));
}

// Pattern.matcher(text).replaceFirst(with) (also jregex Replacer.replace for the anchored
// patterns used with it, which can match only once).
std::string replace_first(int id, const std::string& text, std::string_view with) {
    auto f = find(id, text);
    return f ? f->replace(with) : text;
}

// Utils.sum(Utils.globalFind(pattern, string))
int sum_global_find(int id, const std::string& text) {
    int result = 0;
    std::size_t from = 0;
    while (from <= text.size()) {
        auto m = rx(id).re.find(text, from);
        if (!m) {
            break;
        }
        Found f(text, *m);
        result += f.num(1);
        from = m->end > m->begin ? m->end : m->end + 1;
    }
    return result;
}

std::string i2s(int v) { return std::to_string(v); }

class Modifier {
public:
    Modifier(int position, const std::string& cigar, const std::string& query_sequence,
             const std::string& query_quality, const io::Reference& ref, int indel, int max_read_length,
             const Config& conf)
        : position_(position),
          cigar_str_(cigar),
          query_sequence_(query_sequence),
          query_quality_(query_quality),
          ref_(ref),
          reference_(ref.reference_sequences),
          indel_(indel),
          max_read_length_(max_read_length),
          conf_(conf) {}

    void run();

    int position_;
    std::string cigar_str_;
    std::string query_sequence_;
    std::string query_quality_;

private:
    void combine_begin_dig_m(const Found& matcher);
    void combine_dig_s_dig_m(const Found& matcher);
    void capture_mis_softly_3_mismatches(const Found& matcher);
    void capture_mis_softly_ms(const Found& matcher);
    bool combine_to_close_to_one(const Found& matcher, bool flag);
    bool combine_to_close_to_correct(const Found& matcher, bool flag);
    bool three_indels(const Found& matcher, bool flag);
    bool three_deletions(const Found& matcher, bool flag);
    bool two_deletions_insertion_to_complex(const Found& matcher, bool flag);
    bool begin_digit_m_number_i_or_d_number_m(const Found& matcher);

    // HashSet<Character>.add(reference.get(pos)): null is a distinct element.
    void add_ref(std::set<int>& set, int pos) const {
        auto r = ref_get(reference_, pos);
        set.insert(r ? static_cast<unsigned char>(*r) : -1);
    }

    const io::Reference& ref_;
    const RefMap& reference_;
    int indel_;
    int max_read_length_;
    const Config& conf_;
};

void Modifier::run() {
    // flag is set to true if CIGAR string is modified and should be looked at again
    bool flag = true;
    // if CIGAR starts with deletion cut it off
    if (auto mm = find(BEGIN_NUMBER_D, cigar_str_)) {
        position_ += mm->num(1);
        cigar_str_ = mm->replace("");
    }
    if (auto mm = find(END_NUMBER_D, cigar_str_)) {
        cigar_str_ = mm->replace("");
    }
    // replace insertion at the beginning and end with soft clipping
    if (auto mm = find(BEGIN_NUMBER_I, cigar_str_)) {
        cigar_str_ = mm->replace(i2s(mm->num(1)) + "S");
    }
    if (auto mm = find(END_NUMBER_I, cigar_str_)) {
        cigar_str_ = mm->replace(i2s(mm->num(1)) + "S");
    }
    auto sc5_mm = find(SA_5CLIP_GROUP, cigar_str_);
    auto sc3_mm = find(SA_3CLIP_GROUP, cigar_str_);
    const auto& reference_seed_map = ref_.seed;
    if (sc5_mm) {
        const int cigar_element = sc5_mm->num(1);
        if (!conf_.chimeric && cigar_element >= Config::kSeed2) {
            const std::string sseq = util::substr(query_sequence_, 0, cigar_element);
            const std::string sequence = util::complement(util::reverse(sseq));
            const std::string reverse_complemented_seed = util::substring(sequence, 0, Config::kSeed2);
            auto it = reference_seed_map.find(reverse_complemented_seed);
            if (it != reference_seed_map.end()) {
                const auto& positions = it->second;
                if (positions.size() == 1 && std::abs(position_ - positions[0]) < 2 * max_read_length_) {
                    cigar_str_ = replace_first(SA_5CLIP_REPL, cigar_str_, "");  // $a[5] = ~s / ^\d + S//;
                    query_sequence_ = util::substr(query_sequence_, cigar_element);
                    query_quality_ = util::substr(query_quality_, cigar_element);
                    if (conf_.y) {
                        std::println(stderr, "{} at 5' is a chimeric at {} by SEED {}", sequence, position_,
                                     Config::kSeed2);
                    }
                }
            }
        }
    } else if (sc3_mm) {
        const int cigar_element = sc3_mm->num(1);
        if (!conf_.chimeric && cigar_element >= Config::kSeed2) {
            const std::string sseq = util::substr(query_sequence_, -cigar_element, cigar_element);
            const std::string sequence = util::complement(util::reverse(sseq));
            const std::string reverse_complemented_seed = util::substr(sequence, -Config::kSeed2, Config::kSeed2);
            auto it = reference_seed_map.find(reverse_complemented_seed);
            if (it != reference_seed_map.end()) {
                const auto& positions = it->second;
                if (positions.size() == 1 && std::abs(position_ - positions[0]) < 2 * max_read_length_) {
                    cigar_str_ = replace_first(SA_3CLIP_REPL, cigar_str_, "");  // $a[5] = ~s /\d\d + S$//;
                    query_sequence_ =
                        util::substr(query_sequence_, 0, static_cast<long>(query_sequence_.size()) - cigar_element);
                    query_quality_ =
                        util::substr(query_quality_, 0, static_cast<long>(query_quality_.size()) - cigar_element);
                    if (conf_.y) {
                        std::println(stderr, "{} at 3' is a chimeric at {} by SEED {}", sequence, position_,
                                     Config::kSeed2);
                    }
                }
            }
        }
    }

    while (flag && indel_ > 0) {
        flag = false;
        if (auto mm = find(BEGIN_NUMBER_S_NUMBER_IorD, cigar_str_)) {
            // If CIGAR starts with soft-clipping followed by insertion or deletion: add the
            // inserted sequence to the soft-clipped start, or move position past the deletion.
            const std::string tslen = i2s(mm->num(1) + (mm->is(3, "I") ? mm->num(2) : 0)) + "S";
            position_ += mm->is(3, "D") ? mm->num(2) : 0;
            cigar_str_ = mm->replace(tslen);
            flag = true;
        }
        if (auto mm = find(NUMBER_IorD_NUMBER_S_END, cigar_str_)) {
            // If CIGAR ends with insertion or deletion followed by soft-clipping
            const std::string tslen = i2s(mm->num(3) + (mm->is(2, "I") ? mm->num(1) : 0)) + "S";
            cigar_str_ = mm->replace(tslen);
            flag = true;
        }
        if (auto mm = find(BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD, cigar_str_)) {
            const int tmid = mm->num(2);
            if (tmid <= 10) {  // matched part no longer than 10: everything becomes soft-clipping
                const std::string tslen = i2s(mm->num(1) + tmid + (mm->is(4, "I") ? mm->num(3) : 0)) + "S";
                position_ += tmid + (mm->is(4, "D") ? mm->num(3) : 0);
                cigar_str_ = mm->replace(tslen);
                flag = true;
            }
        }
        if (auto mm = find(NUMBER_IorD_NUMBER_M_NUMBER_S_END, cigar_str_)) {
            const int tmid = mm->num(3);
            if (tmid <= 10) {
                const std::string tslen = i2s(mm->num(4) + tmid + (mm->is(2, "I") ? mm->num(1) : 0)) + "S";
                cigar_str_ = mm->replace(tslen);
                flag = true;
            }
        }

        // The following two clauses make indels at the end of reads soft clipped and let
        // the algorithm identify the indels.
        if (auto mm = find(BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M, cigar_str_)) {
            flag = begin_digit_m_number_i_or_d_number_m(*mm);
        }
        if (auto mm = find(NUMBER_IorD_DIGIT_M_END, cigar_str_)) {
            const int tmid = mm->num(3);
            const std::string tslen = i2s(tmid + (mm->is(2, "I") ? mm->num(1) : 0)) + "S";
            cigar_str_ = replace_first(NUMBER_IorD_NUMBER_M_END, cigar_str_, tslen);
            flag = true;
        }

        // Combine two deletions and insertion into one complex if they are close
        auto mm = find(D_M_D_DD_M_D_I_D_M_D_DD, cigar_str_);
        if (mm) {
            flag = two_deletions_insertion_to_complex(*mm, flag);
        } else if (auto three_deletions_matcher = find(THREE_DELETIONS, cigar_str_)) {
            flag = three_deletions(*three_deletions_matcher, flag);
        } else if (auto three_indels_matcher = find(THREE_INDELS, cigar_str_)) {
            flag = three_indels(*three_indels_matcher, flag);
        }

        if (auto cm = find(DIG_D_DIG_M_DIG_DI_DIGI, cigar_str_)) {
            flag = combine_to_close_to_correct(*cm, flag);
        }
        if (auto cm = find(NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI, cigar_str_); cm && !cm->is(1, "D") && !cm->is(1, "H")) {
            flag = combine_to_close_to_one(*cm, flag);
        }
        if (auto cm = find(DIG_D_DIG_D, cigar_str_)) {
            const int dlen = cm->num(1) + cm->num(2);
            cigar_str_ = cm->replace(i2s(dlen) + "D");
            flag = true;
        }
        if (auto cm = find(DIG_I_DIG_I, cigar_str_)) {
            const int ilen = cm->num(1) + cm->num(2);
            cigar_str_ = cm->replace(i2s(ilen) + "I");
            flag = true;
        }
    }

    if (auto mtch = find(ANY_NUMBER_M_NUMBER_S_END, cigar_str_)) {
        capture_mis_softly_ms(*mtch);
    } else if (auto mtch3 = find(BEGIN_ANY_DIG_M_END, cigar_str_)) {
        capture_mis_softly_3_mismatches(*mtch3);
    }

    if (auto mtch = find(DIG_S_DIG_M, cigar_str_)) {
        combine_dig_s_dig_m(*mtch);
    } else if (auto mtch2 = find(BEGIN_DIG_M, cigar_str_)) {
        combine_begin_dig_m(*mtch2);
    }
}

void Modifier::combine_begin_dig_m(const Found& matcher) {
    int mch = matcher.num(1);
    int rn = 0;
    int rrn = 0;
    int rmch = 0;
    while (rrn < mch && rn < mch) {
        if (!reference_.contains(position_ + rrn)) {
            break;
        }
        if (is_has_and_not_equals(reference_, position_ + rrn, query_sequence_, rrn)) {
            rn = rrn + 1;
            rmch = 0;
        } else if (is_has_and_equals(reference_, position_ + rrn, query_sequence_, rrn)) {
            rmch++;
        }
        rrn++;
        if (rmch >= 3) {
            break;
        }
    }
    if (rn > 0 && rn <= 3) {
        mch -= rn;
        cigar_str_ = matcher.replace(i2s(rn) + "S" + i2s(mch) + "M");
        position_ += rn;
    }
}

void Modifier::combine_dig_s_dig_m(const Found& matcher) {
    int mch = matcher.num(2);   // length of matched sequence
    int soft = matcher.num(1);  // length of soft-clipping
    // number of bases before the matched sequence that match in reference and read
    int rn = 0;
    std::set<int> rn_set;
    while (rn < soft && is_has_and_equals(reference_, position_ - rn - 1, query_sequence_, soft - rn - 1) &&
           util::qual_at(query_quality_, soft - rn - 1) > Config::kLowQual) {
        rn++;
    }
    if (rn > 0) {
        mch += rn;
        soft -= rn;
        if (soft > 0) {
            cigar_str_ = matcher.replace(i2s(soft) + "S" + i2s(mch) + "M");
        } else {
            cigar_str_ = matcher.replace(i2s(mch) + "M");
        }
        position_ -= rn;
        rn = 0;
    }
    if (soft > 0) {
        while (rn + 1 < soft && is_has_and_equals(reference_, position_ - rn - 2, query_sequence_, soft - rn - 2) &&
               util::qual_at(query_quality_, soft - rn - 2) > Config::kLowQual) {
            rn++;
            add_ref(rn_set, position_ - rn - 2);
        }
        const int rn_nt = static_cast<int>(rn_set.size());  // Don't adjust for homopolymers
        if ((rn > 4 && rn_nt > 1) || is_has_and_equals(reference_, position_ - 1, query_sequence_, soft - 1)) {
            mch += rn + 1;
            soft -= rn + 1;
            if (soft > 0) {
                cigar_str_ = replace_first(DIG_S_DIG_M, cigar_str_, i2s(soft) + "S" + i2s(mch) + "M");
            } else {
                cigar_str_ = replace_first(DIG_S_DIG_M, cigar_str_, i2s(mch) + "M");
            }
            position_ -= rn + 1;
        }
        if (rn == 0) {
            int rrn = 0;
            int rmch = 0;
            while (rrn < mch && rn < mch) {
                if (!reference_.contains(position_ + rrn)) {
                    break;
                }
                if (is_has_and_not_equals(reference_, position_ + rrn, query_sequence_, soft + rrn)) {
                    rn = rrn + 1;
                    rmch = 0;
                } else if (is_has_and_equals(reference_, position_ + rrn, query_sequence_, soft + rrn)) {
                    rmch++;
                }
                rrn++;
                if (rmch >= 3) {
                    break;
                }
            }
            if (rn > 0 && rn < mch) {
                soft += rn;
                mch -= rn;
                cigar_str_ = replace_first(DIG_S_DIG_M, cigar_str_, i2s(soft) + "S" + i2s(mch) + "M");
                position_ += rn;
            }
        }
    }
}

void Modifier::capture_mis_softly_3_mismatches(const Found& matcher) {
    const std::string ov5 = matcher.str(1);
    int mch = matcher.num(2);
    int refoff = position_ + matcher.num(2);
    int rdoff = matcher.num(2);
    if (!ov5.empty()) {
        refoff += sum_global_find(ALIGNED_LENGTH_MND, ov5);  // reference position
        rdoff += sum_global_find(SOFT_CLIPPED, ov5);         // read position
    }
    int rn = 0;
    int rrn = 0;
    int rmch = 0;
    while (rrn < mch && rn < mch) {
        if (!reference_.contains(refoff - rrn - 1)) {
            break;
        }
        if (rrn < rdoff && is_has_and_not_equals(reference_, refoff - rrn - 1, query_sequence_, rdoff - rrn - 1)) {
            rn = rrn + 1;
            rmch = 0;
        } else if (rrn < rdoff && is_has_and_equals(reference_, refoff - rrn - 1, query_sequence_, rdoff - rrn - 1)) {
            rmch++;
        }
        rrn++;
        // Stop at three consecutive matches
        if (rmch >= 3) {
            break;
        }
    }
    mch -= rn;
    if (rn > 0 && rn <= 3) {
        cigar_str_ = replace_first(DIG_M_END, cigar_str_, i2s(mch) + "M" + i2s(rn) + "S");
    }
}

void Modifier::capture_mis_softly_ms(const Found& matcher) {
    const std::string ov5 = matcher.str(1);   // CIGAR before the last matched segment
    int mch = matcher.num(2);                 // length of matched sequence
    int soft = matcher.num(3);                // length of soft-clipping
    int refoff = position_ + matcher.num(2);  // reference offset of the soft clip
    int rdoff = matcher.num(2);               // read offset of the soft clip
    if (!ov5.empty()) {
        refoff += sum_global_find(ALIGNED_LENGTH_MND, ov5);
        rdoff += sum_global_find(SOFT_CLIPPED, ov5);
    }
    // number of bases after refoff/rdoff that match in reference and read
    int rn = 0;
    std::set<int> rn_set;
    while (rn < soft && is_has_and_equals(reference_, refoff + rn, query_sequence_, rdoff + rn) &&
           util::qual_at(query_quality_, rdoff + rn) > Config::kLowQual) {
        rn++;
    }
    if (rn > 0) {
        mch += rn;
        soft -= rn;
        if (soft > 0) {
            cigar_str_ = replace_first(DIG_M_DIG_S_END, cigar_str_, i2s(mch) + "M" + i2s(soft) + "S");
        } else {
            cigar_str_ = replace_first(DIG_M_DIG_S_END, cigar_str_, i2s(mch) + "M");
        }
        rn = 0;
    }

    if (soft > 0) {
        while (rn + 1 < soft && is_has_and_equals(reference_, refoff + rn + 1, query_sequence_, rdoff + rn + 1) &&
               util::qual_at(query_quality_, rdoff + rn + 1) > Config::kLowQual) {
            rn++;
            add_ref(rn_set, refoff + rn + 1);
        }
        const int rn_nt = static_cast<int>(rn_set.size());  // don't adjust if homopolymer
        if (rn > 4 && rn_nt > 1) {
            mch += rn + 1;
            soft -= rn + 1;
            if (soft > 0) {
                cigar_str_ = replace_first(DIG_M_DIG_S_END, cigar_str_, i2s(mch) + "M" + i2s(soft) + "S");
            } else {
                cigar_str_ = replace_first(DIG_M_DIG_S_END, cigar_str_, i2s(mch) + "M");
            }
        }
        if (rn == 0) {
            int rrn = 0;
            int rmch = 0;
            while (rrn < mch && rn < mch) {
                if (!reference_.contains(refoff - rrn - 1)) {
                    break;
                }
                if (rrn < rdoff &&
                    is_has_and_not_equals(reference_, refoff - rrn - 1, query_sequence_, rdoff - rrn - 1)) {
                    rn = rrn + 1;
                    rmch = 0;
                } else if (rrn < rdoff &&
                           is_has_and_equals(reference_, refoff - rrn - 1, query_sequence_, rdoff - rrn - 1)) {
                    rmch++;
                }
                rrn++;
                if (rmch >= 3) {
                    break;
                }
            }
            if (rn > 0 && rn < mch) {
                soft += rn;
                mch -= rn;
                cigar_str_ = replace_first(DIG_M_DIG_S_END, cigar_str_, i2s(mch) + "M" + i2s(soft) + "S");
            }
        }
    }
}

bool Modifier::combine_to_close_to_one(const Found& matcher, bool flag) {
    const std::string op = matcher.str(5);
    const int g2 = matcher.num(2);
    const int g3 = matcher.num(3);
    const int g4 = matcher.num(4);
    if (g3 <= 15) {
        int dlen = g3;       // length of matched sequence and deletion
        int ilen = g2 + g3;  // length of first insertion and matched sequence
        if (op == "I") {
            ilen += g4;
        } else if (op == "D") {
            dlen += g4;
            if (auto istr = matcher.group(6)) {  // insertion after the deletion
                ilen += util::parse_int(std::string_view(*istr).substr(0, istr->size() - 1)).value();
            }
        }
        cigar_str_ = replace_first(DIG_I_DIG_M_DIG_DI_DIGI, cigar_str_, i2s(dlen) + "D" + i2s(ilen) + "I");
        flag = true;
    }
    return flag;
}

bool Modifier::combine_to_close_to_correct(const Found& matcher, bool flag) {
    const int g1 = matcher.num(1);
    const int g2 = matcher.num(2);
    const int g3 = matcher.num(3);
    if (g2 <= 15) {
        const std::string op = matcher.str(4);
        int dlen = g1 + g2;  // both deletions and the matched sequence
        int ilen = g2;       // matched sequence length
        if (op == "I") {
            ilen += g3;
        } else if (op == "D") {
            dlen += g3;
            if (auto istr = matcher.group(5)) {  // insertion after the second deletion
                ilen += util::parse_int(std::string_view(*istr).substr(0, istr->size() - 1)).value();
            }
        }
        cigar_str_ = matcher.replace(i2s(dlen) + "D" + i2s(ilen) + "I");
        flag = true;
    }
    return flag;
}

bool Modifier::three_indels(const Found& matcher, bool flag) {
    int tslen = matcher.num(5) + matcher.num(8);
    if (matcher.is(4, "I")) {
        tslen += matcher.num(3);
    }
    if (matcher.is(7, "I")) {
        tslen += matcher.num(6);
    }
    if (matcher.is(10, "I")) {
        tslen += matcher.num(9);
    }

    int dlen = matcher.num(5) + matcher.num(8);
    if (matcher.is(4, "D")) {
        dlen += matcher.num(3);
    }
    if (matcher.is(7, "D")) {
        dlen += matcher.num(6);
    }
    if (matcher.is(10, "D")) {
        dlen += matcher.num(9);
    }

    const int mid = matcher.num(5) + matcher.num(8);
    const std::string ov5 = matcher.str(1);
    int refoff = position_ + matcher.num(2);
    int rdoff = matcher.num(2);
    int RDOFF = matcher.num(2);
    int rm = matcher.num(11);

    if (!ov5.empty()) {  // the complex is not at the start of the CIGAR
        refoff += sum_global_find(ALIGNED_LENGTH_MND, ov5);
        rdoff += sum_global_find(SOFT_CLIPPED, ov5);
    }

    int rn = 0;
    while (rdoff + rn < static_cast<int>(query_sequence_.size()) &&
           is_has_and_equals(util::char_at(query_sequence_, rdoff + rn), reference_, refoff + rn)) {
        rn++;
    }
    RDOFF += rn;
    dlen -= rn;
    tslen -= rn;

    std::string new_cigar_str = i2s(RDOFF) + "M";
    if (tslen <= 0) {
        dlen -= tslen;
        rm += tslen;
        if (dlen == 0) {
            RDOFF = RDOFF + rm;
            new_cigar_str = i2s(RDOFF) + "M";
        } else if (dlen < 0) {
            tslen = -dlen;
            rm += dlen;
            if (rm < 0) {
                RDOFF = RDOFF + rm;
                new_cigar_str = i2s(RDOFF) + "M" + i2s(tslen) + "I";
            } else {
                new_cigar_str += i2s(tslen) + "I" + i2s(rm) + "M";
            }
        } else {
            new_cigar_str += i2s(dlen) + "D" + i2s(rm) + "M";
        }
    } else {
        if (dlen == 0) {
            new_cigar_str += i2s(tslen) + "I" + i2s(rm) + "M";
        } else if (dlen < 0) {
            rm += dlen;
            new_cigar_str += i2s(tslen) + "I" + i2s(rm) + "M";
        } else {
            new_cigar_str += i2s(dlen) + "D" + i2s(tslen) + "I" + i2s(rm) + "M";
        }
    }
    if (mid <= 15) {
        cigar_str_ = replace_first(DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM, cigar_str_, new_cigar_str);
        flag = true;
    }
    return flag;
}

bool Modifier::three_deletions(const Found& matcher, bool flag) {
    int tslen = matcher.num(4) + matcher.num(6);
    int dlen = matcher.num(3) + matcher.num(4) + matcher.num(5) + matcher.num(6) + matcher.num(7);
    const int mid = matcher.num(4) + matcher.num(6);
    const std::string ov5 = matcher.str(1);
    int refoff = position_ + matcher.num(2);
    int rdoff = matcher.num(2);
    int RDOFF = matcher.num(2);
    int rm = matcher.num(8);

    if (!ov5.empty()) {
        refoff += sum_global_find(ALIGNED_LENGTH_MND, ov5);
        rdoff += sum_global_find(SOFT_CLIPPED, ov5);
    }

    int rn = 0;
    while (rdoff + rn < static_cast<int>(query_sequence_.size()) &&
           is_has_and_equals(util::char_at(query_sequence_, rdoff + rn), reference_, refoff + rn)) {
        rn++;
    }
    RDOFF += rn;
    dlen -= rn;
    tslen -= rn;
    std::string new_cigar_str = i2s(RDOFF) + "M";
    if (tslen <= 0) {
        dlen -= tslen;
        rm += tslen;
        new_cigar_str += i2s(dlen) + "D" + i2s(rm) + "M";
    } else {
        new_cigar_str += i2s(dlen) + "D" + i2s(tslen) + "I" + i2s(rm) + "M";
    }
    if (mid <= 15) {
        cigar_str_ = replace_first(DM_DD_DM_DD_DM_DD_DM, cigar_str_, new_cigar_str);
        flag = true;
    }
    return flag;
}

bool Modifier::two_deletions_insertion_to_complex(const Found& matcher, bool flag) {
    int tslen = matcher.num(4) + matcher.num(5) + matcher.num(6);
    int dlen = matcher.num(3) + matcher.num(4) + matcher.num(6) + matcher.num(7);
    const int mid = matcher.num(4) + matcher.num(6);
    const std::string ov5 = matcher.str(1);
    int refoff = position_ + matcher.num(2);
    int rdoff = matcher.num(2);
    int RDOFF = matcher.num(2);
    int rm = matcher.num(8);

    if (!ov5.empty()) {
        refoff += sum_global_find(ALIGNED_LENGTH_MND, ov5);
        rdoff += sum_global_find(SOFT_CLIPPED, ov5);
    }

    int rn = 0;
    while (rdoff + rn < static_cast<int>(query_sequence_.size()) &&
           is_has_and_equals(util::char_at(query_sequence_, rdoff + rn), reference_, refoff + rn)) {
        rn++;
    }
    RDOFF += rn;
    dlen -= rn;
    tslen -= rn;

    std::string new_cigar_str = i2s(RDOFF) + "M";
    if (tslen <= 0) {
        dlen -= tslen;
        rm += tslen;
        new_cigar_str += i2s(dlen) + "D" + i2s(rm) + "M";
    } else {
        new_cigar_str += i2s(dlen) + "D" + i2s(tslen) + "I" + i2s(rm) + "M";
    }
    if (mid <= 15) {
        cigar_str_ = replace_first(D_M_D_DD_M_D_I_D_M_D_DD_prim, cigar_str_, new_cigar_str);
        flag = true;
    }
    return flag;
}

bool Modifier::begin_digit_m_number_i_or_d_number_m(const Found& matcher) {
    const int tmid = matcher.num(1);
    int mlen = matcher.num(4);
    int tn = 0;

    int tslen = tmid + (matcher.is(3, "I") ? matcher.num(2) : 0);
    position_ += tmid + (matcher.is(3, "D") ? matcher.num(2) : 0);
    while (tn < mlen && is_has_and_not_equals(util::char_at(query_sequence_, tslen + tn), reference_, position_ + tn)) {
        tn++;
    }
    tslen += tn;
    mlen -= tn;
    position_ += tn;
    cigar_str_ = replace_first(BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_, cigar_str_, i2s(tslen) + "S" + i2s(mlen) + "M");
    return true;
}

}  // namespace

ModifiedCigar modify_cigar(int position, const std::string& cigar, const std::string& query_sequence,
                           const std::string& query_quality, const io::Reference& reference,
                           int insertion_deletion_length, const Region& region, int max_read_length, const Config& conf,
                           ExceptionCounter& exceptions) {
    Modifier m(position, cigar, query_sequence, query_quality, reference, insertion_deletion_length, max_read_length,
               conf);
    try {
        m.run();
    } catch (const util::RecoverableError& e) {
        print_exception_and_continue(e, "cigar", std::to_string(m.position_) + " " + cigar, &region, exceptions);
    }
    return {m.position_, m.cigar_str_, m.query_sequence_, m.query_quality_};
}

}  // namespace senri::core
