#pragma once

// Per-allele accumulators and per-region maps:
//   VariationMap   (LinkedHashMap<String, Variation> + SV info)
//   VariationData  (+ RealignedVariationData fields)
// Upstream HashMap fields use util::HashOrderMap so that iteration order matches.

#include <map>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include "util/ordered_maps.hpp"

namespace omnivar::core {

struct Variation {
    int vars_count = 0;             // varsCount
    int vars_count_on_forward = 0;  // varsCountOnForward
    int vars_count_on_reverse = 0;  // varsCountOnReverse
    double mean_position = 0;       // sum of read positions until ToVarsBuilder divides
    double mean_quality = 0;
    double mean_mapping_quality = 0;
    double number_of_mismatches = 0;
    int low_quality_reads_count = 0;
    int high_quality_reads_count = 0;
    bool pstd = false;
    bool qstd = false;
    int pp = 0;
    double pq = 0;
    int extracnt = 0;

    void inc_dir(bool dir) { dir ? ++vars_count_on_reverse : ++vars_count_on_forward; }
    void dec_dir(bool dir) { dir ? --vars_count_on_reverse : --vars_count_on_forward; }
    [[nodiscard]] int get_dir(bool dir) const { return dir ? vars_count_on_reverse : vars_count_on_forward; }
    void add_dir(bool dir, int add) { (dir ? vars_count_on_reverse : vars_count_on_forward) += add; }
    void sub_dir(bool dir, int sub) { (dir ? vars_count_on_reverse : vars_count_on_forward) -= sub; }
};

struct Mate {
    int mate_start_ms = 0;
    int mate_end_me = 0;
    int mate_length_mlen = 0;
    int start_s = 0;
    int end_e = 0;
    double pmean_rp = 0;
    double qmean_q = 0;
    double qmean_big_q = 0;
    double nm = 0;
};

struct Sclip : Variation {
    std::map<int, std::map<char, int>> nt;                   // TreeMap<Integer, TreeMap<Character, Integer>>
    std::map<int, util::HashOrderMap<char, Variation>> seq;  // TreeMap<Integer, HashMap<Character, Variation>>
    std::optional<std::string> sequence;                     // null until findconseq
    bool used = false;
    int start = 0;
    int end = 0;
    int mstart = 0;
    int mend = 0;
    int mlen = 0;
    int disc = 0;
    int softp = 0;
    util::LinkedMap<int, int> soft;
    std::vector<Mate> mates;
};

// VariationMap.SV
struct SvInfo {
    std::string type;
    int pairs = 0;
    int splits = 0;
    int clusters = 0;
};

// VariationMap<String, Variation>: alleles at one position, in insertion order.
struct VariationMap : util::LinkedMap<std::string, Variation> {
    std::optional<SvInfo> sv;
};

using PositionMap = util::HashOrderMap<int, VariationMap>;
using CountMap = util::HashOrderMap<int, util::HashOrderMap<std::string, int>>;

struct SVStructures {
    int svdelfend = 0;
    int svdelrend = 0;
    std::vector<Sclip> svfdel, svrdel;
    int svdupfend = 0;
    int svduprend = 0;
    std::vector<Sclip> svfdup, svrdup;
    int svinvfend3 = 0;
    int svinvrend3 = 0;
    int svinvfend5 = 0;
    int svinvrend5 = 0;
    std::vector<Sclip> svfinv3, svfinv5, svrinv3, svrinv5;
    util::HashOrderMap<std::string, std::vector<Sclip>> svffus, svrfus;
    util::HashOrderMap<std::string, int> svfusfend, svfusrend;
};

// Output of CigarParser (VariationData), carried through realignment and SV processing.
struct VariationData {
    PositionMap non_insertion_variants;
    PositionMap insertion_variants;
    CountMap position_to_insertion_count;
    CountMap position_to_deletion_count;
    CountMap mnp;
    util::HashOrderMap<int, int> ref_coverage;
    util::HashOrderMap<int, Sclip> soft_clips_5end;
    util::HashOrderMap<int, Sclip> soft_clips_3end;
    util::HashOrderMap<std::string, int> splice_count;
    SVStructures sv_structures;
    // VariationRealigner.SOFTP2SV: soft-clip position -> the SV clusters pointing at it (into
    // sv_structures, which no stage resizes after CigarParser), by count, highest first.
    util::HashOrderMap<int, std::vector<Sclip*>> softp2sv;
    int max_read_length = 0;
    double duprate = 0;
};

// VariationUtils.incCnt for Map<K, Integer>
template <class Map, class K>
void inc_cnt(Map& counts, const K& key, int add) {
    if (int* v = counts.get(key)) {
        *v += add;
    } else {
        counts.put(key, add);
    }
}

// VariationUtils.incCnt for a TreeMap<Character, Integer>
inline void inc_cnt_tree(std::map<char, int>& counts, char key) { ++counts[key]; }

// VariationUtils.getVariation: creates the position map and the allele entry if absent.
inline Variation& get_variation(PositionMap& hash, int start, const std::string& description) {
    return hash.get_or_create(start).get_or_create(description);
}

// VariationUtils.getVariationMaybe(hash, start, refBase): null if the base is null or absent.
inline Variation* get_variation_maybe(PositionMap& hash, int start, std::optional<char> ref_base) {
    if (!ref_base) {
        return nullptr;
    }
    VariationMap* map = hash.get(start);
    if (map == nullptr) {
        return nullptr;
    }
    return map->get(std::string(1, *ref_base));
}

// VariationUtils.getVariationFromSeq
inline Variation& get_variation_from_seq(Sclip& sclip, int idx, char ch) { return sclip.seq[idx].get_or_create(ch); }

}  // namespace omnivar::core
