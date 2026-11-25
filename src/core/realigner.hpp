#pragma once

//   VariationRealigner.process / adjustMNP / realignIndels (realigndel, realignins,
//   realignlgdel, realignlgins30, realignlgins)
// The SV search itself (StructuralVariantsProcessor) is in core/sv.cpp.

#include <string>
#include <utility>
#include <vector>

#include "core/scope.hpp"
#include "core/variation.hpp"
#include "core/variation_utils.hpp"

namespace omnivar::core {

// VariationRealigner.process without filterAllSVStructures.
void realigner_process(VariationData& data, RegionScope& scope);

// StructuralVariantsProcessor.process; core/sv.cpp.
void sv_process(VariationData& data, RegionScope& scope);

void adjust_mnp(VariationData& data, const RegionScope& scope);
void adj_snv(VariationData& data, const RegionScope& scope);

void realign_indels(VariationData& data, RegionScope& scope);

void filter_all_sv_structures(VariationData& data, const RegionScope& scope);

// new VariationRealigner().initFromScope(previousScope).realigndel(bams, dels): used by the SV
// processor (findINVsub).
void realigndel_with_bams(VariationData& data, RegionScope& scope, const CountMap& dels);

// ReferenceResource.getReference(region, extension, reference); a failure is a RecoverableError
// (an exception in the caller's loop).
void extend_reference(RegionScope& scope, const Region& region, int extension);

// AbstractMode.partialPipeline on `region`, parsing into data's maps.
void partial_pipeline(VariationData& data, RegionScope& scope, const Region& region, int max_read_length);

struct Match {
    int base_position = 0;
    std::string matched_sequence;
};

[[nodiscard]] Match find_match(std::string seq, const io::Reference& ref, int position, int dir, int seed, int mm,
                               const Config& conf);

[[nodiscard]] bool is_overlap(int start1, int end1, int start2, int end2, int rlen);

// StructuralVariantsProcessor.markSV: {cov, cnt, pairs}; markDUPSV: {cnt, pairs}.
struct SvMark {
    int cov = 0;
    int count = 0;
    int pairs = 0;
};
SvMark mark_sv(int start, int end, const std::vector<std::vector<Sclip>*>& svs, int rlen, const Config& conf);
SvMark mark_dup_sv(int start, int end, const std::vector<std::vector<Sclip>*>& svs, int rlen, const Config& conf);

SvInfo& get_sv(PositionMap& hash, int start);
void remove_sv(PositionMap& hash, int start);

// VariationUtils.joinRef(ref, from, double to): from <= i < to.
[[nodiscard]] std::string join_ref_double(const RefMap& ref, int from, double to);
[[nodiscard]] std::string join_ref_for_5_lgins(const RefMap& ref, int from, int to, const std::string& seq,
                                               const std::string& extra);
[[nodiscard]] std::string join_ref_for_3_lgins(const RefMap& ref, int from, int to, int shift5, const std::string& seq,
                                               const std::string& extra);

[[nodiscard]] bool ismatch(const std::string& seq1, const std::string& seq2, int dir, const Config& conf, int mm = 3);

// VariationRealigner.adjInsPos: {bi, ins}
[[nodiscard]] std::pair<int, std::string> adj_ins_pos(int bi, std::string ins, const RefMap& ref);

// Match35 and VariationRealigner.find35match
struct Match35 {
    int matched5end = 0;
    int matched3end = 0;
    int max_matched_length = 0;
};
[[nodiscard]] Match35 find35match(const std::string& seq5, const std::string& seq3, const Config& conf);

}  // namespace omnivar::core
