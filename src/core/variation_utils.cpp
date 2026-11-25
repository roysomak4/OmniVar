#include "core/variation_utils.hpp"

#include <cstdio>
#include <print>

#include "core/patterns.hpp"
#include "util/checked_string.hpp"
#include "util/sequence.hpp"

namespace omnivar::core {

char ref_at(const RefMap& ref, int pos) {
    auto r = ref_get(ref, pos);
    if (!r) {
        util::throw_npe("reference base missing at " + std::to_string(pos));
    }
    return *r;
}

bool is_has_and_equals(const RefMap& ref, int index1, std::string_view str, long index2) {
    auto r = ref_get(ref, index1);
    if (!r) {
        return false;
    }
    return *r == util::char_at(str, index2);
}

bool is_has_and_not_equals(const RefMap& ref, int index1, std::string_view str, long index2) {
    auto r = ref_get(ref, index1);
    if (!r) {
        return false;
    }
    return *r != util::char_at(str, index2);
}

int strand_bias(int forward_count, int reverse_count, const Config& conf) {
    // using p=0.01, because prop.test(1,12) = 0.01
    if (forward_count + reverse_count <= 12) {
        return forward_count * reverse_count > 0 ? 2 : 0;
    }
    const double total = forward_count + reverse_count;
    return (forward_count / total >= conf.bias && reverse_count / total >= conf.bias &&
            forward_count >= conf.min_bias_reads && reverse_count >= conf.min_bias_reads)
               ? 2
               : 1;
}

std::string findconseq(Sclip& soft_clip, int dir, const RunContext& run, const Config& conf) {
    if (soft_clip.sequence) {
        return *soft_clip.sequence;
    }
    int total = 0;
    int match = 0;
    std::string seq;
    bool flag = false;
    for (const auto& [position_in_sclip, nv] : soft_clip.nt) {
        int max_count = 0;
        double max_quality = 0;
        std::optional<char> chosen_base;
        int total_count = 0;
        for (const auto& [current_base, current_count] : nv) {
            total_count += current_count;
            bool better = current_count > max_count;
            if (!better) {
                auto sit = soft_clip.seq.find(position_in_sclip);
                if (sit != soft_clip.seq.end()) {
                    if (const Variation* v = sit->second.get(current_base);
                        v != nullptr && v->mean_quality > max_quality) {
                        better = true;
                    }
                }
            }
            if (better) {
                max_count = current_count;
                chosen_base = current_base;
                auto sit = soft_clip.seq.find(position_in_sclip);
                if (sit == soft_clip.seq.end() || sit->second.get(current_base) == nullptr) {
                    util::throw_npe("Sclip.seq missing in findconseq");
                }
                max_quality = sit->second.get(current_base)->mean_quality;
            }
        }
        if (position_in_sclip == 3 && soft_clip.nt.size() >= 6 &&
            total_count / static_cast<double>(soft_clip.vars_count) < 0.2 && total_count <= 2) {
            break;
        }
        if ((total_count - max_count > 2 || max_count <= total_count - max_count) &&
            max_count / static_cast<double>(total_count) < 0.8) {
            if (flag) {
                break;
            }
            flag = true;
        }
        total += total_count;
        match += max_count;
        if (chosen_base) {
            seq += *chosen_base;
        }
    }
    std::string result;
    const auto nt_size = static_cast<int>(soft_clip.nt.size());
    const auto seq_len = static_cast<int>(seq.size());
    if (total != 0 && match / static_cast<double>(total) > 0.9 && seq_len / 1.5 > nt_size - seq_len &&
        (seq_len / static_cast<double>(nt_size) > 0.8 || nt_size - seq_len < 10 || seq_len > 25)) {
        result = seq;
    }
    if (!result.empty() && result.size() > static_cast<std::size_t>(Config::kSeed2)) {
        if (patterns::begin_any_7(result, 'A') || patterns::begin_any_7(result, 'T')) {
            soft_clip.used = true;
        }
        if (islowcomplexseq(result)) {
            soft_clip.used = true;
        }
    }
    if (!result.empty() && result.size() >= static_cast<std::size_t>(Config::kAdSeed)) {
        if (dir == 3) {
            if (run.adaptor_forward.contains(util::substr(result, 0, Config::kAdSeed))) {
                result.clear();
            }
        } else if (dir == 5) {
            if (run.adaptor_reverse.contains(util::reverse(util::substr(result, 0, Config::kAdSeed)))) {
                result.clear();
            }
        }
    }
    soft_clip.sequence = result;
    if (conf.y) {
        std::println(stderr, "  Candidate consensus: {} Reads: {} M: {} T: {} Final: {}", seq, soft_clip.vars_count,
                     match, total, result);
    }
    return result;
}

std::string join_ref(const RefMap& ref, int from, int to) {
    std::string sb;
    for (int i = from; i <= to; ++i) {
        if (auto c = ref_get(ref, i)) {
            sb += *c;
        }
    }
    return sb;
}

void adj_cnt(Variation& to_add, const Variation& variant, Variation* reference_var) {
    to_add.vars_count += variant.vars_count;
    to_add.extracnt += variant.vars_count;
    to_add.high_quality_reads_count += variant.high_quality_reads_count;
    to_add.low_quality_reads_count += variant.low_quality_reads_count;
    to_add.mean_position += variant.mean_position;
    to_add.mean_quality += variant.mean_quality;
    to_add.mean_mapping_quality += variant.mean_mapping_quality;
    to_add.number_of_mismatches += variant.number_of_mismatches;
    to_add.pstd = true;
    to_add.qstd = true;
    to_add.add_dir(true, variant.get_dir(true));
    to_add.add_dir(false, variant.get_dir(false));
    if (reference_var == nullptr) {
        return;
    }
    reference_var->vars_count -= variant.vars_count;
    reference_var->high_quality_reads_count -= variant.high_quality_reads_count;
    reference_var->low_quality_reads_count -= variant.low_quality_reads_count;
    reference_var->mean_position -= variant.mean_position;
    reference_var->mean_quality -= variant.mean_quality;
    reference_var->mean_mapping_quality -= variant.mean_mapping_quality;
    reference_var->number_of_mismatches -= variant.number_of_mismatches;
    reference_var->sub_dir(true, variant.get_dir(true));
    reference_var->sub_dir(false, variant.get_dir(false));
    correct_cnt(*reference_var);
}

void correct_cnt(Variation& v) {
    if (v.vars_count < 0) {
        v.vars_count = 0;
    }
    if (v.high_quality_reads_count < 0) {
        v.high_quality_reads_count = 0;
    }
    if (v.low_quality_reads_count < 0) {
        v.low_quality_reads_count = 0;
    }
    if (v.mean_position < 0) {
        v.mean_position = 0;
    }
    if (v.mean_quality < 0) {
        v.mean_quality = 0;
    }
    if (v.mean_mapping_quality < 0) {
        v.mean_mapping_quality = 0;
    }
    if (v.get_dir(true) < 0) {
        v.add_dir(true, -v.get_dir(true));
    }
    if (v.get_dir(false) < 0) {
        v.add_dir(false, -v.get_dir(false));
    }
}

bool islowcomplexseq(std::string_view seq) {
    const auto len = static_cast<double>(seq.size());
    if (seq.empty()) {
        return true;
    }
    int ntcnt = 0;
    for (char base : {'A', 'T', 'G', 'C'}) {
        int n = 0;
        for (char c : seq) {
            n += c == base;
        }
        if (n > 0) {
            ++ntcnt;
        }
        if (n / len > 0.75) {
            return true;
        }
    }
    return ntcnt < 3;
}

bool ismatchref(std::string_view sequence, const RefMap& ref, int position, int dir, int mm_limit) {
    int mm = 0;
    for (int n = 0; n < static_cast<int>(sequence.size()); ++n) {
        auto ref_ch = ref_get(ref, position + dir * n);
        if (!ref_ch) {
            return false;
        }
        // Utils.charAt(sequence, dir == 1 ? n : dir * n - 1)
        if (util::utils_char_at(sequence, dir == 1 ? n : dir * n - 1) != static_cast<unsigned char>(*ref_ch)) {
            ++mm;
        }
    }
    return mm <= mm_limit && mm / static_cast<double>(sequence.size()) < 0.15;
}

}  // namespace omnivar::core
