// synthetic: writes the synthetic TSV rows for the VCF annotation checks to stdout.
//
// --single (tests/golden/vcf_synthetic.tsv): real Simple rows from the k0_U expectations, re-placed
// to reach var2vcf_valid.pl paths the single-region golden cases do not:
//   - contig order (reorder): numeric names, chrX/chrY/chrM, names with '_' or no digit;
//   - several alleles at one start, with equal AFs (stable sort) and duplicate keys (%seen);
//   - SV alleles (<DEL>, <INV>, <DUP>) with and without split-read support, SV_info "a-b-c";
//   - PASS SNVs a few bp apart (Cluster with -c); a sample name with whitespace.
// --paired (tests/golden/vcf_synthetic_paired.tsv): real Somatic rows (55 columns) from the golden
// Somatic cases: clustered PASS somatic SNVs, contig order, alleles at one start, every STATUS.
//
// Deterministic: the generator is the Mersenne Twister with the seeding and choice/sample
// algorithms of CPython's `random` module, so the checked-in files regenerate byte for byte.

#include <array>
#include <bit>
#include <cstdint>
#include <print>
#include <set>

#include "testkit.hpp"

namespace testkit {

namespace {

// CPython's random.Random(int seed): MT19937 seeded with init_by_array(32-bit chunks of the seed).
class PyRandom {
public:
    explicit PyRandom(std::uint32_t seed) {
        init_genrand(19650218U);
        const std::uint32_t key[] = {seed};
        std::size_t i = 1, j = 0;
        for (std::size_t k = std::max<std::size_t>(N, 1); k; --k) {
            mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1664525U)) + key[j] + static_cast<std::uint32_t>(j);
            ++i;
            ++j;
            if (i >= N) {
                mt_[0] = mt_[N - 1];
                i = 1;
            }
            if (j >= 1) j = 0;
        }
        for (std::size_t k = N - 1; k; --k) {
            mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1566083941U)) - static_cast<std::uint32_t>(i);
            ++i;
            if (i >= N) {
                mt_[0] = mt_[N - 1];
                i = 1;
            }
        }
        mt_[0] = 0x80000000U;
    }

    // random._randbelow_with_getrandbits
    std::size_t below(std::size_t n) {
        const int k = static_cast<int>(std::bit_width(n));
        auto r = getrandbits(k);
        while (r >= n) r = getrandbits(k);
        return r;
    }

    template <class T>
    const T& choice(const std::vector<T>& v) {
        return v[below(v.size())];
    }

    // random.sample(population, 2) for a population larger than 21 (the set-based branch)
    template <class T>
    std::array<T, 2> sample2(const std::vector<T>& v) {
        std::array<T, 2> out;
        std::set<std::size_t> selected;
        for (std::size_t i = 0; i < 2; ++i) {
            if (v.size() <= 21) die("sample2: population too small for the set-based branch");
            auto j = below(v.size());
            while (selected.contains(j)) j = below(v.size());
            selected.insert(j);
            out[i] = v[j];
        }
        return out;
    }

private:
    static constexpr std::size_t N = 624, M = 397;

    void init_genrand(std::uint32_t s) {
        mt_[0] = s;
        for (std::size_t i = 1; i < N; ++i)
            mt_[i] = 1812433253U * (mt_[i - 1] ^ (mt_[i - 1] >> 30)) + static_cast<std::uint32_t>(i);
        mti_ = N;
    }

    std::uint32_t genrand() {
        if (mti_ >= N) {
            for (std::size_t k = 0; k < N; ++k) {
                const std::uint32_t y = (mt_[k] & 0x80000000U) | (mt_[(k + 1) % N] & 0x7fffffffU);
                mt_[k] = mt_[(k + M) % N] ^ (y >> 1) ^ ((y & 1U) ? 0x9908b0dfU : 0U);
            }
            mti_ = 0;
        }
        std::uint32_t y = mt_[mti_++];
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680U;
        y ^= (y << 15) & 0xefc60000U;
        y ^= y >> 18;
        return y;
    }

    std::size_t getrandbits(int k) {
        if (k <= 0) return 0;
        if (k > 32) die("getrandbits: more than 32 bits");
        return genrand() >> (32 - k);
    }

    std::array<std::uint32_t, N> mt_{};
    std::size_t mti_ = N + 1;
};

using Row = std::vector<std::string>;

double num(const std::string& s) {
    const auto v = as_number(s);
    if (!v) die("synthetic: not a number: " + s);
    return *v;
}

Row place(Row r, const std::string& chrom, long start, const std::string& sample, std::size_t region_col) {
    const long length = std::stol(r[4]) - std::stol(r[3]);
    r[0] = sample;
    r[2] = chrom;
    r[3] = std::to_string(start);
    r[4] = std::to_string(start + length);
    r[region_col] = std::format("{}:{}-{}", chrom, start - 100, start + 100);
    return r;
}

std::vector<Row> single_rows() {
    std::vector<Row> rows;
    for (const auto& path : case_files(golden_dir() / "expected_k0_U", false)) {
        const auto n = path.filename().string();
        if (!n.starts_with("Simple") || n.contains("-J") || n.contains("fisher")) continue;
        for (const auto& line : lines(read_file(path))) {
            auto f = split(line, '\t');
            if (f.size() == 36 && !f[5].empty()) rows.push_back(std::move(f));
        }
    }
    std::vector<Row> typed, refrows;
    for (const auto& r : rows) (r[33].empty() ? refrows : typed).push_back(r);
    PyRandom rnd(7);
    const std::vector<std::string> contigs = {"chr10",      "chr2",  "chrX",  "chrM", "chrY", "chr1_gl000191_random",
                                              "GL000192.1", "HLA-A", "chrUn", "1",    "12",   "chr21",
                                              "M"};
    // Not both MT and chrM (or X and chrX, Y and chrY): reorder() would drop the second name, which
    // the omnivar VCF deliberately keeps (FILTER annotation never removes records).
    const std::vector<long> steps = {1, 2, 3, 5, 8, 15, 40};
    const std::string sample = "synthetic sample 1";
    std::vector<Row> out;
    long pos = 1000;
    for (const auto& chrom : contigs) {
        for (int i = 0; i < 12; ++i) {
            pos += rnd.choice(steps);
            out.push_back(place(rnd.choice(typed), chrom, pos, sample, 32));
        }
        // several alleles at one start, including an exact duplicate and an AF tie
        const auto [a, b] = rnd.sample2(typed);
        pos += 50;
        auto a2 = a;
        a2[14] = b[14];
        out.push_back(place(a, chrom, pos, sample, 32));
        out.push_back(place(b, chrom, pos, sample, 32));
        out.push_back(place(a2, chrom, pos, sample, 32));
        out.push_back(place(a, chrom, pos, sample, 32));
        out.push_back(place(rnd.choice(refrows), chrom, pos + 1, sample, 32));
    }
    // PASS SNVs close together (Cluster), on chr2
    std::vector<Row> passing;
    for (const auto& r : typed) {
        if (r[33] == "SNV" && num(r[14]) >= 0.3 && num(r[16]) >= 8 && num(r[18]) >= 22.5 && num(r[20]) >= 10 &&
            num(r[21]) >= 1.5 && std::stol(r[28]) >= 2)
            passing.push_back(r);
    }
    for (const long d : {0, 3, 4, 9, 30, 31, 45, 52})
        out.push_back(place(rnd.choice(passing), "chr2", 500000 + d, sample, 32));
    // SV alleles: ref base + symbolic ALT; SV_info "splitreads-spanpairs-cluster"
    long spos = 700000;
    const std::tuple<const char*, const char*, long> svs[] = {{"<DEL>", "2-3-1", 1500},
                                                              {"<INV>", "1-0-0", 2400},
                                                              {"<DUP>", "0-5-0", 900},
                                                              {"<DEL>", "0-0-0", 3000},
                                                              {"<INV>", "3-4-2", 12000}};
    for (const auto& [alt, info, length] : svs) {
        auto r = place(rnd.choice(typed), "chr21", spos, sample, 32);
        r[5] = r[5].substr(0, 1);
        r[6] = alt;
        r[4] = std::to_string(spos + length);
        const std::string a(alt);
        r[33] = a.substr(1, a.size() - 2);
        r[35] = info;
        out.push_back(std::move(r));
        spos += 20000;
    }
    return out;
}

// passes every tumour check of var2vcf_paired.pl with the default thresholds
bool good_tumour(const Row& r) {
    return r[50] == "SNV" && r[49] == "StrongSomatic" && std::stol(r[7]) >= 5 && std::stol(r[8]) >= 3 &&
           num(r[14]) >= 0.02 && num(r[16]) >= 8 && r[17] == "1" && num(r[18]) >= 22.5 && num(r[20]) >= 10 &&
           num(r[21]) >= 1.5 && num(r[24]) < 5.25 && num(r[44]) <= 12;
}

std::vector<Row> paired_rows() {
    std::vector<Row> rows;
    for (const auto& path : case_files(golden_cases(), false)) {
        const auto n = path.filename().string();
        if (!n.starts_with("Somatic") || n.contains("fisher")) continue;
        const auto ls = lines(read_file(path));
        for (std::size_t i = 1; i < ls.size(); ++i) {
            auto f = split(ls[i], '\t');
            if (f.size() == 55 && !f[5].empty()) rows.push_back(std::move(f));
        }
    }
    PyRandom rnd(11);
    const std::vector<std::string> contigs = {
        "chr3", "chr11", "chrX", "chrY", "chrM", "chr6_GL000250v2_alt", "chrUn_KI270742v1", "EBV", "chr1"};
    const std::vector<long> steps = {1, 2, 3, 5, 8, 15, 40};
    const std::string sample = "synthetic tumour";
    std::vector<Row> out;
    long pos = 1000;
    for (const auto& chrom : contigs) {
        for (int i = 0; i < 10; ++i) {
            pos += rnd.choice(steps);
            out.push_back(place(rnd.choice(rows), chrom, pos, sample, 48));
        }
        const auto [a, b] = rnd.sample2(rows);
        pos += 50;
        auto a2 = a;
        a2[14] = b[14];
        out.push_back(place(a, chrom, pos, sample, 48));
        out.push_back(place(b, chrom, pos, sample, 48));
        out.push_back(place(a2, chrom, pos, sample, 48));
        out.push_back(place(a, chrom, pos, sample, 48));
    }
    std::vector<Row> passing;
    for (const auto& r : rows) {
        if (good_tumour(r)) passing.push_back(r);
    }
    if (passing.empty()) die("no PASS somatic SNV rows in the golden cases");
    for (const long d : {0, 3, 4, 9, 30, 31, 45, 52, 53})
        out.push_back(place(rnd.choice(passing), "chr11", 400000 + d, sample, 48));
    return out;
}

}  // namespace

int cmd_synthetic(Args& args) {
    const bool single = args.flag("--single"), paired = args.flag("--paired");
    args.finish();
    if (single == paired) die("synthetic: give exactly one of --single or --paired");
    for (const auto& r : single ? single_rows() : paired_rows()) std::println("{}", join(r, "\t"));
    return 0;
}

}  // namespace testkit
