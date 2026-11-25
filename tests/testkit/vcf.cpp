// vcf: the omnivar VCF against the reference chain on the 1.8.3 jar's output,
//
//     VarDict-1.8.3.jar | Rscript teststrandbias.R | perl var2vcf_valid.pl -A      (Simple)
//     VarDict-1.8.3.jar | Rscript testsomatic.R    | perl var2vcf_paired.pl -A     (Somatic)
//
// Expected VCFs live in tests/golden/expected_vcf_<stage>/ and are made with --generate from the
// jar's TSV for that stage (tests/golden/expected_<stage>/, or docs/golden for "full"). A case with
// --fisher is re-run through the jar without it (the chain reads the default 36-column layout).
// Skipped: -i (splicing output) and -J (37 columns, which teststrandbias.R rejects).
// Comparison: vcfcmp.hpp. Every omnivar VCF must also parse with `bcftools view` cleanly.
//
// --synthetic also checks tests/golden/vcf_synthetic.tsv (and vcf_synthetic_paired.tsv) through
// `omnivar --dev-tsv-to-vcf` against the chain run live, for every option variant.
// -p cases (M10): the expectation is the chain on the jar's output without -p; rejected records
// (FILTER REJ_*) are left out.

#include <algorithm>
#include <map>
#include <print>
#include <stdexcept>

#include "testkit.hpp"
#include "vcfcmp.hpp"

namespace testkit {

namespace {

using Variant = std::pair<std::string, std::string>;  // (perl options, omnivar options)

const std::map<std::string, Variant>& variants(bool paired) {
    static const std::string thresholds =
        "--filter-min-depth 10 --filter-min-var-depth 4 --filter-min-af 0.05 --filter-min-pmean 10 "
        "--filter-min-qual 30 --filter-min-mq 20 --filter-min-sn 3 --filter-max-nm 2.5 --filter-msi 6 "
        "--filter-hom-af 0.3 ";
    static const std::map<std::string, Variant> single = {
        {"default", {"-A", ""}},
        {"top", {"", "--vcf-top-allele-only"}},
        {"cluster", {"-A -c 10", "--filter-cluster-bp 10"}},
        {"thresholds",
         {"-A -d 10 -v 4 -f 0.05 -p 10 -q 30 -Q 20 -o 3 -m 2.5 -I 6 -F 0.3 -P 0", thresholds + "--no-filter-pstd"}},
        {"somatic-only", {"-A", ""}},  // no such option for the single-sample script
    };
    static const std::map<std::string, Variant> paired_v = {
        {"default", {"-A", ""}},
        {"top", {"", "--vcf-top-allele-only"}},
        {"cluster", {"-A -c 10", "--filter-cluster-bp 10"}},
        {"thresholds",
         {"-A -d 10 -v 4 -f 0.05 -p 10 -q 30 -Q 20 -o 3 -m 2.5 -I 6 -F 0.3 -P 0.01",
          thresholds + "--filter-max-pvalue 0.01"}},
        {"somatic-only", {"-A -M", "--vcf-somatic-only"}},
    };
    return paired ? paired_v : single;
}

struct StageDef {
    std::optional<fs::path> tsv_dir;
    std::string extra;
};

StageDef stage_def(const std::string& stage) {
    if (stage == "k0_U") return {golden_dir() / "expected_k0_U", "-k 0 -U"};
    if (stage == "U") return {golden_dir() / "expected_U", "-U"};
    if (stage == "full") return {std::nullopt, ""};
    die("vcf: --stage must be k0_U, U or full");
}

std::optional<std::string> skip_reason(const Case& c) {
    const auto a = words(c.args);
    if (contains(a, "-i") || contains(a, "--splice")) return "splicing output (-i)";
    if (contains(a, "-J")) return "CRISPR column: teststrandbias.R rejects 37 columns (-J)";
    return std::nullopt;
}

std::vector<std::string> case_args(const Case& c, const std::string& bed, const std::string& extra,
                                   const std::vector<std::string>& drop = {}) {
    std::vector<std::string> v = {bed, "-b", case_bams(c, golden_bams())};
    for (auto& x : words(c.args)) {
        if (!contains(drop, x)) v.push_back(std::move(x));
    }
    for (auto& x : words(extra)) v.push_back(std::move(x));
    for (auto& x : words("-z -c 1 -S 2 -E 3 -g 4 -G")) v.push_back(std::move(x));
    v.push_back(case_fasta_csv(c, golden_fastas()));
    return v;
}

std::string drop_rejected(std::string_view vcf) {
    std::string out;
    for (const auto& l : lines_keepends(vcf)) {
        if (l.starts_with('#') || !l.contains("\tREJ_")) out += l;
    }
    return out;
}

struct Opts {
    std::string binary, stage, variant;
    fs::path expected;
    int timeout = 600;
};

// The jar's TSV for the case at this stage, in the layout the chain reads (no --fisher, no -p).
std::string java_tsv(const Case& c, const std::string& stage, int timeout) {
    auto [tsv_dir, extra] = stage_def(stage);
    if (case_has_arg(c, "-p")) tsv_dir = golden_dir() / "expected_nopileup" / stage;
    if (!case_has_arg(c, "--fisher")) {
        if (!tsv_dir) return c.expected;
        return read_file(*tsv_dir / c.name);
    }
    const TempFile bed(".bed", case_bed_line(c, false));
    std::vector<std::string> cmd = {"java", "-cp", jar_classpath(), (root() / "tests/oracle/OracleMain.java").string(),
                                    case_fasta_csv(c, golden_fastas())};
    for (auto& x : case_args(c, bed.path(), extra, {"--fisher", "-p"})) cmd.push_back(std::move(x));
    const auto p = run(cmd, timeout_opt(timeout));
    if (p.rc || p.timed_out) throw std::runtime_error("jar failed: " + tail(p.err, 300));
    return p.out;
}

struct Result {
    std::string name, status;
    std::vector<std::string> detail;
};

std::size_t record_lines(std::string_view vcf) {
    std::size_t n = 0;
    for (const auto& l : lines(vcf)) n += !l.empty() && !l.starts_with('#');
    return n;
}

Result generate(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    if (c.mode != "Simple" && c.mode != "Somatic") return {c.name, "skip", {c.mode + " mode"}};
    if (auto r = skip_reason(c)) return {c.name, "skip", {*r}};
    try {
        const auto vcf = chain(java_tsv(c, a.stage, a.timeout), a.timeout, "-A", c.mode == "Somatic");
        write_file(a.expected / c.name, vcf);
        return {c.name, "ok", {std::format("{} records", record_lines(vcf))}};
    } catch (const std::exception& e) {
        return {c.name, "fail", {e.what()}};
    }
}

std::vector<std::string> known_lines(const std::vector<std::string>& known) {
    std::vector<std::string> v;
    for (const auto& k : known) v.push_back("known R artefact: " + k);
    return v;
}

Result compare(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    if (c.mode != "Simple" && c.mode != "Somatic") return {c.name, "skip", {c.mode + " mode"}};
    if (auto r = skip_reason(c)) return {c.name, "skip", {*r}};
    const auto& [perl_args, cpp_args] = variants(c.mode == "Somatic").at(a.variant);
    std::string expected;
    try {
        if (a.variant == "default") {
            const auto exp = a.expected / c.name;
            if (!fs::exists(exp))
                return {c.name, "fail", {"no expected VCF " + exp.string() + " (run with --generate)"}};
            expected = read_file(exp);
        } else {
            expected = chain(java_tsv(c, a.stage, a.timeout), a.timeout, perl_args, c.mode == "Somatic");
        }
    } catch (const std::exception& e) {
        return {c.name, "fail", {e.what()}};
    }
    const TempFile bed(".bed", case_bed_line(c, false));
    std::vector<std::string> cmd = {a.binary};
    for (auto& x : case_args(c, bed.path(), stage_def(a.stage).extra + " " + cpp_args)) cmd.push_back(std::move(x));
    const auto p = run(cmd, timeout_opt(a.timeout));
    if (p.rc || p.timed_out) return {c.name, "fail", {std::format("exit {}: {}", p.rc, tail(p.err, 300))}};
    std::vector<std::string> problems;
    if (auto bad = bcftools_check(p.out)) problems.push_back("bcftools: " + *bad);
    const auto want = parse_records(expected);
    const auto got = parse_records(case_has_arg(c, "-p") ? drop_rejected(p.out) : p.out);
    std::vector<std::string> known;
    for (auto& d : compare_records(want, got, known)) problems.push_back(std::move(d));
    if (!problems.empty()) return {c.name, "fail", problems};
    std::vector<std::string> detail = {std::format("{} records", want.count())};
    for (auto& k : known_lines(known)) detail.push_back(std::move(k));
    return {c.name, "pass", detail};
}

// vcf_synthetic[_paired].tsv through --dev-tsv-to-vcf against the live chain.
Result synthetic(const Opts& a, const std::string& variant, bool paired) {
    const auto& [perl_args, cpp_args] = variants(paired).at(variant);
    const std::string file = paired ? "vcf_synthetic_paired.tsv" : "vcf_synthetic.tsv";
    const auto name = file + " [" + variant + "]";
    const auto tsv = read_file(golden_dir() / file);
    std::string expected;
    try {
        expected = chain(tsv, a.timeout, perl_args, paired);
    } catch (const std::exception& e) {
        return {name, "fail", {e.what()}};
    }
    std::vector<std::string> cmd = {a.binary, "--dev-tsv-to-vcf", "-"};
    for (auto& x : words(cpp_args)) cmd.push_back(std::move(x));
    const auto p = run(cmd, input_opt(tsv, a.timeout));
    if (p.rc || p.timed_out) return {name, "fail", {std::format("exit {}: {}", p.rc, tail(p.err, 300))}};
    std::vector<std::string> problems;
    if (auto bad = bcftools_check(p.out)) problems.push_back("bcftools: " + *bad);
    const auto want = parse_records(expected), got = parse_records(p.out);
    std::vector<std::string> known;
    for (auto& d : compare_records(want, got, known)) problems.push_back(std::move(d));
    // Contig order: names with a digit and no '_' are sorted by their digits, then X, Y, MT (or
    // chrX, chrY, chrM); the rest come in Perl hash order, which changes from run to run, so only
    // their set is compared. Records within a contig must be in the same order.
    auto ordered = [](const std::string& ctg) {
        static const std::vector<std::string> special = {"X", "chrX", "Y", "chrY", "MT", "chrM"};
        return contains(special, ctg) ||
               (std::ranges::any_of(ctg, [](char ch) { return ch >= '0' && ch <= '9'; }) && !ctg.contains('_'));
    };
    auto contigs = [](const std::vector<RecordKey>& order) {
        std::vector<std::string> v;
        for (const auto& k : order) {
            if (!contains(v, std::get<0>(k))) v.push_back(std::get<0>(k));
        }
        return v;
    };
    const auto wc = contigs(want.order), gc = contigs(got.order);
    auto only_ordered = [&](const std::vector<std::string>& v) {
        std::vector<std::string> o;
        for (const auto& x : v) {
            if (ordered(x)) o.push_back(x);
        }
        return o;
    };
    auto as_set = [](std::vector<std::string> v) {
        std::ranges::sort(v);
        return v;
    };
    if (only_ordered(wc) != only_ordered(gc) || as_set(wc) != as_set(gc))
        problems.push_back("contig order: want [" + join(wc, ", ") + "] got [" + join(gc, ", ") + "]");
    auto per_contig = [](const std::vector<RecordKey>& order, const std::string& ctg) {
        std::vector<RecordKey> v;
        for (const auto& k : order) {
            if (std::get<0>(k) == ctg) v.push_back(k);
        }
        return v;
    };
    for (const auto& ctg : wc) {
        if (per_contig(want.order, ctg) != per_contig(got.order, ctg))
            problems.push_back("record order differs on " + ctg);
    }
    if (!problems.empty()) return {name, "fail", problems};
    return {name, "pass", {std::format("{} records", want.order.size())}};
}

}  // namespace

int cmd_vcf(Args& args) {
    Opts a;
    const auto binary = args.value("--binary");
    const bool gen = args.flag("--generate");
    a.stage = args.value_or("--stage", "k0_U");
    auto modes = args.values("--mode");
    a.variant = args.value_or("--variant", "default");
    const auto expected = args.value("--expected");
    const bool synth = args.flag("--synthetic");
    const auto match = args.value_or("-k", "");
    const int diff = args.int_or("--diff", 3);
    a.timeout = args.int_or("--timeout", 600);
    const int jobs = args.int_or("-j", default_jobs());
    args.finish();
    stage_def(a.stage);
    if (!variants(false).contains(a.variant)) die("vcf: unknown --variant " + a.variant);
    a.expected = expected ? fs::path(*expected) : golden_dir() / ("expected_vcf_" + a.stage);
    if (modes.empty()) modes = {"Simple", "Somatic"};
    if (!gen && !binary) die("vcf: --binary is required unless --generate");
    if (binary) a.binary = fs::absolute(*binary).string();
    if (gen) fs::create_directories(a.expected);

    std::vector<fs::path> paths;
    for (const auto& p : case_files(golden_cases(), false)) {
        const auto n = p.filename().string();
        if (std::ranges::any_of(modes, [&](const std::string& m) { return n.starts_with(m); }) && n.contains(match))
            paths.push_back(p);
    }
    std::vector<std::pair<std::string, bool>> synth_jobs;
    if (synth && !gen) {
        if (contains(modes, "Simple")) {
            for (const auto& [v, _] : variants(false)) {
                if (v != "somatic-only") synth_jobs.emplace_back(v, false);
            }
        }
        if (contains(modes, "Somatic")) {
            for (const auto& [v, _] : variants(true)) synth_jobs.emplace_back(v, true);
        }
    }
    std::vector<Result> results(paths.size() + synth_jobs.size());
    parallel_for(results.size(), jobs, [&](std::size_t i) {
        if (i < paths.size())
            results[i] = gen ? generate(paths[i], a) : compare(paths[i], a);
        else
            results[i] = synthetic(a, synth_jobs[i - paths.size()].first, synth_jobs[i - paths.size()].second);
    });
    std::ranges::sort(results, [](const Result& x, const Result& y) {
        return std::tie(x.name, x.status, x.detail) < std::tie(y.name, y.status, y.detail);
    });

    std::map<std::string, int> counts;
    int artefacts = 0;
    for (const auto& r : results) {
        ++counts[r.status];
        artefacts += static_cast<int>(
            std::ranges::count_if(r.detail, [](const std::string& d) { return d.starts_with("known R artefact"); }));
        if (gen) {
            std::println("{:<5} {}: {}", r.status, r.name, join(r.detail, "; "));
        } else if (r.status != "pass") {
            std::string up = r.status;
            std::ranges::transform(up, up.begin(), ::toupper);
            std::println("{:<5} {}", up, r.name);
            for (std::size_t i = 0; i < r.detail.size() && i < static_cast<std::size_t>(diff); ++i)
                std::println("        {}", r.detail[i]);
            if (r.detail.size() > static_cast<std::size_t>(diff))
                std::println("        ... {} more", r.detail.size() - static_cast<std::size_t>(diff));
        }
    }
    if (gen) {
        std::println("generated {}, skipped {}, failed {} -> {}", counts["ok"], counts["skip"], counts["fail"],
                     a.expected.string());
        return counts["fail"] ? 1 : 0;
    }
    const int total = counts["pass"] + counts["fail"];
    const auto label = a.stage + " [" + a.variant + "]" + (synth ? " + synthetic [all variants]" : "");
    std::println("VCF {}: {}/{} pass, {} skipped{}", label, counts["pass"], total, counts["skip"],
                 artefacts ? std::format(", {} known R artefact(s) accepted", artefacts) : "");
    return counts["fail"] == 0 ? 0 : 1;
}

}  // namespace testkit
