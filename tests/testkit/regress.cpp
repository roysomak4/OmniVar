// regress: omnivar's output on the regression cases against its own checked-in output.
//
//   omnivar-testkit regress --binary B [--part tsv,vcf,synthetic] [--mode Simple|Somatic] [-k substr]
//       [--update] [--diff N] [-v] [-j N] [--timeout S]
//
// Parts:
//   tsv        tests/regression/cases/<case>.txt: line 1 is the case spec, the rest the expected
//              `--tsv` output. Pass: the same multiset of non-empty rows (EXACT = byte-identical,
//              SORTED = identical after sorting).
//   vcf        tests/regression/vcf/<case>.vcf: the expected VCF (default output) of the same case.
//              Pass: identical lines, except ##source (version) and ##reference (a local path).
//   synthetic  tests/regression/synthetic/{single,paired}.tsv: crafted rows (many contigs, SVs,
//              clusters, edge values) through `omnivar --dev-tsv-to-vcf`, for each VCF option set;
//              expected synthetic/<input>.<variant>.vcf.
//
// Every case runs as
//   <binary> [--tsv] <bed> -b <bam> <case args> -z -c 1 -S 2 -E 3 -g 4 -G <fasta slice>.csv
// and every VCF must also pass `bcftools view` when bcftools is installed.
//
// The expectations are omnivar's own output, so a deliberate change of behaviour is recorded with
// --update (review the diff before committing it). Version 0.1.0 produced exactly the output of
// VarDictJava 1.8.3 and its R and Perl scripts on these cases (see docs/release-notes/v0.1.0.md).

#include <algorithm>
#include <functional>
#include <map>
#include <print>
#include <set>
#include <tuple>
#include <unistd.h>

#include "testkit.hpp"

namespace testkit {

namespace {

struct Opts {
    std::string binary;
    bool update = false, verbose = false;
    int diff = 0, timeout = 600;
};

struct Result {
    std::string part, name, status;  // EXACT SORTED DIFF ERROR TIMEOUT UPDATED
    std::vector<std::string> detail;
};

bool passing(const Result& r) { return r.status == "EXACT" || r.status == "SORTED" || r.status == "UPDATED"; }

std::vector<std::string> set_minus(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    const std::set<std::string> sb(b.begin(), b.end());
    std::set<std::string> out;
    for (const auto& x : a) {
        if (!sb.contains(x)) out.insert(x);
    }
    return {out.begin(), out.end()};
}

void add_row_diff(Result& r, const std::vector<std::string>& want, const std::vector<std::string>& got, int n) {
    const auto missing = set_minus(want, got), extra = set_minus(got, want);
    r.detail.push_back(std::format("expected {} rows, got {}; {} missing, {} unexpected", want.size(), got.size(),
                                   missing.size(), extra.size()));
    for (std::size_t i = 0; i < missing.size() && i < static_cast<std::size_t>(n); ++i)
        r.detail.push_back("  - " + missing[i]);
    for (std::size_t i = 0; i < extra.size() && i < static_cast<std::size_t>(n); ++i)
        r.detail.push_back("  + " + extra[i]);
}

// The VCF without the lines that depend on the version or the machine.
std::string normalise_vcf(std::string_view vcf) {
    std::string out;
    for (const auto& l : lines_keepends(vcf)) {
        if (!l.starts_with("##source=") && !l.starts_with("##reference=")) out += l;
    }
    return out;
}

std::optional<std::string> bcftools_check(const std::string& vcf) {
    static const auto bcftools = which("bcftools");
    if (!bcftools) return std::nullopt;
    const auto p = run({*bcftools, "view", "-Ov", "-o", "/dev/null", "-"}, input_opt(vcf));
    const auto msg = strip(p.err);
    if (p.rc || !msg.empty()) return msg.empty() ? std::format("exit {}", p.rc) : lines(msg).front();
    return std::nullopt;
}

std::vector<std::string> case_cmd(const Opts& a, const Case& c, const std::string& bed, bool tsv) {
    std::vector<std::string> cmd = {a.binary};
    if (tsv) cmd.emplace_back("--tsv");
    for (auto& w : words(bed + " -b " + case_bams(c, bams_dir()) + " " + c.args + " -z -c 1 -S 2 -E 3 -g 4 -G " +
                         case_fasta_csv(c, fastas_dir())))
        cmd.push_back(std::move(w));
    return cmd;
}

// A failed run: ERROR or TIMEOUT with the last line of stderr.
std::optional<Result> failed(const std::string& part, const std::string& name, const RunResult& p) {
    if (p.timed_out) return Result{part, name, "TIMEOUT", {}};
    if (p.rc != 0) return Result{part, name, "ERROR", {std::format("rc={} {}", p.rc, last_line(p.err))}};
    return std::nullopt;
}

Result tsv_case(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    const TempFile bed(".bed", case_bed_line(c, false));
    const auto p = run(case_cmd(a, c, bed.path(), true), timeout_opt(a.timeout));
    if (auto f = failed("tsv", c.name, p)) return *f;
    if (a.update) {
        const auto text = read_file(path);
        write_file(path, text.substr(0, text.find('\n') + 1) + p.out);
        return {"tsv", c.name, "UPDATED", {}};
    }
    Result r{"tsv", c.name, "", {}};
    const auto got = sorted_rows(p.out), want = sorted_rows(c.expected);
    r.status = p.out == c.expected ? "EXACT" : got == want ? "SORTED" : "DIFF";
    if (r.status == "DIFF") add_row_diff(r, want, got, a.diff);
    return r;
}

Result vcf_compare(const std::string& part, const std::string& name, const fs::path& expected_path,
                   const std::string& out, const Opts& a) {
    const auto vcf = normalise_vcf(out);
    if (a.update) {
        fs::create_directories(expected_path.parent_path());
        write_file(expected_path, vcf);
        return {part, name, "UPDATED", {}};
    }
    Result r{part, name, "", {}};
    if (!fs::exists(expected_path)) {
        r.status = "DIFF";
        r.detail.push_back("no expected VCF " + expected_path.string() + " (run with --update)");
        return r;
    }
    const auto want = read_file(expected_path);
    r.status = vcf == want ? "EXACT" : "DIFF";
    if (r.status == "DIFF") add_row_diff(r, lines(want), lines(vcf), a.diff);
    // (-i writes splice-junction counts instead of a VCF)
    if (auto bad = out.starts_with("##fileformat=VCF") ? bcftools_check(out) : std::nullopt) {
        r.status = "DIFF";
        r.detail.push_back("bcftools: " + *bad);
    }
    return r;
}

Result vcf_case(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    const TempFile bed(".bed", case_bed_line(c, false));
    const auto p = run(case_cmd(a, c, bed.path(), false), timeout_opt(a.timeout));
    if (auto f = failed("vcf", c.name, p)) return *f;
    auto name = c.name;
    name.replace(name.size() - 4, 4, ".vcf");
    return vcf_compare("vcf", c.name, vcf_dir() / name, p.out, a);
}

const std::string kThresholds =
    "--filter-min-depth 10 --filter-min-var-depth 4 --filter-min-af 0.05 --filter-min-pmean 10 "
    "--filter-min-qual 30 --filter-min-mq 20 --filter-min-sn 3 --filter-max-nm 2.5 --filter-msi 6 "
    "--filter-hom-af 0.3 ";

// (input, variant) -> omnivar VCF options
std::vector<std::tuple<std::string, std::string, std::string>> synthetic_jobs() {
    return {
        {"single", "default", ""},
        {"single", "top", "--vcf-top-allele-only"},
        {"single", "cluster", "--filter-cluster-bp 10"},
        {"single", "thresholds", kThresholds + "--no-filter-pstd"},
        {"paired", "default", ""},
        {"paired", "top", "--vcf-top-allele-only"},
        {"paired", "cluster", "--filter-cluster-bp 10"},
        {"paired", "thresholds", kThresholds + "--filter-max-pvalue 0.01"},
        {"paired", "somatic-only", "--vcf-somatic-only"},
    };
}

Result synthetic_case(const std::string& input, const std::string& variant, const std::string& options, const Opts& a) {
    const auto name = input + "." + variant;
    std::vector<std::string> cmd = {a.binary, "--dev-tsv-to-vcf", "-"};
    for (auto& w : words(options)) cmd.push_back(std::move(w));
    const auto p = run(cmd, input_opt(read_file(synthetic_dir() / (input + ".tsv")), a.timeout));
    if (auto f = failed("synthetic", name, p)) return *f;
    return vcf_compare("synthetic", name, synthetic_dir() / (name + ".vcf"), p.out, a);
}

}  // namespace

int cmd_regress(Args& args) {
    Opts a;
    a.binary = fs::absolute(args.required("--binary")).string();
    const auto parts_arg = args.value_or("--part", "tsv,vcf,synthetic");
    auto modes = args.values("--mode");
    const auto match = args.values("-k");
    a.update = args.flag("--update");
    a.diff = args.int_or("--diff", 0);
    a.verbose = args.flag("-v") || args.flag("--verbose");
    a.timeout = args.int_or("--timeout", 600);
    const int jobs = args.int_or("-j", default_jobs());
    args.finish();
    if (::access(a.binary.c_str(), X_OK) != 0) die("not an executable: " + a.binary);
    const auto parts = split(parts_arg, ',');
    for (const auto& p : parts) {
        if (p != "tsv" && p != "vcf" && p != "synthetic") die("regress: --part takes tsv, vcf and synthetic");
    }
    if (modes.empty()) modes = {"Simple", "Somatic"};

    std::vector<fs::path> files;
    for (const auto& f : case_files(cases_dir())) {
        const auto c = parse_case(f);
        if (!contains(modes, c.mode)) continue;
        if (!match.empty() && std::ranges::none_of(match, [&](const std::string& m) { return c.name.contains(m); }))
            continue;
        files.push_back(f);
    }
    std::ranges::sort(files);

    // one job per (part, case)
    std::vector<std::function<Result()>> work;
    for (const auto& f : files) {
        if (contains(parts, "tsv")) work.emplace_back([&a, f] { return tsv_case(f, a); });
        if (contains(parts, "vcf")) work.emplace_back([&a, f] { return vcf_case(f, a); });
    }
    if (contains(parts, "synthetic") && match.empty()) {
        for (const auto& [input, variant, options] : synthetic_jobs()) {
            const bool paired = input == "paired";
            if (!contains(modes, paired ? "Somatic" : "Simple")) continue;
            work.emplace_back([&a, input, variant, options] { return synthetic_case(input, variant, options, a); });
        }
    }
    if (work.empty()) die("regress: nothing selected");

    std::vector<Result> results(work.size());
    parallel_for(work.size(), jobs, [&](std::size_t i) { results[i] = work[i](); });

    for (const auto& r : results) {
        if (passing(r) && !a.verbose) continue;
        std::println("{:<8} {:<9} {}", r.status, r.part, r.name);
        for (const auto& d : r.detail) std::println("         {}", d);
    }
    std::println("");
    std::map<std::string, std::pair<int, int>> by_part;
    int failures = 0;
    for (const auto& r : results) {
        auto& [pass, total] = by_part[r.part];
        ++total;
        pass += passing(r);
        failures += !passing(r);
    }
    for (const auto& [part, t] : by_part)
        std::println("{:<10} {:3}/{:<3} {}", part, t.first, t.second, a.update ? "updated" : "pass");
    if (a.update) std::println("Expectations rewritten: review `git diff tests/regression` before committing.");
    return failures == 0 ? 0 : 1;
}

}  // namespace testkit
