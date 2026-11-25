// oracle: run VarDictJava 1.8.3 (through tests/oracle/OracleMain.java) on the golden cases.
//
// Each case file: line 1 is a CSV header
//     Mode,fasta,bam,chr,start,end[,args]                  (Simple / Somatic)
//     Mode,fasta,bam,chr,start,end,ampStart,ampEnd[,args]  (Amplicon)
// Somatic cases give two BAMs as "tumor.bam|normal.bam"; the remaining lines are the jar's
// expected stdout. The argument list mirrors IntegrationTest.buildArgs:
//     <bed> -b <bam> [args] -z -c 1 -S 2 -E 3 -g 4 -G <fasta>.csv
// --check compares with the expected lines: EXACT (byte-identical) or SORTED (same non-empty rows
// in any order, the parity criterion). --out writes each run's stdout to DIR/<case file name>.
//
// nopileup: expected output for the golden cases that use -p (M10). Here -p means "also report the
// rejected variants", not the jar's pileup, so the calls of a -p run must equal the jar's output
// for the same case without -p. Runs the jar on those cases with -p removed, at every stage, into
// tests/golden/expected_nopileup/<stage>/<case file name>.

#include <map>
#include <print>

#include "testkit.hpp"

namespace testkit {

namespace {

struct OracleOpts {
    fs::path cases, bams, fastas, jars, runner, out;
    std::string extra;
    bool check = false;
    int timeout = 600, jobs = 4;
};

struct OracleResult {
    std::string name, mode;
    int rc = 0;
    bool exact = false, sorted_eq = false;
    std::string err;
};

OracleResult oracle_case(const fs::path& path, const OracleOpts& a) {
    const auto c = parse_case(path);
    const TempFile bed(".bed", case_bed_line(c, true));
    const auto csv = case_fasta_csv(c, a.fastas);
    std::vector<std::string> cmd = {"java", "-cp", (a.jars / "*").string(), a.runner.string(), csv};
    for (auto& w : words(bed.path() + " -b " + case_bams(c, a.bams) + " " + c.args + " " + a.extra +
                         " -z -c 1 -S 2 -E 3 -g 4 -G " + csv))
        cmd.push_back(std::move(w));
    auto p = run(cmd, timeout_opt(a.timeout));
    OracleResult r{c.name, c.mode, p.timed_out ? -1 : p.rc, false, false, ""};
    if (p.timed_out) {
        p.out.clear();
        r.err = "timeout";
    } else if (p.rc) {
        r.err = last_line(p.err);
    }
    if (!a.out.empty()) write_file(a.out / c.name, p.out);
    r.exact = p.out == c.expected;
    r.sorted_eq = sorted_rows(p.out) == sorted_rows(c.expected);
    return r;
}

int run_oracle(const OracleOpts& a) {
    if (!a.out.empty()) fs::create_directories(a.out);
    const auto files = case_files(a.cases);
    std::vector<OracleResult> results(files.size());
    parallel_for(files.size(), a.jobs, [&](std::size_t i) { results[i] = oracle_case(files[i], a); });
    int failed = 0;
    for (const auto& r : results) {
        std::string status = r.exact ? "EXACT" : (r.sorted_eq ? "SORTED" : "DIFF");
        if (r.rc != 0) status = std::format("ERROR rc={} {}", r.rc, r.err);
        if (a.check && status != "EXACT" && status != "SORTED") ++failed;
        if (a.check || r.rc != 0) std::println("{:<8} {:<8} {}", status, r.mode, r.name);
    }
    if (a.check) {
        std::map<std::string, std::array<int, 3>> by;
        for (const auto& r : results) {
            auto& t = by[r.mode];
            ++t[0];
            t[1] += r.exact;
            t[2] += r.sorted_eq;
        }
        for (const auto& [mode, t] : by)
            std::println("{}: {} cases, {} exact, {} identical after sorting", mode, t[0], t[1], t[2]);
    }
    return failed ? 1 : 0;
}

}  // namespace

int cmd_oracle(Args& args) {
    OracleOpts a;
    a.cases = args.value_or("--cases", golden_cases().string());
    a.bams = args.value_or("--bams", golden_bams().string());
    a.fastas = args.value_or("--fastas", golden_fastas().string());
    a.jars = args.value_or("--jars", jars_dir().string());
    a.runner = args.value_or("--runner", (root() / "tests/oracle/OracleMain.java").string());
    a.out = args.value_or("--out", "");
    a.extra = args.value_or("--extra-args", "");
    a.check = args.flag("--check");
    a.timeout = args.int_or("--timeout", 600);
    a.jobs = args.int_or("-j", default_jobs());
    args.finish();
    return run_oracle(a);
}

int cmd_nopileup(Args& args) {
    const int jobs = args.int_or("-j", 8);
    args.finish();
    const TempDir tmp;
    int n = 0;
    for (const auto& path : case_files(golden_cases(), false)) {
        const auto c = parse_case(path);
        if ((c.mode != "Simple" && c.mode != "Somatic") || !case_has_arg(c, "-p")) continue;
        const auto text = read_file(path);
        auto header = split(text.substr(0, text.find('\n')), ',');
        std::vector<std::string> kept;  // the args are field 7 for Simple / Somatic
        for (auto& w : words(header.at(6))) {
            if (w != "-p") kept.push_back(std::move(w));
        }
        header[6] = join(kept, " ");
        write_file(tmp.path() / c.name, join(header, ",") + "\n");
        ++n;
    }
    const std::pair<const char*, const char*> stages[] = {{"full", ""}, {"U", "-U"}, {"k0_U", "-k 0 -U"}};
    for (const auto& [stage, extra] : stages) {
        const auto out = golden_dir() / "expected_nopileup" / stage;
        fs::remove_all(out);
        fs::create_directories(out);
        OracleOpts a;
        a.cases = tmp.path();
        a.bams = golden_bams();
        a.fastas = golden_fastas();
        a.jars = jars_dir();
        a.runner = root() / "tests/oracle/OracleMain.java";
        a.out = out;
        a.extra = extra;
        a.jobs = jobs;
        run_oracle(a);
        const auto written = std::distance(fs::directory_iterator(out), fs::directory_iterator{});
        std::println("{}: {}/{} cases -> {}", stage, written, n, out.string());
    }
    return 0;
}

}  // namespace testkit
