// threads: `senri -th N` must write exactly what `-th 1` writes.
//
// The regression cases have one region each, which never exercises the worker pool. Every Simple and
// Somatic case whose region is at least --min-length bp is split into --chunk bp BED regions, and
// the binary runs on that BED with -th 1 and -th N, for the TSV (--tsv) and the default VCF. The
// outputs must be byte-identical (output is in region order). Then a fatal region (outside the
// reference slice) is inserted halfway through: -th N must stop where -th 1 stops, with the same
// stdout, stderr and exit status.

#include <print>

#include "testkit.hpp"

namespace testkit {

namespace {

struct Opts {
    std::string binary;
    int threads = 4, chunk = 500, min_length = 2000, timeout = 900;
};

struct Result {
    std::string name, status, detail;
    long regions = 0;
};

Result check(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    if (c.mode != "Simple" && c.mode != "Somatic") return {c.name, "skip", c.mode + " mode"};
    const auto args = words(c.args);
    if (contains(args, "-i") || contains(args, "-J")) return {c.name, "skip", "splicing / CRISPR"};
    const long start = std::stol(c.start), end = std::stol(c.end);
    if (end - start < a.min_length) return {c.name, "skip", "short region"};
    std::vector<std::string> bed_lines;
    for (long s = start; s < end; s += a.chunk)
        bed_lines.push_back(std::format("{}\t{}\t{}\tchunk{}\n", c.chr, s, std::min(s + a.chunk, end), s));
    const TempFile bed(".bed", join(bed_lines, ""));
    std::vector<std::string> base = {bed.path(), "-b", case_bams(c, bams_dir())};
    base.insert(base.end(), args.begin(), args.end());
    for (auto& w : words("-c 1 -S 2 -E 3 -g 4 -G")) base.push_back(std::move(w));
    base.push_back(case_fasta_csv(c, fastas_dir()));
    auto cmd = [&](bool tsv, int th) {
        std::vector<std::string> v = {a.binary};
        if (tsv) v.emplace_back("--tsv");
        v.insert(v.end(), base.begin(), base.end());
        v.emplace_back("-th");
        v.push_back(std::to_string(th));
        return v;
    };
    for (const bool tsv : {true, false}) {
        const auto r1 = run(cmd(tsv, 1), timeout_opt(a.timeout));
        const auto rn = run(cmd(tsv, a.threads), timeout_opt(a.timeout));
        if (r1.rc != rn.rc || r1.out != rn.out)
            return {c.name, "fail",
                    std::format("{}: -th 1 exit {}, -th {} exit {}, outputs differ", tsv ? "TSV" : "VCF", r1.rc,
                                a.threads, rn.rc)};
    }
    // A fatal region halfway through ("No region: ..."): -th N must stop at the same region.
    auto with_fatal = bed_lines;
    with_fatal.insert(with_fatal.begin() + static_cast<long>(bed_lines.size() / 2),
                      std::format("{}\t{}\t{}\tfatal\n", c.chr, end + 50000, end + 50300));
    bed.write(join(with_fatal, ""));
    const auto r1 = run(cmd(true, 1), timeout_opt(a.timeout));
    const auto rn = run(cmd(true, a.threads), timeout_opt(a.timeout));
    if (r1.rc == 0) return {c.name, "fail", "fatal region: -th 1 did not fail"};
    if (r1.rc != rn.rc || r1.out != rn.out || r1.err != rn.err)
        return {c.name, "fail",
                std::format("fatal region: -th 1 exit {}, -th {} exit {}, output differs", r1.rc, a.threads, rn.rc)};
    return {c.name, "pass", "", static_cast<long>(bed_lines.size())};
}

}  // namespace

int cmd_threads(Args& args) {
    Opts a;
    a.binary = fs::absolute(args.required("--binary")).string();
    a.threads = args.int_or("-th", 4);
    a.chunk = args.int_or("--chunk", 500);
    a.min_length = args.int_or("--min-length", 2000);
    a.timeout = args.int_or("--timeout", 900);
    const int jobs = args.int_or("-j", 4);
    args.finish();
    const auto paths = case_files(cases_dir(), false);
    std::vector<Result> results(paths.size());
    parallel_for(paths.size(), jobs, [&](std::size_t i) { results[i] = check(paths[i], a); });
    int pass = 0, fail = 0, skip = 0;
    long regions = 0;
    for (const auto& r : results) {
        if (r.status == "pass") {
            ++pass;
            regions += r.regions;
        } else if (r.status == "fail") {
            ++fail;
            std::println("FAIL  {}: {}", r.name, r.detail);
        } else {
            ++skip;
        }
    }
    std::println("threads -th {} vs -th 1: {}/{} cases identical ({} regions), {} skipped", a.threads, pass,
                 pass + fail, regions, skip);
    return fail == 0 ? 0 : 1;
}

}  // namespace testkit
