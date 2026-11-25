// golden: run omnivar on the golden cases and compare with the reference output.
//
//   <binary> --tsv <bed> -b <bam> [case args] [extra args] -z -c 1 -S 2 -E 3 -g 4 -G <fasta>.csv
//
// Expected output comes from the case file (lines after the header) or, with --expected DIR, from
// DIR/<case file name> (oracle output, `omnivar-testkit oracle --out`).
//
// Stages (shorthand for --expected + --extra-args):
//     full   docs/golden/cases                  no extra args
//     U      tests/golden/expected_U            -U        (realignment, no SV)
//     k0_U   tests/golden/expected_k0_U         -k 0 -U   (no realignment, no SV)
//
// Pass: exit code 0 and the same multiset of non-empty rows (order and trailing blank lines are
// not significant). EXACT = byte-identical, SORTED = identical after sorting.
//
// -p cases (M10): calls must equal the reference output without -p (tests/golden/expected_nopileup,
// `omnivar-testkit nopileup`); each rejected row carries its reasons in an extra last column, and
// the rejected rows are compared with a snapshot (tests/golden/expected_rejected; --update-rejected).

#include <algorithm>
#include <map>
#include <print>
#include <set>
#include <unistd.h>

#include "testkit.hpp"

namespace testkit {

namespace {

struct Stage {
    std::optional<fs::path> expected;
    std::string extra;
};

std::map<std::string, Stage> stages() {
    return {{"full", {std::nullopt, ""}},
            {"U", {golden_dir() / "expected_U", "-U"}},
            {"k0_U", {golden_dir() / "expected_k0_U", "-k 0 -U"}}};
}

const std::set<std::string>& reject_ids() {
    static const std::set<std::string> ids = [] {
        const char* reasons[] = {"REJ_AF",     "REJ_MinReads", "REJ_PMean", "REJ_Qual", "REJ_RefMQ", "REJ_Splice",
                                 "REJ_QRatio", "REJ_MQ",       "REJ_MSI1",  "REJ_MSI",  "REJ_Bias"};
        std::set<std::string> s;
        for (const char* r : reasons) {
            s.insert(r);
            s.insert(std::string(r) + "_T");
            s.insert(std::string(r) + "_N");
        }
        s.insert({"REJ_Masked_T", "REJ_Combine_T", "REJ_Combine_N"});
        return s;
    }();
    return ids;
}

struct Split {
    std::string calls;
    std::vector<std::string> rejected, bad;
};

// -p single-sample TSV -> (calls without the Filter column, rejected rows, rows with an unknown Filter)
Split split_rejected(std::string_view out) {
    Split s;
    for (const auto& line : lines_keepends(out)) {
        if (is_blank(line)) {
            s.calls += line;
            continue;
        }
        const auto l = rstrip_newline(line);
        const auto tab = l.rfind('\t');
        const std::string row = tab == std::string::npos ? "" : l.substr(0, tab);
        const std::string reason = tab == std::string::npos ? l : l.substr(tab + 1);
        bool known = !reason.empty();
        for (const auto& r : split(reason, ';')) known = known && reject_ids().contains(r);
        if (reason == "PASS")
            s.calls += row + "\n";
        else if (known)
            s.rejected.push_back(row + "\t" + reason);
        else
            s.bad.push_back(l);
    }
    return s;
}

struct Opts {
    std::string binary, stage, extra_args, out;
    fs::path cases, bams, fastas;
    std::optional<fs::path> expected;
    bool stage_expectations = false, verbose = false, update_rejected = false;
    int diff = 0, timeout = 600;
};

struct Result {
    std::string name, mode, status;
    std::vector<std::string> detail, cmd;
};

std::vector<std::string> set_minus(const std::vector<std::string>& a, const std::vector<std::string>& b) {
    const std::set<std::string> sb(b.begin(), b.end());
    std::set<std::string> out;
    for (const auto& x : a) {
        if (!sb.contains(x)) out.insert(x);
    }
    return {out.begin(), out.end()};
}

Result run_case(const fs::path& path, const Opts& a) {
    const auto c = parse_case(path);
    const bool pileup = case_has_arg(c, "-p") && a.stage_expectations;
    std::string expected = c.expected;
    if (pileup)
        expected = read_file(golden_dir() / "expected_nopileup" / a.stage / c.name);
    else if (a.expected)
        expected = read_file(*a.expected / c.name);

    const TempFile bed(".bed", case_bed_line(c, true));
    std::vector<std::string> cmd = {a.binary, "--tsv"};
    for (auto& w : words(bed.path() + " -b " + case_bams(c, a.bams) + " " + c.args + " " + a.extra_args +
                         " -z -c 1 -S 2 -E 3 -g 4 -G " + case_fasta_csv(c, a.fastas)))
        cmd.push_back(std::move(w));
    auto p = run(cmd, timeout_opt(a.timeout));
    if (p.timed_out) p.out.clear();
    if (!a.out.empty()) write_file(fs::path(a.out) / c.name, p.out);

    Split sp;
    std::string out = p.out;
    if (pileup && p.rc == 0 && !p.timed_out) {
        sp = split_rejected(out);
        out = sp.calls;
    }
    const auto got = sorted_rows(out), want = sorted_rows(expected);
    Result r{c.name, c.mode, "", {}, cmd};
    if (p.timed_out)
        r.status = "TIMEOUT";
    else if (p.rc != 0)
        r.status = "ERROR";
    else if (out == expected)
        r.status = "EXACT";
    else if (got == want)
        r.status = "SORTED";
    else
        r.status = "DIFF";
    if (!sp.bad.empty()) r.status = "DIFF";

    std::string snapshot_diff;
    if (pileup && p.rc == 0 && !p.timed_out) {
        const auto snap = golden_dir() / "expected_rejected" / a.stage / c.name;
        auto rej = sp.rejected;
        std::ranges::sort(rej);
        if (a.update_rejected) {
            fs::create_directories(snap.parent_path());
            std::string text;
            for (const auto& x : rej) text += x + "\n";
            write_file(snap, text);
        } else if (!fs::exists(snap)) {
            snapshot_diff = "no rejected-row snapshot " + snap.string() + " (run with --update-rejected)";
        } else {
            const auto want_rej = sorted_rows(read_file(snap));
            if (want_rej != rej)
                snapshot_diff = std::format("rejected rows differ from the snapshot: {} missing, {} unexpected",
                                            set_minus(want_rej, rej).size(), set_minus(rej, want_rej).size());
        }
        if (!snapshot_diff.empty()) r.status = "DIFF";
    }

    if (r.status == "ERROR") {
        auto d = std::format("rc={} {}", p.rc, last_line(p.err));
        while (!d.empty() && d.back() == ' ') d.pop_back();
        r.detail.push_back(d);
    }
    if ((r.status == "DIFF" || r.status == "ERROR" || r.status == "TIMEOUT") && (!want.empty() || !got.empty())) {
        const auto missing = set_minus(want, got), extra = set_minus(got, want);
        r.detail.push_back(std::format("expected {} rows, got {}; {} missing, {} unexpected", want.size(), got.size(),
                                       missing.size(), extra.size()));
        for (std::size_t i = 0; i < missing.size() && i < static_cast<std::size_t>(a.diff); ++i)
            r.detail.push_back("  - " + missing[i]);
        for (std::size_t i = 0; i < extra.size() && i < static_cast<std::size_t>(a.diff); ++i)
            r.detail.push_back("  + " + extra[i]);
    }
    if (!snapshot_diff.empty()) r.detail.push_back(snapshot_diff);
    if (!sp.bad.empty()) {
        const auto& b = sp.bad.front();
        r.detail.push_back(std::format("{} rows with an unknown Filter column, e.g. {}", sp.bad.size(),
                                       repr(b.size() > 80 ? b.substr(b.size() - 80) : b)));
    }
    if (pileup && (a.verbose || r.status != "EXACT"))
        r.detail.push_back(
            std::format("-p: {} rejected rows (calls compared with the jar without -p)", sp.rejected.size()));
    return r;
}

}  // namespace

int cmd_golden(Args& args) {
    Opts a;
    a.binary = fs::absolute(args.required("--binary")).string();
    a.cases = args.value_or("--cases", golden_cases().string());
    a.bams = args.value_or("--bams", golden_bams().string());
    a.fastas = args.value_or("--fastas", golden_fastas().string());
    a.stage = args.value_or("--stage", "full");
    const auto expected = args.value("--expected");
    const auto extra = args.value("--extra-args");
    auto modes = args.values("--mode");
    const auto match = args.values("-k");
    a.out = args.value_or("--out", "");
    a.diff = args.int_or("--diff", 0);
    a.verbose = args.flag("-v") || args.flag("--verbose");
    a.timeout = args.int_or("--timeout", 600);
    a.update_rejected = args.flag("--update-rejected");
    const int jobs = args.int_or("-j", default_jobs());
    args.finish();

    const auto st = stages();
    if (!st.contains(a.stage)) die("golden: --stage must be full, U or k0_U");
    a.stage_expectations = !expected && !extra;
    a.expected = expected ? std::optional<fs::path>(*expected) : st.at(a.stage).expected;
    a.extra_args = extra.value_or(st.at(a.stage).extra);
    if (::access(a.binary.c_str(), X_OK) != 0) die("not an executable: " + a.binary);
    if (!a.out.empty()) fs::create_directories(a.out);

    const bool modes_given = !modes.empty();
    if (modes.empty()) modes = {"Simple", "Somatic"};
    std::vector<fs::path> files;
    int skipped = 0;
    for (const auto& f : case_files(a.cases)) {
        const auto c = parse_case(f);
        if (!contains(modes, c.mode)) {
            ++skipped;
            continue;
        }
        if (!match.empty() && std::ranges::none_of(match, [&](const std::string& m) { return c.name.contains(m); }))
            continue;
        files.push_back(f);
    }
    if (files.empty()) die("no cases selected");

    std::vector<Result> results(files.size());
    parallel_for(files.size(), jobs, [&](std::size_t i) { results[i] = run_case(files[i], a); });

    auto passing = [](const Result& r) { return r.status == "EXACT" || r.status == "SORTED"; };
    for (const auto& r : results) {
        if (passing(r) && !a.verbose) continue;
        std::println("{:<8} {:<8} {}", r.status, r.mode, r.name);
        for (const auto& d : r.detail) std::println("         {}", d);
        if (a.verbose && !passing(r)) std::println("         cmd: {}", join(r.cmd, " "));
    }
    std::println("");
    std::map<std::string, std::array<int, 3>> by_mode;
    for (const auto& r : results) {
        auto& t = by_mode[r.mode];
        ++t[0];
        t[1] += passing(r);
        t[2] += r.status == "EXACT";
    }
    int total = 0;
    for (const auto& [mode, t] : by_mode) {
        std::println("{:<8} {:3}/{:<3} pass ({} exact)", mode, t[1], t[0], t[2]);
        total += t[1];
    }
    const bool custom = expected.has_value() || extra.has_value();
    std::println("{:<8} {:3}/{:<3} pass   [stage {}, extra args: '{}']", "Total", total, results.size(),
                 custom ? "custom" : a.stage, a.extra_args);
    if (skipped && !modes_given)
        std::println("({} Amplicon cases skipped: amplicon mode is out of scope; --mode Amplicon runs them)", skipped);
    return total == static_cast<int>(results.size()) ? 0 : 1;
}

}  // namespace testkit
