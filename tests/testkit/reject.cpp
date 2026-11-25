// reject-invariant (M10): adding -p must not change the calls.
//
// For every Simple and Somatic golden case (except splicing / CRISPR), the binary runs with and
// without -p, for the TSV (--tsv) and the default VCF:
//     TSV  the rows whose last column is PASS, without that column, equal the run without -p;
//          every rejected row has a REF allele
//     both stderr is the same as without -p (no exception is raised for a rejected variant)
//     VCF  the records whose FILTER does not start with REJ_ equal the run without -p, the header
//          differs only by ##FILTER=<ID=REJ_ lines, and every rejected FILTER ID is declared
// Also prints how often each reason occurs over all cases.

#include <algorithm>
#include <map>
#include <print>
#include <set>

#include "testkit.hpp"

namespace testkit {

namespace {

struct Result {
    std::string name;
    std::optional<bool> ok;  // nothing: skipped
    std::vector<std::string> problems;
    std::map<std::string, int> reasons;
};

// "##FILTER=<ID=REJ_\w+," -> the ID
std::optional<std::string> declared_reject(std::string_view l) {
    constexpr std::string_view head = "##FILTER=<ID=";
    if (!l.starts_with(std::string(head) + "REJ_")) return std::nullopt;
    auto i = head.size() + 4;
    while (i < l.size() && (std::isalnum(static_cast<unsigned char>(l[i])) || l[i] == '_')) ++i;
    if (i == head.size() + 4 || i >= l.size() || l[i] != ',') return std::nullopt;
    return std::string(l.substr(head.size(), i - head.size()));
}

std::set<std::string> line_set(std::string_view s) {
    const auto v = lines(s);
    return {v.begin(), v.end()};
}

Result check(const fs::path& path, const std::string& binary, const std::string& extra, int timeout) {
    const auto c = parse_case(path);
    std::vector<std::string> args;
    for (auto& w : words(c.args)) {
        if (w != "-p") args.push_back(std::move(w));
    }
    if ((c.mode != "Simple" && c.mode != "Somatic") || contains(args, "-i") || contains(args, "-J"))
        return {c.name, {}, {}, {}};
    const TempFile bed(".bed", c.chr + "\t" + c.start + "\t" + c.end + "\tt\n");
    std::vector<std::string> base = {bed.path(), "-b", case_bams(c, golden_bams())};
    base.insert(base.end(), args.begin(), args.end());
    for (auto& w : words(extra)) base.push_back(std::move(w));
    for (auto& w : words("-z -c 1 -S 2 -E 3 -g 4 -G")) base.push_back(std::move(w));
    base.push_back(case_fasta_csv(c, golden_fastas()));
    auto go = [&](std::vector<std::string> pre) {
        std::vector<std::string> v = {binary};
        v.insert(v.end(), pre.begin(), pre.end());
        v.insert(v.end(), base.begin(), base.end());
        return run(v, timeout_opt(timeout));
    };
    const auto t0 = go({"--tsv"}), t1 = go({"--tsv", "-p"}), v0 = go({}), v1 = go({"-p"});
    Result r{c.name, true, {}, {}};
    for (const auto* p : {&t0, &t1, &v0, &v1}) {
        if (p->rc || p->timed_out) {
            r.ok = false;
            r.problems.push_back(std::format("exit {}: {}", p->rc, tail(p->err, 200)));
            return r;
        }
    }
    std::string calls;
    for (const auto& line : lines_keepends(t1.out)) {
        const auto l = rstrip_newline(line);
        const auto tab = l.rfind('\t');
        const std::string row = tab == std::string::npos ? "" : l.substr(0, tab);
        const std::string flt = tab == std::string::npos ? l : l.substr(tab + 1);
        if (flt == "PASS") {
            calls += row + "\n";
        } else if (flt.starts_with("REJ_")) {
            for (const auto& x : split(flt, ';')) ++r.reasons[x];
            const auto cols = split(row, '\t');
            if (cols.size() <= 5 || cols[5].empty()) {
                r.problems.push_back("TSV: rejected row without REF: " + repr(row.substr(0, 120)));
                break;
            }
        } else {
            r.problems.push_back("TSV row with Filter " + repr(flt));
            break;
        }
    }
    if (t1.err != t0.err || v1.err != v0.err) {
        std::string first = "(order)";
        const auto base_err = line_set(t0.err);
        for (const auto& l : lines(t1.err)) {
            if (!base_err.contains(l)) {
                first = l;
                break;
            }
        }
        r.problems.push_back("stderr differs from the run without -p: " + first.substr(0, 150));
    }
    if (calls != t0.out) r.problems.push_back("TSV: calls with -p differ from the run without -p");
    std::vector<std::string> head0, head1, recs0, recs1;
    std::set<std::string> declared;
    for (auto& l : lines(v0.out)) (l.starts_with('#') ? head0 : recs0).push_back(l);
    for (auto& l : lines(v1.out)) {
        if (!l.starts_with('#')) continue;
        if (auto id = declared_reject(l)) declared.insert(*id);
        if (!l.starts_with("##FILTER=<ID=REJ_")) head1.push_back(l);
    }
    if (head1 != head0) r.problems.push_back("VCF: header differs beyond the REJ_ FILTER lines");
    for (const auto& l : lines(v1.out)) {
        if (l.starts_with('#')) continue;
        const auto cols = split(l, '\t');
        const std::string flt = cols.size() > 6 ? cols[6] : "";
        if (flt.starts_with("REJ_")) {
            std::vector<std::string> undeclared;
            for (const auto& t : split(flt, ';')) {
                if (t.starts_with("REJ_") && !declared.contains(t) && !contains(undeclared, t)) undeclared.push_back(t);
            }
            if (!undeclared.empty()) {
                std::ranges::sort(undeclared);
                r.problems.push_back("VCF: undeclared FILTER [" + join(undeclared, ", ") + "]");
                break;
            }
        } else {
            recs1.push_back(l);
        }
    }
    if (recs1 != recs0)
        r.problems.push_back("VCF: records with -p (rejected ones left out) differ from the run without -p");
    r.ok = r.problems.empty();
    return r;
}

}  // namespace

int cmd_reject(Args& args) {
    const auto binary = fs::absolute(args.required("--binary")).string();
    const auto extra = args.value_or("--extra-args", "");
    const int timeout = args.int_or("--timeout", 600);
    const int jobs = args.int_or("-j", 8);
    args.finish();
    const auto paths = case_files(golden_cases(), false);
    std::vector<Result> results(paths.size());
    parallel_for(paths.size(), jobs, [&](std::size_t i) { results[i] = check(paths[i], binary, extra, timeout); });
    std::map<std::string, int> total;
    int n = 0, failed = 0;
    for (const auto& r : results) {
        if (!r.ok) continue;
        ++n;
        for (const auto& [k, v] : r.reasons) total[k] += v;
        if (!*r.ok) {
            ++failed;
            std::println("FAIL  {}", r.name);
            for (const auto& p : r.problems) std::println("      {}", p);
        }
    }
    std::vector<std::string> parts;
    for (const auto& [k, v] : total) parts.push_back(std::format("{} {}", k, v));
    std::println("reasons: {}", join(parts, ", "));
    std::println("-p invariant: {}/{} cases keep their calls [extra args: '{}']", n - failed, n, extra);
    return failed ? 1 : 0;
}

}  // namespace testkit
