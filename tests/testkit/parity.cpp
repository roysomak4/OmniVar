// Stage parity against the 1.8.3 jar (needs Java).
//
// cli-parity: the C++ command-line parser builds the same Configuration as the jar. Each command
// line goes through the jar's CmdParser (tests/oracle/ConfigDumpMain.java) and through
// `omnivar --dev-dump-config`; the sorted "field=value" dumps must be identical, and command lines
// the jar rejects must be rejected with the same message. Command lines come from
// tests/cli/cases.tsv (arguments separated by TAB) and, with --golden, from the golden cases.
//
// inputs-parity (M2): what the jar and omnivar derive from their inputs, through
// tests/oracle/InputsDumpMain.java and `omnivar --dev-dump-inputs`: sample names, amplicon
// parameters, regions, reference bases and seed index (digests), and per region/BAM the reads
// returned by RecordPreprocessor (counts and an order-sensitive digest). Extra runs on other data
// are given after "--" as the caller argument list.

#include <algorithm>
#include <map>
#include <print>

#include "testkit.hpp"

namespace testkit {

namespace {

std::vector<std::vector<std::string>> golden_command_lines() {
    std::vector<std::vector<std::string>> out;
    for (const auto& path : case_files(golden_cases())) {
        const auto c = parse_case(path);
        std::vector<std::string> bams;
        for (const auto& b : split(c.bam, '|')) bams.push_back("docs/golden/bams/" + b);
        const auto csv = "docs/golden/fastas/" + c.fasta + ".csv";
        for (const char* extra : {"", "-k 0 -U", "-U"})
            out.push_back(words("case.bed -b " + join(bams, "|") + " " + c.args + " " + extra +
                                " -z -c 1 -S 2 -E 3 -g 4 -G " + csv));
    }
    return out;
}

std::vector<std::string> java_dumps(const std::vector<std::vector<std::string>>& cmdlines) {
    std::string input;
    for (const auto& a : cmdlines) input += join(a, "\t") + "\n";
    const auto p = run({"java", "-XX:-UsePerfData", "-Xlog:disable", "-cp", jar_classpath(),
                        (root() / "tests/oracle/ConfigDumpMain.java").string()},
                       input_opt(input, 3600));
    if (p.rc) die("ConfigDumpMain failed: " + tail(p.err, 300));
    std::vector<std::string> out;
    std::size_t pos = p.out.find("### ");
    while (pos != std::string::npos) {
        const auto next = p.out.find("### ", pos + 4);
        const auto block = p.out.substr(pos + 4, next == std::string::npos ? std::string::npos : next - pos - 4);
        const auto nl = block.find('\n');
        out.push_back(nl == std::string::npos ? "" : block.substr(nl + 1));
        pos = next;
    }
    return out;
}

std::string cpp_dump(const std::string& binary, const std::vector<std::string>& args) {
    std::vector<std::string> cmd = {binary, "--dev-dump-config"};
    cmd.insert(cmd.end(), args.begin(), args.end());
    const auto p = run(cmd);
    if (p.rc != 0) {
        auto msg = strip(p.err).empty() ? std::format("rc={}", p.rc) : last_line(p.err);
        if (msg.starts_with("omnivar: ")) msg = msg.substr(9);
        return "ERROR " + msg + "\n";
    }
    return p.out;
}

// "ERROR SomeException: message" -> "ERROR message"
std::string normalise_java(const std::string& text) {
    if (!text.starts_with("ERROR ")) return text;
    const auto rest = text.substr(6);
    const auto colon = rest.find(": ");
    return "ERROR " + (colon == std::string::npos ? rest : rest.substr(colon + 2));
}

std::string field(const std::string& line) { return line.substr(0, line.find('=')); }

}  // namespace

int cmd_cli(Args& args) {
    const auto binary = fs::absolute(args.required("--binary")).string();
    const auto cases = args.value_or("--cases", (root() / "tests/cli/cases.tsv").string());
    const bool golden = args.flag("--golden");
    const bool verbose = args.flag("-v") || args.flag("--verbose");
    args.finish();
    std::vector<std::vector<std::string>> cmdlines;
    for (const auto& l : lines(read_file(cases))) {
        if (!l.empty()) cmdlines.push_back(split(l, '\t'));
    }
    if (golden) {
        for (auto& g : golden_command_lines()) cmdlines.push_back(std::move(g));
    }
    // Deliberate difference from upstream: -p keeps -f and -r as given, where the jar's
    // pileup sets freq = -1 and minr = 0; for those lines freq and minr must match the jar's dump of
    // the same command line without -p.
    auto all = cmdlines;
    for (const auto& a : cmdlines) {
        std::vector<std::string> no_p;
        for (const auto& x : a) {
            if (x != "-p") no_p.push_back(x);
        }
        all.push_back(std::move(no_p));
    }
    const auto java = java_dumps(all);
    if (java.size() != all.size())
        die(std::format("ConfigDumpMain wrote {} dumps for {} command lines", java.size(), all.size()));
    const auto n = cmdlines.size();
    std::vector<std::string> got(n);
    parallel_for(n, default_jobs(), [&](std::size_t i) { got[i] = cpp_dump(binary, cmdlines[i]); });
    std::size_t failed = 0;
    for (std::size_t i = 0; i < n; ++i) {
        auto want = normalise_java(java[i]);
        if (contains(cmdlines[i], "-p") && !want.starts_with("ERROR")) {
            std::map<std::string, std::string> keep;
            for (const auto& l : lines(normalise_java(java[n + i]))) keep[field(l)] = l;
            std::string w;
            for (const auto& l : lines(want)) {
                const auto f = field(l);
                w += ((f == "freq" || f == "minr") && keep.contains(f) ? keep[f] : l) + "\n";
            }
            want = w;
        }
        if (got[i] != want) {
            ++failed;
            std::println("DIFF   {}", join(cmdlines[i], " "));
            const auto wl = lines(want), gl = lines(got[i]);
            for (const auto& w : wl) {
                if (!contains(gl, w)) std::println("   java: {}", w);
            }
            for (const auto& g : gl) {
                if (!contains(wl, g)) std::println("    c++: {}", g);
            }
        } else if (verbose) {
            std::println("OK     {}", join(cmdlines[i], " "));
        }
    }
    std::println("{}/{} command lines parse identically", n - failed, n);
    return failed ? 1 : 0;
}

int cmd_inputs(Args& args) {
    const auto binary = fs::absolute(args.required("--binary")).string();
    const fs::path cases = args.value_or("--cases", golden_cases().string());
    const fs::path bams = args.value_or("--bams", golden_bams().string());
    const fs::path fastas = args.value_or("--fastas", golden_fastas().string());
    const auto match = args.values("-k");
    const bool no_golden = args.flag("--no-golden");
    const int jobs = args.int_or("-j", default_jobs());
    const int timeout = args.int_or("--timeout", 3600);
    const auto heap = args.value_or("--java-heap", "8g");
    const int diff = args.int_or("--diff", 5);
    const bool verbose = args.flag("-v") || args.flag("--verbose");
    const auto extra_run = args.finish();

    struct Job {
        std::string name;
        std::vector<std::string> args;
    };
    std::vector<Job> jobs_v;
    std::vector<std::unique_ptr<TempFile>> beds;
    if (!no_golden) {
        for (const auto& f : case_files(cases)) {
            const auto name = f.filename().string();
            if (!match.empty() && std::ranges::none_of(match, [&](const std::string& m) { return name.contains(m); }))
                continue;
            const auto c = parse_case(f);
            beds.push_back(std::make_unique<TempFile>(".bed", case_bed_line(c, true)));
            jobs_v.push_back({name, words(beds.back()->path() + " -b " + case_bams(c, bams) + " " + c.args +
                                          " -z -c 1 -S 2 -E 3 -g 4 -G " + case_fasta_csv(c, fastas))});
        }
        // Hand-made BAM covering every RecordPreprocessor/SamView branch (tests/data/filters.sam).
        const auto data = (root() / "tests/data").string();
        for (const char* extra : {"", "-t", "-t -F 0x500", "-t -F 0", "-Q 10 -t"}) {
            auto n = std::string("tests/data/filters.bam ") + extra;
            while (n.ends_with(' ')) n.pop_back();
            jobs_v.push_back({n, words(data + "/filters.bed -b " + data + "/filters.bam -G " + data +
                                       "/filters.fa.csv -c 1 -S 2 -E 3 -g 4 " + extra)});
        }
    }
    if (!extra_run.empty()) jobs_v.push_back({"--run " + join(extra_run, " "), extra_run});

    struct Res {
        std::string status;
        std::vector<std::string> diffs;
        long reads = 0;
    };
    std::vector<Res> res(jobs_v.size());
    parallel_for(jobs_v.size(), jobs, [&](std::size_t i) {
        std::vector<std::string> jc = {"java",
                                       "-XX:-UsePerfData",
                                       "-Xlog:disable",
                                       "-Xmx" + heap,
                                       "-cp",
                                       jar_classpath(),
                                       (root() / "tests/oracle/InputsDumpMain.java").string()};
        jc.insert(jc.end(), jobs_v[i].args.begin(), jobs_v[i].args.end());
        std::vector<std::string> cc = {binary, "--dev-dump-inputs"};
        cc.insert(cc.end(), jobs_v[i].args.begin(), jobs_v[i].args.end());
        const auto j = run(jc, timeout_opt(timeout));
        const auto c = run(cc, timeout_opt(timeout));
        if (j.rc) {
            res[i] = {"JAVA-ERROR", {last_line(j.err)}, 0};
            return;
        }
        if (c.rc) {
            res[i] = {"CPP-ERROR", {std::string(strip(c.err))}, 0};
            return;
        }
        const auto jl = lines(j.out), cl = lines(c.out);
        Res r;
        for (std::size_t k = 0; k < std::max(jl.size(), cl.size()); ++k) {
            const auto x = k < jl.size() ? jl[k] : "<missing>";
            const auto y = k < cl.size() ? cl[k] : "<missing>";
            if (x != y) r.diffs.push_back(std::format("line {}:\n      java: {}\n      c++:  {}", k + 1, x, y));
        }
        for (const auto& l : jl) {
            const auto f = split(l, '\t');
            if (l.starts_with("READS\t") && f.size() > 3 && !f[3].empty() && std::ranges::all_of(f[3], ::isdigit))
                r.reads += std::stol(f[3]);
        }
        r.status = r.diffs.empty() ? "OK" : "DIFF";
        res[i] = std::move(r);
    });
    std::size_t ok = 0;
    long total_reads = 0;
    for (std::size_t i = 0; i < res.size(); ++i) {
        total_reads += res[i].reads;
        if (res[i].status == "OK") {
            ++ok;
            if (verbose) std::println("OK         {}  ({} reads)", jobs_v[i].name, res[i].reads);
            continue;
        }
        std::println("{:<10} {}", res[i].status, jobs_v[i].name);
        for (std::size_t k = 0; k < res[i].diffs.size() && k < static_cast<std::size_t>(diff); ++k)
            std::println("    {}", res[i].diffs[k]);
        if (res[i].diffs.size() > static_cast<std::size_t>(diff))
            std::println("    ... {} more", res[i].diffs.size() - static_cast<std::size_t>(diff));
    }
    std::println("\n{}/{} runs identical; {} reads returned by RecordPreprocessor (Java)", ok, res.size(), total_reads);
    return ok == res.size() ? 0 : 1;
}

}  // namespace testkit
