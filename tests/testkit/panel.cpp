// panel (M8): omnivar against the 1.8.3 jar on real data. Runs, in parallel:
//     jar      java -cp tests/oracle/lib/* com.astrazeneca.vardict.Main -th N <args>   (TSV)
//     omnivar  --tsv <args>                                                          (TSV)
//     omnivar  <args> [vcf args]                                                     (VCF)
// then the annotation chain on the jar's TSV (single: teststrandbias.R | var2vcf_valid.pl -A;
// paired: testsomatic.R | var2vcf_paired.pl -A -N "T|N") and compares:
//     TSV  the same multiset of non-empty rows
//     VCF  records as in `vcf` (vcfcmp.hpp); bcftools view must parse the omnivar VCF.
// Data paths come from the command line (user data is never committed). Outputs and summary.txt
// go to --out. Exit 0 only when both comparisons are clean.

#include <future>
#include <print>
#include <set>

#include "testkit.hpp"
#include "vcfcmp.hpp"

namespace testkit {

int cmd_panel(Args& args) {
    const auto binary = fs::absolute(args.required("--binary")).string();
    const fs::path out = args.required("--out");
    const auto reference = args.required("--reference");
    const auto bed = args.required("--bed");
    const auto bam = args.required("-b");
    const auto name = args.value("-N");
    const int java_threads = args.int_or("--java-threads", 16);
    const auto vcf_args = args.value_or("--vcf-args", "");
    auto perl_args = args.value_or("--perl-args", "-A");
    const int diff = args.int_or("--diff", 10);
    const auto extra = args.finish();
    const bool paired = bam.contains('|');
    fs::create_directories(out);
    auto o = [&](const char* f) { return (out / f).string(); };

    std::vector<std::string> common = {"-G", reference, "-b", bam};
    if (name) {
        common.emplace_back("-N");
        common.push_back(*name);
    }
    for (auto& w : words("-c 1 -S 2 -E 3 -g 4")) common.push_back(std::move(w));
    common.insert(common.end(), extra.begin(), extra.end());
    common.push_back(bed);
    std::vector<std::string> jar = {
        "java", "-Xmx16g", "-cp", jar_classpath(), "com.astrazeneca.vardict.Main", "-th", std::to_string(java_threads)};
    jar.insert(jar.end(), common.begin(), common.end());
    std::vector<std::string> tsv = {binary, "--tsv"};
    tsv.insert(tsv.end(), common.begin(), common.end());
    std::vector<std::string> vcf = {binary};
    vcf.insert(vcf.end(), common.begin(), common.end());
    for (auto& w : words(vcf_args)) vcf.push_back(std::move(w));

    constexpr int day = 86400;
    auto fj = std::async(std::launch::async, [&] { return run(jar, files_opt(o("java.tsv"), o("java.err"), day)); });
    auto ft = std::async(std::launch::async, [&] { return run(tsv, files_opt(o("cpp.tsv"), o("cpp.err"), day)); });
    auto fv = std::async(std::launch::async, [&] { return run(vcf, files_opt(o("cpp.vcf"), o("cpp_vcf.err"), day)); });
    const auto rj = fj.get(), rt = ft.get(), rv = fv.get();

    std::vector<std::string> report = {
        std::format("mode: {}  args: {}  vcf: {}", paired ? "paired" : "single", join(extra, " "), vcf_args),
        std::format("time: jar {:.0f} s (-th {}), C++ TSV {:.0f} s, C++ VCF {:.0f} s (1 thread)", rj.seconds,
                    java_threads, rt.seconds, rv.seconds)};
    const bool ok = rj.rc == 0 && rt.rc == 0 && rv.rc == 0;
    if (!ok) report.push_back(std::format("exit codes: jar {}, C++ TSV {}, C++ VCF {}", rj.rc, rt.rc, rv.rc));

    const auto want = sorted_rows(read_file(o("java.tsv"))), got = sorted_rows(read_file(o("cpp.tsv")));
    const std::set<std::string> ws(want.begin(), want.end()), gs(got.begin(), got.end());
    std::vector<std::string> only_j, only_c;
    for (const auto& r : want) {
        if (!gs.contains(r)) only_j.push_back(r);
    }
    for (const auto& r : got) {
        if (!ws.contains(r)) only_c.push_back(r);
    }
    const bool tsv_ok = want == got;
    report.push_back(std::format(
        "TSV: jar {} rows, C++ {} rows, {}", want.size(), got.size(),
        tsv_ok ? "identical" : std::format("{} only in jar, {} only in C++", only_j.size(), only_c.size())));
    for (std::size_t i = 0; i < only_j.size() && i < static_cast<std::size_t>(diff); ++i)
        report.push_back("  jar only: " + only_j[i].substr(0, 200));
    for (std::size_t i = 0; i < only_c.size() && i < static_cast<std::size_t>(diff); ++i)
        report.push_back("  C++ only: " + only_c[i].substr(0, 200));

    if (paired && name && name->contains('|') && !contains(words(perl_args), "-N")) perl_args += " -N " + *name;
    if (name && name->contains(' '))
        die("sample names with spaces are not supported (the chain splits options on spaces)");
    std::string expected;
    try {
        expected = chain(read_file(o("java.tsv")), 7200, perl_args, paired);
    } catch (const std::exception& e) {
        die(e.what());
    }
    write_file(o("chain.vcf"), expected);
    const auto got_vcf = read_file(o("cpp.vcf"));
    const auto bad = bcftools_check(got_vcf);
    const auto w = parse_records(expected), g = parse_records(got_vcf);
    std::vector<std::string> known;
    const auto problems = compare_records(w, g, known);
    const bool vcf_ok = problems.empty() && !bad;
    report.push_back(std::format("VCF: chain {} records, C++ {}, bcftools {}, {} differences{}", w.count(), g.count(),
                                 bad ? *bad : "clean", problems.size(),
                                 known.empty() ? "" : std::format(", {} known R artefacts", known.size())));
    for (std::size_t i = 0; i < problems.size() && i < static_cast<std::size_t>(diff); ++i)
        report.push_back("  " + problems[i]);
    const bool pass = ok && tsv_ok && vcf_ok;
    report.push_back(std::string("RESULT: ") + (pass ? "PASS" : "FAIL"));
    const auto text = join(report, "\n");
    std::println("{}", text);
    write_file(o("summary.txt"), text + "\n");
    return pass ? 0 : 1;
}

}  // namespace testkit
