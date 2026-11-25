#include <csignal>
#include <print>
#include <string_view>

#include "testkit.hpp"

namespace {

constexpr std::string_view usage = R"(usage: omnivar-testkit <command> [options]

Golden cases (no Java needed; the ctest gates):
  golden            omnivar TSV vs the golden / stage expectations
                    --binary B [--stage full|U|k0_U] [--mode M]... [-k substr]... [--diff N] [-v]
                    [--expected DIR] [--extra-args=ARGS] [--out DIR] [--update-rejected] [-j N]
  vcf               omnivar VCF vs the R + Perl chain (checked-in expectations)
                    --binary B [--stage k0_U|U|full] [--mode M]... [--variant V] [--synthetic]
                    [-k substr] [--diff N] [-j N];  --generate [--stage S] (needs java, Rscript, perl)
  threads           -th N output equals -th 1   --binary B [-th 4] [--chunk 500] [-j N]
  reject-invariant  -p never changes the calls  --binary B [--extra-args=ARGS] [-j N]

Reference implementation (needs Java 21 and tests/oracle/lib):
  oracle            run the 1.8.3 jar on the golden cases  [--check] [--out DIR] [--extra-args=ARGS] [-j N]
  nopileup          regenerate tests/golden/expected_nopileup  [-j N]
  cli-parity        command-line parser vs the jar  --binary B [--golden] [-v]
  inputs-parity     samples, regions, reference, reads vs the jar  --binary B [-j N] [-- caller args]
  panel             panel-scale parity on user data  --binary B --out DIR --reference FA --bed BED
                    -b BAM|"T|N" [-N name|"T|N"] [--vcf-args=..] [--perl-args=..] [-- caller args]

Generators and benchmarks:
  synthetic         --single | --paired: synthetic TSV rows for the VCF checks (stdout)
  bench             truth-set benchmark (aardvark), optional call snapshot
                    (--binary B -b BAM|"T|N" [-N name] | --vcf FILE) --reference FA --bed BED --truth VCF
                    [--snapshot FILE [--update-snapshot]] [--min-recall R] [--somatic] [--min-af X]
                    [--out DIR] [-th 16] [-- caller args]
  release-notes     draft release notes  --tag vX.Y.Z [--previous TAG] [--checksums SHA256SUMS]
                    [--bench bench.txt] [--build TEXT] [--repo owner/name]

Option values may follow with a space or "=" (pass values starting with "-" as --extra-args="-U").
)";

}  // namespace

int main(int argc, char** argv) {
    (void)std::signal(SIGPIPE, SIG_IGN);
    if (argc < 2 || std::string_view(argv[1]) == "-h" || std::string_view(argv[1]) == "--help") {
        std::print("{}", usage);
        return argc < 2 ? 2 : 0;
    }
    const std::string_view cmd = argv[1];
    testkit::Args a(argc - 1, argv + 1);
    using F = int (*)(testkit::Args&);
    const std::pair<std::string_view, F> table[] = {
        {"golden", testkit::cmd_golden},   {"vcf", testkit::cmd_vcf},
        {"threads", testkit::cmd_threads}, {"reject-invariant", testkit::cmd_reject},
        {"oracle", testkit::cmd_oracle},   {"nopileup", testkit::cmd_nopileup},
        {"cli-parity", testkit::cmd_cli},  {"inputs-parity", testkit::cmd_inputs},
        {"panel", testkit::cmd_panel},     {"synthetic", testkit::cmd_synthetic},
        {"bench", testkit::cmd_bench},     {"release-notes", testkit::cmd_release_notes},
    };
    for (const auto& [name, f] : table) {
        if (name == cmd) return f(a);
    }
    std::print(stderr, "omnivar-testkit: unknown command {}\n\n{}", cmd, usage);
    return 2;
}
