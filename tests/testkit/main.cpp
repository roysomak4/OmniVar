#include <csignal>
#include <print>
#include <string_view>

#include "testkit.hpp"

namespace {

constexpr std::string_view usage = R"(usage: senri-testkit <command> [options]

Regression and invariants (the ctest gates):
  regress           senri output vs its checked-in output (tests/regression)
                    --binary B [--part tsv,vcf,synthetic] [--mode Simple|Somatic]... [-k substr]...
                    [--update] [--diff N] [-v] [-j N]
  threads           -th N output equals -th 1   --binary B [-th 4] [--chunk 500] [-j N]
  reject-invariant  -p never changes the calls  --binary B [--extra-args=ARGS] [-j N]

Benchmarks and releases:
  bench             truth-set benchmark (aardvark), optional call snapshot
                    (--binary B -b BAM|"T|N" [-N name] | --vcf FILE) --reference FA --bed BED --truth VCF
                    [--snapshot FILE [--update-snapshot]] [--min-recall R] [--somatic] [--min-af X]
                    [--out DIR] [-th 16] [-- caller args]
  release-notes     draft release notes  --tag vX.Y.Z [--previous TAG] [--checksums SHA256SUMS]
                    [--bench bench.txt] [--build TEXT] [--repo owner/name] [--highlights FILE]

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
        {"regress", testkit::cmd_regress},
        {"threads", testkit::cmd_threads},
        {"reject-invariant", testkit::cmd_reject},
        {"bench", testkit::cmd_bench},
        {"release-notes", testkit::cmd_release_notes},
    };
    for (const auto& [name, f] : table) {
        if (name == cmd) return f(a);
    }
    std::print(stderr, "senri-testkit: unknown command {}\n\n{}", cmd, usage);
    return 2;
}
