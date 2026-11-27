#pragma once

// senri-testkit: the test and validation harness (regression cases, thread and -p invariants,
// truth-set benchmark, release notes). Test tooling only; not installed.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace testkit {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------------------------
// Repository layout
// ---------------------------------------------------------------------------------------------
fs::path root();  // source tree (compiled in; SENRI_SOURCE_DIR env overrides)
inline fs::path regression_dir() { return root() / "tests/regression"; }
inline fs::path cases_dir() { return regression_dir() / "cases"; }  // case spec + expected TSV
inline fs::path vcf_dir() { return regression_dir() / "vcf"; }      // expected VCF per case
inline fs::path bams_dir() { return regression_dir() / "bams"; }
inline fs::path fastas_dir() { return regression_dir() / "fastas"; }  // reference slices (CSV)
inline fs::path synthetic_dir() { return regression_dir() / "synthetic"; }

// ---------------------------------------------------------------------------------------------
// Strings and files
// ---------------------------------------------------------------------------------------------
std::vector<std::string> split(std::string_view s, char sep);  // like str.split(sep)
std::vector<std::string> words(std::string_view s);            // like str.split()
std::vector<std::string> lines(std::string_view s);            // like str.splitlines()
std::vector<std::string> lines_keepends(std::string_view s);   // like splitlines(keepends=True)
bool is_blank(std::string_view s);                             // not s.strip()
std::string_view strip(std::string_view s);
std::string rstrip_newline(std::string_view s);
std::string join(const std::vector<std::string>& v, std::string_view sep);
std::vector<std::string> sorted_rows(std::string_view text);  // sorted non-blank lines
bool contains(const std::vector<std::string>& v, std::string_view x);
std::string repr(std::string_view s);                 // 'text' (for messages)
std::string last_line(std::string_view s);            // last line of strip(s), or ""
std::string tail(std::string_view s, std::size_t n);  // last n chars of strip(s)
std::optional<double> as_number(std::string_view s);  // float(s), or nothing

std::string read_file(const fs::path& p);
void write_file(const fs::path& p, std::string_view text);
std::optional<std::string> which(std::string_view program);

// Unlinked on destruction.
class TempFile {
public:
    explicit TempFile(std::string_view suffix = "", std::string_view content = "");
    ~TempFile();
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;
    [[nodiscard]] const std::string& path() const { return path_; }
    void write(std::string_view content) const;

private:
    std::string path_;
};

// Removed (recursively) on destruction.
class TempDir {
public:
    TempDir();
    ~TempDir();
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    [[nodiscard]] const fs::path& path() const { return path_; }

private:
    fs::path path_;
};

// ---------------------------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------------------------
struct RunOptions {
    std::optional<std::string> input;  // stdin text (none: /dev/null)
    int timeout_s = 600;
    std::string stdout_path;  // redirect stdout to this file instead of capturing it
    std::string stderr_path;
};

inline RunOptions timeout_opt(int timeout_s) {
    RunOptions o;
    o.timeout_s = timeout_s;
    return o;
}
inline RunOptions input_opt(std::string input, int timeout_s = 600) {
    RunOptions o;
    o.input = std::move(input);
    o.timeout_s = timeout_s;
    return o;
}
inline RunOptions files_opt(std::string out, std::string err, int timeout_s) {
    RunOptions o;
    o.stdout_path = std::move(out);
    o.stderr_path = std::move(err);
    o.timeout_s = timeout_s;
    return o;
}

struct RunResult {
    int rc = 0;  // exit status, or -signal
    bool timed_out = false;
    std::string out, err;
    double seconds = 0;
};

RunResult run(const std::vector<std::string>& argv, const RunOptions& opt = {});

// Calls f(i) for i in [0, n) on `jobs` threads.
template <class F>
void parallel_for(std::size_t n, int jobs, F&& f) {
    std::atomic<std::size_t> next{0};
    auto work = [&] {
        for (std::size_t i; (i = next.fetch_add(1)) < n;) f(i);
    };
    std::vector<std::jthread> pool;
    const auto k = std::min<std::size_t>(n, static_cast<std::size_t>(std::max(1, jobs)));
    for (std::size_t t = 1; t < k; ++t) pool.emplace_back(work);
    work();
}

int default_jobs();

// ---------------------------------------------------------------------------------------------
// Command-line options of a subcommand
// ---------------------------------------------------------------------------------------------
class Args {
public:
    Args(int argc, char** argv);  // argv[0] is the subcommand
    // Last value of an option ("--x v", "--x=v", and "-xV" for one-letter options).
    std::optional<std::string> value(std::string_view name);
    std::string value_or(std::string_view name, std::string fallback);
    int int_or(std::string_view name, int fallback);
    std::vector<std::string> values(std::string_view name);  // every occurrence
    bool flag(std::string_view name);
    std::string required(std::string_view name);
    // Arguments after "--" (when allow_rest) and positional arguments; errors on unknown options.
    std::vector<std::string> finish(bool allow_positional = false);

private:
    std::vector<std::string> tok_;
    std::vector<bool> used_;
    std::size_t end_;  // index of "--" or size
    std::string cmd_;
};

[[noreturn]] void die(std::string_view msg);

// ---------------------------------------------------------------------------------------------
// Regression cases (tests/regression/cases/*.txt): line 1 is the case spec
// (mode,fasta,bam,chr,start,end[,amp_start,amp_end][,args]), the rest is the expected TSV
// ---------------------------------------------------------------------------------------------
struct Case {
    std::string name;  // file name
    std::string mode, fasta, bam, chr, start, end;
    std::optional<std::pair<std::string, std::string>> amp;
    std::string args;
    std::string expected;  // the lines after the header
};

Case parse_case(const fs::path& path);
std::vector<fs::path> case_files(const fs::path& dir, bool txt_only = true);
std::string case_bams(const Case& c, const fs::path& bams_dir);  // "a.bam|b.bam" with directory
std::string case_fasta_csv(const Case& c, const fs::path& fastas_dir);
std::string case_bed_line(const Case& c, bool with_amp, std::string_view name = "testbed");
bool case_has_arg(const Case& c, std::string_view arg);

// ---------------------------------------------------------------------------------------------
// Subcommands
// ---------------------------------------------------------------------------------------------
int cmd_regress(Args& a);
int cmd_threads(Args& a);
int cmd_reject(Args& a);
int cmd_bench(Args& a);
int cmd_release_notes(Args& a);

}  // namespace testkit
