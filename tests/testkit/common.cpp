#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <poll.h>
#include <print>
#include <signal.h>
#include <spawn.h>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>

#include "testkit.hpp"

extern char** environ;

namespace testkit {

fs::path root() {
    if (const char* env = std::getenv("SENRI_SOURCE_DIR"); env && *env) return env;
    return SENRI_SOURCE_DIR;
}

[[noreturn]] void die(std::string_view msg) {
    std::println(stderr, "senri-testkit: {}", msg);
    std::exit(2);
}

// ---------------------------------------------------------------------------------------------
std::vector<std::string> split(std::string_view s, char sep) {
    std::vector<std::string> out;
    std::size_t b = 0;
    for (;;) {
        const auto e = s.find(sep, b);
        if (e == std::string_view::npos) {
            out.emplace_back(s.substr(b));
            return out;
        }
        out.emplace_back(s.substr(b, e - b));
        b = e + 1;
    }
}

static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f'; }

std::vector<std::string> words(std::string_view s) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && is_space(s[i])) ++i;
        const auto b = i;
        while (i < s.size() && !is_space(s[i])) ++i;
        if (i > b) out.emplace_back(s.substr(b, i - b));
    }
    return out;
}

std::vector<std::string> lines_keepends(std::string_view s) {
    std::vector<std::string> out;
    std::size_t b = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\n' || s[i] == '\r') {
            auto e = i + 1;
            if (s[i] == '\r' && e < s.size() && s[e] == '\n') ++e;
            out.emplace_back(s.substr(b, e - b));
            b = e;
            i = e - 1;
        }
    }
    if (b < s.size()) out.emplace_back(s.substr(b));
    return out;
}

std::vector<std::string> lines(std::string_view s) {
    auto v = lines_keepends(s);
    for (auto& l : v) {
        while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
    }
    return v;
}

bool is_blank(std::string_view s) { return std::ranges::all_of(s, is_space); }

std::string_view strip(std::string_view s) {
    while (!s.empty() && is_space(s.front())) s.remove_prefix(1);
    while (!s.empty() && is_space(s.back())) s.remove_suffix(1);
    return s;
}

std::string rstrip_newline(std::string_view s) {
    while (!s.empty() && s.back() == '\n') s.remove_suffix(1);
    return std::string(s);
}

std::string join(const std::vector<std::string>& v, std::string_view sep) {
    std::string out;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) out += sep;
        out += v[i];
    }
    return out;
}

std::vector<std::string> sorted_rows(std::string_view text) {
    std::vector<std::string> out;
    for (auto& l : lines(text)) {
        if (!is_blank(l)) out.push_back(std::move(l));
    }
    std::ranges::sort(out);
    return out;
}

bool contains(const std::vector<std::string>& v, std::string_view x) { return std::ranges::find(v, x) != v.end(); }

std::string repr(std::string_view s) {
    std::string out = "'";
    for (const char c : s) {
        switch (c) {
            case '\t': out += "\\t"; break;
            case '\n': out += "\\n"; break;
            case '\'': out += "\\'"; break;
            case '\\': out += "\\\\"; break;
            default: out += c;
        }
    }
    return out + "'";
}

std::string last_line(std::string_view s) {
    const auto v = lines(strip(s));
    return v.empty() ? std::string() : v.back();
}

std::string tail(std::string_view s, std::size_t n) {
    s = strip(s);
    return std::string(s.size() > n ? s.substr(s.size() - n) : s);
}

std::optional<double> as_number(std::string_view s) {
    s = strip(s);
    if (s.empty() || s.find_first_of("xXpP(") != std::string_view::npos) return std::nullopt;
    std::string t;
    for (std::size_t i = 0; i < s.size(); ++i) {  // Python allows '_' between digits
        if (s[i] == '_') {
            if (i == 0 || i + 1 == s.size() || !std::isdigit(static_cast<unsigned char>(s[i - 1])) ||
                !std::isdigit(static_cast<unsigned char>(s[i + 1])))
                return std::nullopt;
            continue;
        }
        t += s[i];
    }
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(t.c_str(), &end);
    if (end != t.c_str() + t.size()) return std::nullopt;
    return v;
}

std::string read_file(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) die("cannot read " + p.string());
    std::ostringstream ss;
    ss << in.rdbuf();
    return std::move(ss).str();
}

void write_file(const fs::path& p, std::string_view text) {
    std::ofstream out(p, std::ios::binary | std::ios::trunc);
    if (!out) die("cannot write " + p.string());
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
}

std::optional<std::string> which(std::string_view program) {
    const char* path = std::getenv("PATH");
    if (!path) return std::nullopt;
    for (const auto& dir : split(path, ':')) {
        const auto p = fs::path(dir.empty() ? "." : dir) / program;
        if (::access(p.c_str(), X_OK) == 0 && !fs::is_directory(p)) return p.string();
    }
    return std::nullopt;
}

static std::string temp_base() {
    const char* t = std::getenv("TMPDIR");
    return (t && *t) ? t : "/tmp";
}

TempFile::TempFile(std::string_view suffix, std::string_view content) {
    std::string tmpl = temp_base() + "/senri-testkit-XXXXXX" + std::string(suffix);
    const int fd = ::mkstemps(tmpl.data(), static_cast<int>(suffix.size()));
    if (fd < 0) die("mkstemps failed: " + std::string(std::strerror(errno)));
    ::close(fd);
    path_ = tmpl;
    if (!content.empty()) write(content);
}

TempFile::~TempFile() { ::unlink(path_.c_str()); }

void TempFile::write(std::string_view content) const { write_file(path_, content); }

TempDir::TempDir() {
    std::string tmpl = temp_base() + "/senri-testkit-XXXXXX";
    if (!::mkdtemp(tmpl.data())) die("mkdtemp failed: " + std::string(std::strerror(errno)));
    path_ = tmpl;
}

TempDir::~TempDir() {
    std::error_code ec;
    fs::remove_all(path_, ec);
}

// ---------------------------------------------------------------------------------------------
RunResult run(const std::vector<std::string>& argv, const RunOptions& opt) {
    RunResult r;
    const auto t0 = std::chrono::steady_clock::now();
    int in_pipe[2] = {-1, -1}, out_pipe[2] = {-1, -1}, err_pipe[2] = {-1, -1};
    auto mkpipe = [](int p[2]) {
        if (::pipe2(p, O_CLOEXEC) != 0) die("pipe2 failed");
    };
    if (opt.input) mkpipe(in_pipe);
    if (opt.stdout_path.empty()) mkpipe(out_pipe);
    if (opt.stderr_path.empty()) mkpipe(err_pipe);

    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    if (opt.input)
        posix_spawn_file_actions_adddup2(&fa, in_pipe[0], 0);
    else
        posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
    if (opt.stdout_path.empty())
        posix_spawn_file_actions_adddup2(&fa, out_pipe[1], 1);
    else
        posix_spawn_file_actions_addopen(&fa, 1, opt.stdout_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (opt.stderr_path.empty())
        posix_spawn_file_actions_adddup2(&fa, err_pipe[1], 2);
    else
        posix_spawn_file_actions_addopen(&fa, 2, opt.stderr_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);

    std::vector<char*> cargv;
    for (const auto& s : argv) cargv.push_back(const_cast<char*>(s.c_str()));
    cargv.push_back(nullptr);
    pid_t pid = 0;
    const int rc = ::posix_spawnp(&pid, cargv[0], &fa, nullptr, cargv.data(), environ);
    posix_spawn_file_actions_destroy(&fa);
    for (const int fd : {in_pipe[0], out_pipe[1], err_pipe[1]}) {
        if (fd >= 0) ::close(fd);
    }
    if (rc != 0) {
        for (const int fd : {in_pipe[1], out_pipe[0], err_pipe[0]}) {
            if (fd >= 0) ::close(fd);
        }
        r.rc = 127;
        r.err = "cannot run " + argv[0] + ": " + std::strerror(rc);
        return r;
    }

    std::size_t written = 0;
    const std::string_view input = opt.input ? std::string_view(*opt.input) : std::string_view();
    int in_fd = in_pipe[1];
    if (in_fd >= 0) {
        ::fcntl(in_fd, F_SETFL, O_NONBLOCK);
        if (input.empty()) {
            ::close(in_fd);
            in_fd = -1;
        }
    }
    int fds[2] = {out_pipe[0], err_pipe[0]};
    std::string* sinks[2] = {&r.out, &r.err};
    const auto deadline = t0 + std::chrono::seconds(opt.timeout_s);
    char buf[65536];
    while (in_fd >= 0 || fds[0] >= 0 || fds[1] >= 0) {
        std::vector<pollfd> p;
        if (in_fd >= 0) p.push_back({in_fd, POLLOUT, 0});
        for (const int fd : fds) {
            if (fd >= 0) p.push_back({fd, POLLIN, 0});
        }
        const auto left =
            std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::steady_clock::now());
        if (left.count() <= 0) {
            ::kill(pid, SIGKILL);
            r.timed_out = true;
            break;
        }
        const int n = ::poll(p.data(), p.size(), static_cast<int>(std::min<long long>(left.count(), 1000)));
        if (n < 0 && errno == EINTR) continue;
        for (const auto& q : p) {
            if (!q.revents) continue;
            if (q.fd == in_fd) {
                const auto w =
                    ::write(in_fd, input.data() + written, std::min<std::size_t>(input.size() - written, 65536));
                if (w > 0) written += static_cast<std::size_t>(w);
                if (w < 0 && errno != EAGAIN) written = input.size();  // EPIPE: the child stopped reading
                if (written >= input.size()) {
                    ::close(in_fd);
                    in_fd = -1;
                }
                continue;
            }
            for (int k = 0; k < 2; ++k) {
                if (q.fd != fds[k]) continue;
                const auto got = ::read(fds[k], buf, sizeof buf);
                if (got > 0) {
                    sinks[k]->append(buf, static_cast<std::size_t>(got));
                } else if (got == 0 || errno != EAGAIN) {
                    ::close(fds[k]);
                    fds[k] = -1;
                }
            }
        }
    }
    for (const int fd : {in_fd, fds[0], fds[1]}) {
        if (fd >= 0) ::close(fd);
    }
    int status = 0;
    while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    r.rc = WIFEXITED(status) ? WEXITSTATUS(status) : (WIFSIGNALED(status) ? -WTERMSIG(status) : -1);
    r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return r;
}

int default_jobs() {
    const auto n = std::thread::hardware_concurrency();
    return n ? static_cast<int>(n) : 4;
}

// ---------------------------------------------------------------------------------------------
Args::Args(int argc, char** argv) : cmd_(argv[0]) {
    for (int i = 1; i < argc; ++i) tok_.emplace_back(argv[i]);
    used_.assign(tok_.size(), false);
    end_ = static_cast<std::size_t>(std::ranges::find(tok_, "--") - tok_.begin());
}

std::vector<std::string> Args::values(std::string_view name) {
    std::vector<std::string> out;
    const bool short_opt = name.size() == 2 && name[0] == '-' && name[1] != '-';
    for (std::size_t i = 0; i < end_; ++i) {
        if (used_[i]) continue;
        const std::string& t = tok_[i];
        if (t == name) {
            if (i + 1 >= end_) die(cmd_ + ": " + std::string(name) + " needs a value");
            used_[i] = used_[i + 1] = true;
            out.push_back(tok_[i + 1]);
            ++i;
        } else if (t.starts_with(name) && t.size() > name.size() && t[name.size()] == '=' && !short_opt) {
            used_[i] = true;
            out.push_back(t.substr(name.size() + 1));
        } else if (short_opt && t.starts_with(name) && t.size() > 2) {
            used_[i] = true;
            out.push_back(t.substr(2));
        }
    }
    return out;
}

std::optional<std::string> Args::value(std::string_view name) {
    auto v = values(name);
    if (v.empty()) return std::nullopt;
    return v.back();
}

std::string Args::value_or(std::string_view name, std::string fallback) {
    return value(name).value_or(std::move(fallback));
}

int Args::int_or(std::string_view name, int fallback) {
    const auto v = value(name);
    if (!v) return fallback;
    try {
        std::size_t pos = 0;
        const int n = std::stoi(*v, &pos);
        if (pos != v->size()) throw std::invalid_argument("trailing");
        return n;
    } catch (const std::exception&) {
        die(cmd_ + ": " + std::string(name) + " expects an integer, got " + *v);
    }
}

bool Args::flag(std::string_view name) {
    bool seen = false;
    for (std::size_t i = 0; i < end_; ++i) {
        if (!used_[i] && tok_[i] == name) {
            used_[i] = true;
            seen = true;
        }
    }
    return seen;
}

std::string Args::required(std::string_view name) {
    auto v = value(name);
    if (!v) die(cmd_ + ": " + std::string(name) + " is required");
    return *v;
}

std::vector<std::string> Args::finish(bool allow_positional) {
    std::vector<std::string> rest;
    for (std::size_t i = 0; i < end_; ++i) {
        if (used_[i]) continue;
        if (!allow_positional || (tok_[i].starts_with("-") && tok_[i].size() > 1))
            die(cmd_ + ": unknown argument " + tok_[i]);
        rest.push_back(tok_[i]);
    }
    for (std::size_t i = end_ + 1; i < tok_.size(); ++i) rest.push_back(tok_[i]);
    return rest;
}

// ---------------------------------------------------------------------------------------------
Case parse_case(const fs::path& path) {
    const auto text = read_file(path);
    const auto nl = text.find('\n');
    const std::string header = text.substr(0, nl);
    Case c;
    c.name = path.filename().string();
    c.expected = nl == std::string::npos ? std::string() : text.substr(nl + 1);
    const auto f = split(header, ',');
    if (f.size() < 6) die("bad case header in " + path.string());
    c.mode = f[0];
    c.fasta = f[1];
    c.bam = f[2];
    c.chr = f[3];
    c.start = f[4];
    c.end = f[5];
    if (f.size() < 8) {
        c.args = f.size() == 7 ? f[6] : "";
    } else {
        c.amp = std::pair{f[6], f[7]};
        c.args = f.size() == 9 ? f[8] : "";
    }
    return c;
}

std::vector<fs::path> case_files(const fs::path& dir, bool txt_only) {
    std::vector<fs::path> out;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (!txt_only || e.path().extension() == ".txt") out.push_back(e.path());
    }
    std::ranges::sort(out, {}, [](const fs::path& p) { return p.filename().string(); });
    return out;
}

std::string case_bams(const Case& c, const fs::path& bams_dir) {
    std::vector<std::string> v;
    for (const auto& b : split(c.bam, '|')) v.push_back((bams_dir / b).string());
    return join(v, "|");
}

std::string case_fasta_csv(const Case& c, const fs::path& fastas_dir) {
    return (fastas_dir / (c.fasta + ".csv")).string();
}

std::string case_bed_line(const Case& c, bool with_amp, std::string_view name) {
    std::string line = c.chr + "\t" + c.start + "\t" + c.end + "\t" + std::string(name);
    if (with_amp && c.amp) line += "\t0\t.\t" + c.amp->first + "\t" + c.amp->second;
    return line;
}

bool case_has_arg(const Case& c, std::string_view arg) { return contains(words(c.args), arg); }

}  // namespace testkit
