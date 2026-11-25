// bench: score omnivar calls against a truth set with `aardvark compare` (haplotype-aware, the
// comparator GIAB recommends), and optionally check the call set against a snapshot.
//
//   omnivar-testkit bench --binary build/omnivar --reference ref.fa --bed targets.bed -b s.bam
//       [-N name] --truth truth.vcf.gz [--snapshot FILE [--update-snapshot]] [--min-recall 1]
//       [--somatic] [--min-af X] [--out DIR] [-th 16] [-- extra caller options]
//   omnivar-testkit bench --vcf calls.vcf --reference ref.fa --bed eval.bed --truth truth.vcf.gz ...
//
// 1. omnivar writes the VCF (default output) for the targets (-b "T|N" for a tumour/normal pair),
//    or --vcf names an existing VCF.
// 2. The query is the PASS records; with --somatic only those whose INFO STATUS contains "Somatic"
//    (paired VCF); with --min-af only those with INFO AF >= X.
// 3. The truth VCF is copied with a repaired header (trailing whitespace removed, Number=. added
//    to INFO/FORMAT lines without one); both VCFs are bgzipped and tabix-indexed (htslib), as
//    aardvark requires. The BED is reduced to its first three columns.
// 4. aardvark compare scores the query (GT metrics per variant type). Fail when recall is below
//    --min-recall (default 0: report only).
// 5. Snapshot (--snapshot): record count, counts per FILTER tag and per TYPE x PASS, an FNV-1a
//    digest of the sorted records, and the aardvark GT rows. Any change fails until the snapshot
//    is regenerated with --update-snapshot (a deliberate behaviour change). Counts and digests
//    only, no data.

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <print>

#include "testkit.hpp"

namespace testkit {

namespace {

std::string read_any(const fs::path& p) {  // plain, gzip or BGZF
    BGZF* f = bgzf_open(p.c_str(), "r");
    if (!f) die("cannot open " + p.string());
    std::string out;
    char buf[65536];
    for (ssize_t n; (n = bgzf_read(f, buf, sizeof buf)) > 0;) out.append(buf, static_cast<std::size_t>(n));
    bgzf_close(f);
    return out;
}

void write_bgzf_indexed(const fs::path& p, std::string_view text) {
    BGZF* f = bgzf_open(p.c_str(), "w");
    if (!f) die("cannot write " + p.string());
    if (bgzf_write(f, text.data(), text.size()) < 0 || bgzf_close(f) != 0) die("bgzf write failed: " + p.string());
    if (tbx_index_build(p.c_str(), 0, &tbx_conf_vcf) != 0) die("tabix indexing failed: " + p.string());
}

std::string repair_header(std::string_view vcf) {
    std::string out;
    for (auto l : lines(vcf)) {
        if (l.starts_with("##")) {
            while (!l.empty() && (l.back() == ' ' || l.back() == '\t')) l.pop_back();
            if ((l.starts_with("##INFO=<") || l.starts_with("##FORMAT=<")) && !l.contains("Number=")) {
                const auto at = l.find('<') + 1;
                const auto comma = l.find(',', at);
                if (comma != std::string::npos) l.insert(comma + 1, "Number=.,");
            }
        }
        out += l + "\n";
    }
    return out;
}

std::optional<std::string> info_value(const std::string& info, std::string_view key) {
    for (const auto& kv : split(info, ';')) {
        if (kv.size() > key.size() && kv.starts_with(key) && kv[key.size()] == '=') return kv.substr(key.size() + 1);
    }
    return std::nullopt;
}

// The records aardvark scores: PASS, optionally Somatic STATUS and a minimum AF.
std::string query_records(std::string_view vcf, bool somatic, double min_af, std::size_t& kept) {
    std::string out;
    kept = 0;
    for (const auto& l : lines(vcf)) {
        if (l.starts_with('#')) {
            out += l + "\n";
            continue;
        }
        const auto f = split(l, '\t');
        if (f.size() < 8 || f[6] != "PASS") continue;
        if (somatic && !info_value(f[7], "STATUS").value_or("").contains("Somatic")) continue;
        if (min_af > 0) {
            const auto af = as_number(split(info_value(f[7], "AF").value_or(""), ',').front());
            if (!af || *af < min_af) continue;
        }
        out += l + "\n";
        ++kept;
    }
    return out;
}

std::uint64_t fnv1a(std::uint64_t h, std::string_view s) {
    for (const unsigned char c : s) {
        h ^= c;
        h *= 0x100000001b3ULL;
    }
    return h;
}

using Snapshot = std::vector<std::pair<std::string, std::string>>;

Snapshot call_snapshot(std::string_view vcf) {
    std::vector<std::string> recs;
    std::map<std::string, long> filters, types;
    long pass = 0;
    for (auto& l : lines(vcf)) {
        if (l.empty() || l.starts_with('#')) continue;
        const auto f = split(l, '\t');
        const std::string flt = f.size() > 6 ? f[6] : "";
        const std::string type = f.size() > 7 ? info_value(f[7], "TYPE").value_or("none") : "none";
        pass += flt == "PASS";
        for (const auto& t : split(flt, ';')) ++filters[t];
        ++types[type + (flt == "PASS" ? ".PASS" : ".filtered")];
        recs.push_back(std::move(l));
    }
    std::ranges::sort(recs);
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (const auto& r : recs) h = fnv1a(fnv1a(h, r), "\n");
    Snapshot s = {{"records", std::to_string(recs.size())}, {"pass", std::to_string(pass)}};
    for (const auto& [k, v] : filters) s.emplace_back("filter." + k, std::to_string(v));
    for (const auto& [k, v] : types) s.emplace_back("type." + k, std::to_string(v));
    s.emplace_back("digest", std::format("fnv1a64:{:016x}", h));
    return s;
}

bool check_snapshot(const fs::path& snapshot, const Snapshot& snap, bool update) {
    std::string text =
        "# omnivar-testkit bench snapshot: counts and digests only (regenerate with --update-snapshot)\n";
    for (const auto& [k, v] : snap) text += k + "\t" + v + "\n";
    if (update) {
        fs::create_directories(fs::absolute(snapshot).parent_path());
        write_file(snapshot, text);
        std::println("snapshot written: {}", snapshot.string());
        return true;
    }
    if (!fs::exists(snapshot)) {
        std::println("FAIL: no snapshot {} (run with --update-snapshot)", snapshot.string());
        return false;
    }
    std::map<std::string, std::string> want, got(snap.begin(), snap.end());
    for (const auto& l : lines(read_file(snapshot))) {
        if (l.empty() || l.starts_with('#')) continue;
        const auto tab = l.find('\t');
        want[l.substr(0, tab)] = tab == std::string::npos ? "" : l.substr(tab + 1);
    }
    std::vector<std::string> diffs;
    for (const auto& [k, v] : want) {
        const auto it = got.find(k);
        if (it == got.end())
            diffs.push_back(std::format("{}: want {} got (none)", k, v));
        else if (it->second != v)
            diffs.push_back(std::format("{}: want {} got {}", k, v, it->second));
    }
    for (const auto& [k, v] : got) {
        if (!want.contains(k)) diffs.push_back(std::format("{}: want (none) got {}", k, v));
    }
    if (diffs.empty()) {
        std::println("snapshot: identical ({} records, {})", got["records"], got["digest"]);
        return true;
    }
    std::println("FAIL: the call set differs from the snapshot {} ({} keys):", snapshot.string(), diffs.size());
    for (const auto& d : diffs) std::println("  {}", d);
    return false;
}

}  // namespace

int cmd_bench(Args& args) {
    const auto binary = args.value("--binary");
    const auto vcf_in = args.value("--vcf");
    const auto reference = args.required("--reference");
    const auto bed = args.required("--bed");
    const auto bam = args.value("-b");
    const auto name = args.value("-N");
    const auto truth = args.required("--truth");
    const auto snapshot = args.value("--snapshot");
    const bool update = args.flag("--update-snapshot");
    const auto min_recall = as_number(args.value_or("--min-recall", "0"));
    const bool somatic = args.flag("--somatic");
    const auto min_af = as_number(args.value_or("--min-af", "0"));
    const auto out_opt = args.value("--out");
    const int threads = args.int_or("-th", 16);
    const auto aardvark = args.value_or("--aardvark", "aardvark");
    const int timeout = args.int_or("--timeout", 7200);
    const auto extra = args.finish();
    if (!min_recall || !min_af) die("bench: --min-recall and --min-af expect numbers");
    if (vcf_in.has_value() == binary.has_value()) die("bench: give either --binary (with -b) or --vcf");
    if (binary && !bam) die("bench: --binary needs -b");
    if (update && !snapshot) die("bench: --update-snapshot needs --snapshot");

    std::optional<TempDir> tmp;
    fs::path out;
    if (out_opt) {
        out = *out_opt;
        fs::create_directories(out);
    } else {
        tmp.emplace();
        out = tmp->path();
    }

    // 1. calls
    std::string calls;
    if (binary) {
        std::vector<std::string> cmd = {fs::absolute(*binary).string(), "-G", reference, "-b", *bam};
        if (name) {
            cmd.emplace_back("-N");
            cmd.push_back(*name);
        }
        for (auto& w : words("-c 1 -S 2 -E 3 -g 4")) cmd.push_back(std::move(w));
        cmd.insert(cmd.end(), extra.begin(), extra.end());
        cmd.emplace_back("-th");
        cmd.push_back(std::to_string(threads));
        cmd.push_back(bed);
        const auto calls_path = out / "calls.vcf";
        const auto rc = run(cmd, files_opt(calls_path.string(), (out / "calls.err").string(), timeout));
        std::println("omnivar: exit {}, {:.1f} s ({} threads)", rc.rc, rc.seconds, threads);
        if (rc.rc != 0 || rc.timed_out) {
            std::println("FAIL: omnivar failed: {}", tail(read_file(out / "calls.err"), 300));
            return 1;
        }
        calls = read_file(calls_path);
    } else {
        if (!extra.empty()) die("bench: caller options after -- need --binary");
        calls = read_any(*vcf_in);
    }

    // 2-3. aardvark inputs
    std::size_t kept = 0;
    write_bgzf_indexed(out / "query.vcf.gz", query_records(calls, somatic, *min_af, kept));
    std::println("query: {} PASS records{}{}", kept, somatic ? " with a Somatic STATUS" : "",
                 *min_af > 0 ? std::format(", AF >= {}", *min_af) : "");
    write_bgzf_indexed(out / "truth.vcf.gz", repair_header(read_any(truth)));
    std::string bed3;
    for (const auto& l : lines(read_file(bed))) {
        const auto f = split(l, '\t');
        if (l.empty() || l.starts_with('#') || l.starts_with("track") || l.starts_with("browser") || f.size() < 3)
            continue;
        bed3 += f[0] + "\t" + f[1] + "\t" + f[2] + "\n";
    }
    write_file(out / "regions.bed", bed3);

    // 4. aardvark
    const auto adir = out / "aardvark";
    fs::remove_all(adir);
    const auto ar = run({aardvark, "compare", "-r", reference, "-t", (out / "truth.vcf.gz").string(), "-q",
                         (out / "query.vcf.gz").string(), "-b", (out / "regions.bed").string(), "-o", adir.string(),
                         "--threads", std::to_string(threads)},
                        timeout_opt(timeout));
    if (ar.rc != 0 || !fs::exists(adir / "summary.tsv")) {
        std::println("FAIL: aardvark failed (exit {}): {}", ar.rc, last_line(ar.err));
        return 1;
    }
    const auto summary = lines(read_file(adir / "summary.tsv"));
    if (summary.empty()) die("bench: empty aardvark summary");
    const auto head = split(summary[0], '\t');
    auto col = [&](std::string_view c) {
        const auto it = std::ranges::find(head, c);
        if (it == head.end()) die("bench: aardvark summary has no column " + std::string(c));
        return static_cast<std::size_t>(it - head.begin());
    };
    const auto c_cmp = col("comparison"), c_reg = col("region_label"), c_flt = col("filter"),
               c_type = col("variant_type");
    Snapshot gt;
    std::optional<std::vector<std::string>> all;
    std::println("{:<12} {:>6} {:>6} {:>6} {:>7} {:>7} {:>9}", "aardvark GT", "truth", "TP", "FN", "query", "FP",
                 "recall");
    for (std::size_t i = 1; i < summary.size(); ++i) {
        const auto f = split(summary[i], '\t');
        if (f.size() < head.size() || f[c_cmp] != "GT" || f[c_reg] != "ALL" || f[c_flt] != "ALL") continue;
        const auto& type = f[c_type];
        std::println("{:<12} {:>6} {:>6} {:>6} {:>7} {:>7} {:>9}", type, f[col("truth_total")], f[col("truth_tp")],
                     f[col("truth_fn")], f[col("query_total")], f[col("query_fp")], f[col("metric_recall")]);
        for (const char* k : {"truth_total", "truth_tp", "truth_fn", "query_total", "query_tp", "query_fp"})
            gt.emplace_back("aardvark.GT." + type + "." + k, f[col(k)]);
        if (type == "ALL") all = f;
    }
    if (!all) die("bench: no GT/ALL row in the aardvark summary");
    const double t_total = as_number((*all)[col("truth_total")]).value_or(0);
    const double t_tp = as_number((*all)[col("truth_tp")]).value_or(0);
    const double recall = t_total > 0 ? t_tp / t_total : 0;
    std::println(
        "(FP: scored calls not in the truth set; a false-positive count only when the truth set is exhaustive)");
    bool ok = recall >= *min_recall;
    std::println("{}recall {}/{} = {:.4f} (minimum {})", ok ? "" : "FAIL: ", t_tp, t_total, recall, *min_recall);

    // 5. snapshot
    if (snapshot) {
        auto snap = call_snapshot(calls);
        snap.insert(snap.end(), gt.begin(), gt.end());
        ok = check_snapshot(*snapshot, snap, update) && ok;
    }
    if (out_opt) std::println("outputs in {}", out.string());
    std::println("RESULT: {}", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

}  // namespace testkit
