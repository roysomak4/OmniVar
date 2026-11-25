#include "vcfcmp.hpp"

#include <cmath>
#include <set>
#include <stdexcept>

#include "testkit.hpp"

namespace testkit {

Records parse_records(std::string_view vcf) {
    Records r;
    for (auto& line : lines(vcf)) {
        if (line.empty() || line.starts_with('#')) continue;
        auto f = split(line, '\t');
        if (f.size() < 8) f.resize(8);
        std::optional<std::string> end;
        for (const auto& kv : split(f[7], ';')) {
            const auto eq = kv.find('=');
            if (eq != std::string::npos && kv.substr(0, eq) == "END")
                end = kv.substr(eq + 1);
            else if (eq == std::string::npos && kv == "END")
                end = std::nullopt;
        }
        RecordKey key{f[0], f[1], f[3], f[4], end};
        r.by_key[key].push_back(std::move(f));
        r.order.push_back(std::move(key));
    }
    return r;
}

bool same_value(std::string_view x, std::string_view y) {
    if (x == y) return true;
    const auto a = as_number(x), b = as_number(y);
    if (!a || !b) return false;
    if (std::isinf(*a) || std::isinf(*b) || std::isnan(*a) || std::isnan(*b)) return *a == *b;
    return std::abs(*a - *b) <= 1e-12 * std::max({1.0, std::abs(*a), std::abs(*b)});
}

namespace {

struct Kv {
    std::string key, last;  // last = value, or the key when there is no '='
};

std::vector<Kv> info_items(const std::string& info) {
    std::vector<Kv> out;
    for (const auto& kv : split(info, ';')) {
        const auto eq = kv.find('=');
        if (eq == std::string::npos)
            out.push_back({kv, kv});
        else
            out.push_back({kv.substr(0, eq), kv.substr(eq + 1)});
    }
    return out;
}

std::string list_repr(const std::vector<Kv>& v) {
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + repr(v[i].key);
    return s + "]";
}

}  // namespace

// (INFO key, expected value, C++ value) pairs that are R type-conversion artefacts: read.table
// types a column that is empty in every row as logical NA, so Perl sees TYPE "NA" instead of "REF".
static bool known_r_artefact(std::string_view key, std::string_view want, std::string_view got) {
    return key == "TYPE" && want == "NA" && got == "REF";
}

std::vector<std::string> record_diffs(const Record& want, const Record& got, std::vector<std::string>& known) {
    std::vector<std::string> d;
    auto at = [](const Record& r, std::size_t i) { return i < r.size() ? r[i] : std::string(); };
    for (const auto& [i, col] : {std::pair{2, "ID"}, {5, "QUAL"}, {6, "FILTER"}, {8, "FORMAT"}}) {
        const auto ui = static_cast<std::size_t>(i);
        if (at(want, ui) != at(got, ui))
            d.push_back(std::format("{}: want {} got {}", col, repr(at(want, ui)), repr(at(got, ui))));
    }
    const auto wi = info_items(at(want, 7)), gi = info_items(at(got, 7));
    bool same_keys = wi.size() == gi.size();
    for (std::size_t i = 0; same_keys && i < wi.size(); ++i) same_keys = wi[i].key == gi[i].key;
    if (!same_keys) {
        d.push_back(std::format("INFO keys: want {} got {}", list_repr(wi), list_repr(gi)));
    } else {
        for (std::size_t i = 0; i < wi.size(); ++i) {
            if (same_value(wi[i].last, gi[i].last)) continue;
            if (known_r_artefact(wi[i].key, wi[i].last, gi[i].last))
                known.push_back(wi[i].key + "=" + wi[i].last);
            else
                d.push_back(std::format("INFO {}: want {} got {}", wi[i].key, repr(wi[i].last), repr(gi[i].last)));
        }
    }
    if (want.size() != got.size()) d.push_back(std::format("columns: want {} got {}", want.size(), got.size()));
    const auto names = split(at(want, 8), ':');
    for (std::size_t col = 9; col < std::min(want.size(), got.size()); ++col) {
        const auto ws = split(want[col], ':'), gs = split(got[col], ':');
        if (ws.size() != gs.size()) {
            d.push_back(std::format("sample {}: want {} got {}", col - 8, repr(want[col]), repr(got[col])));
            continue;
        }
        const auto n = std::min(names.size(), ws.size());
        for (std::size_t k = 0; k < n; ++k) {
            const auto wv = split(ws[k], ','), gv = split(gs[k], ',');
            bool ok = wv.size() == gv.size();
            for (std::size_t j = 0; ok && j < wv.size(); ++j) ok = same_value(wv[j], gv[j]);
            if (!ok)
                d.push_back(
                    std::format("FORMAT {} (sample {}): want {} got {}", names[k], col - 8, repr(ws[k]), repr(gs[k])));
        }
    }
    return d;
}

std::string key_label(const RecordKey& k) {
    return std::get<0>(k) + ":" + std::get<1>(k) + " " + std::get<2>(k) + ">" + std::get<3>(k);
}

std::vector<std::string> compare_records(const Records& want, const Records& got, std::vector<std::string>& known) {
    std::set<RecordKey> keys;
    for (const auto& [k, _] : want.by_key) keys.insert(k);
    for (const auto& [k, _] : got.by_key) keys.insert(k);
    static const std::vector<Record> none;
    std::vector<std::string> problems;
    for (const auto& key : keys) {
        const auto wi = want.by_key.find(key), gi = got.by_key.find(key);
        const auto& w = wi == want.by_key.end() ? none : wi->second;
        const auto& g = gi == got.by_key.end() ? none : gi->second;
        if (w.size() != g.size()) {
            problems.push_back(
                std::format("{}: {} expected record(s), {} written", key_label(key), w.size(), g.size()));
            continue;
        }
        for (std::size_t i = 0; i < w.size(); ++i) {
            for (auto& d : record_diffs(w[i], g[i], known)) problems.push_back(key_label(key) + ": " + d);
        }
    }
    return problems;
}

std::optional<std::string> bcftools_check(const std::string& vcf) {
    static const auto bcftools = which("bcftools");
    if (!bcftools) return std::nullopt;
    const auto p = run({*bcftools, "view", "-Ov", "-o", "/dev/null", "-"}, input_opt(vcf));
    const auto msg = strip(p.err);
    if (p.rc || !msg.empty()) return msg.empty() ? std::format("exit {}", p.rc) : lines(msg).front();
    return std::nullopt;
}

std::string chain(std::string_view tsv, int timeout, std::string_view perl_args, bool paired) {
    std::string rows;
    for (const auto& l : lines(tsv)) {
        if (!is_blank(l)) rows += l + "\n";
    }
    const std::string rscript = paired ? "testsomatic.R" : "teststrandbias.R";
    const std::string perl = paired ? "var2vcf_paired.pl" : "var2vcf_valid.pl";
    const auto r = run({"Rscript", (scripts_dir() / rscript).string()}, input_opt(rows, timeout));
    if (r.rc || r.timed_out) throw std::runtime_error(rscript + " failed: " + tail(r.err, 300));
    std::vector<std::string> cmd = {"perl", (scripts_dir() / perl).string()};
    for (auto& w : words(perl_args)) cmd.push_back(std::move(w));
    const auto p = run(cmd, input_opt(r.out, timeout));
    if (p.rc || p.timed_out) throw std::runtime_error(perl + " failed: " + tail(p.err, 300));
    return p.out;
}

}  // namespace testkit
