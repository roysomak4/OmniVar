#include <doctest.h>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "stats/fisher.hpp"
#include "util/num_format.hpp"
#include "util/number_parse.hpp"
#include "util/ordered_maps.hpp"

using namespace omnivar;

namespace {
std::vector<std::string> fields(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        auto p = s.find(sep, start);
        out.push_back(s.substr(start, p - start));
        if (p == std::string::npos) {
            break;
        }
        start = p + 1;
    }
    return out;
}
}  // namespace

TEST_CASE("FisherExact matches the jar (tests/golden/fisher_vectors.tsv)") {
    std::ifstream in(OMNIVAR_SOURCE_DIR "/tests/golden/fisher_vectors.tsv");
    REQUIRE(in);
    std::string line;
    int rows = 0;
    int failures = 0;
    while (std::getline(in, line)) {
        auto f = fields(line, '\t');
        REQUIRE(f.size() == 8);
        const int a = util::parse_int(f[0]).value();
        const int b = util::parse_int(f[1]).value();
        const int c = util::parse_int(f[2]).value();
        const int d = util::parse_int(f[3]).value();
        stats::FisherExact fe(a, b, c, d);
        const std::string got = util::double_to_string(fe.p_value()) + "\t" + fe.odd_ratio() + "\t" +
                                util::double_to_string(fe.p_value_less()) + "\t" +
                                util::double_to_string(fe.p_value_greater());
        const std::string want = f[4] + "\t" + f[5] + "\t" + f[6] + "\t" + f[7];
        if (got != want) {
            ++failures;
            if (failures <= 20) {
                FAIL_CHECK(a << " " << b << " " << c << " " << d << ": got " << got << " want " << want);
            }
        }
        ++rows;
    }
    CHECK(rows > 4000);
    CHECK(failures == 0);
}

TEST_CASE("FisherExact (R flavour) matches R fisher.test (tests/golden/fisher_r_vectors.tsv)") {
    std::ifstream in(OMNIVAR_SOURCE_DIR "/tests/golden/fisher_r_vectors.tsv");
    REQUIRE(in);
    // R prints 15 significant digits, which identifies a value rounded to 5 dp exactly.
    const auto num = [](const std::string& s) {
        return s == "Inf" ? std::numeric_limits<double>::infinity() : util::parse_double(s).value();
    };
    std::string line;
    int rows = 0;
    int failures = 0;
    while (std::getline(in, line)) {
        auto f = fields(line, '\t');
        REQUIRE(f.size() == 10);
        const int a = util::parse_int(f[0]).value();
        const int b = util::parse_int(f[1]).value();
        const int c = util::parse_int(f[2]).value();
        const int d = util::parse_int(f[3]).value();
        stats::FisherExact fe(a, b, c, d, stats::FisherFlavor::r);
        const double got[] = {fe.p_value(), fe.odds_ratio(), fe.p_value_greater(), fe.p_value_less()};
        const double want[] = {num(f[4]), num(f[5]), num(f[6]), num(f[8])};
        const bool estimates_agree = num(f[5]) == num(f[7]) && num(f[5]) == num(f[9]);
        if (got[0] != want[0] || got[1] != want[1] || got[2] != want[2] || got[3] != want[3] || !estimates_agree) {
            ++failures;
            if (failures <= 20) {
                FAIL_CHECK(a << " " << b << " " << c << " " << d << ": got " << got[0] << " " << got[1] << " " << got[2]
                             << " " << got[3] << " want " << line);
            }
        }
        ++rows;
    }
    CHECK(rows > 4000);
    CHECK(failures == 0);
}

TEST_CASE("HashOrderMap iterates like the JDK HashMap (tests/golden/hashmap_order_vectors.txt)") {
    std::ifstream in(OMNIVAR_SOURCE_DIR "/tests/golden/hashmap_order_vectors.txt");
    REQUIRE(in);
    std::string line;
    int checked = 0;
    while (std::getline(in, line)) {
        auto f = fields(line, '\t');
        auto keys = fields(f[1], ',');
        std::vector<std::string> want;
        for (const auto& i : fields(f[2], ',')) {
            want.push_back(keys[static_cast<std::size_t>(util::parse_int(i).value())]);
        }
        std::vector<std::string> got;
        bool tree = false;
        if (f[0] == "S") {
            util::HashOrderMap<std::string, int> m;
            for (std::size_t i = 0; i < keys.size(); ++i) {
                m.put(keys[i], static_cast<int>(i));
            }
            tree = m.tree_bins();
            m.for_each([&](const std::string& k, int) { got.push_back(k); });
        } else {
            util::HashOrderMap<int, int> m;
            for (std::size_t i = 0; i < keys.size(); ++i) {
                m.put(util::parse_int(keys[i]).value(), static_cast<int>(i));
            }
            tree = m.tree_bins();
            m.for_each([&](const int& k, int) { got.push_back(std::to_string(k)); });
        }
        if (tree) {
            continue;  // tree bins are not modelled (reported by tree_bins())
        }
        CHECK(got == want);
        ++checked;
    }
    CHECK(checked > 300);
}

TEST_CASE("HashOrderMap: put keeps position, erase unlinks, values are stable") {
    util::HashOrderMap<int, std::string> m;
    for (int i = 0; i < 40; ++i) {
        m.put(i * 16, "v" + std::to_string(i));  // same low bits: long chains, forces resizes
    }
    std::string* p = m.get(16);
    REQUIRE(p != nullptr);
    for (int i = 40; i < 200; ++i) {
        m.put(i, "w");
    }
    CHECK(*p == "v1");  // address stable across resizes
    auto before = m.keys();
    m.put(16, "x");
    CHECK(m.keys() == before);
    CHECK(m.erase(16));
    CHECK_FALSE(m.contains(16));
    CHECK(m.size() == before.size() - 1);
    m.put(16, "y");  // re-inserted at the end of its bucket chain
    CHECK(m.contains(16));
}

TEST_CASE("LinkedMap keeps insertion order; removed keys re-insert at the end") {
    util::LinkedMap<std::string, int> m;
    m.put("b", 1);
    m.put("a", 2);
    m.put("c", 3);
    m.put("a", 4);
    CHECK(m.keys() == std::vector<std::string>{"b", "a", "c"});
    m.erase("b");
    m.put("b", 5);
    CHECK(m.keys() == std::vector<std::string>{"a", "c", "b"});
    CHECK(*m.get("a") == 4);
}

TEST_CASE("DoubleStream.sum uses Kahan summation (values from the JVM)") {
    const std::vector<double> tenths(10, 0.1);
    CHECK(stats::stream_sum(tenths) == 1.0);                     // JVM: 1.0; a naive loop gives 0.9999999999999999
    CHECK(stats::stream_sum({1.0, 1e100, 1.0, -1e100}) == 0.0);  // JVM: 0.0
}
