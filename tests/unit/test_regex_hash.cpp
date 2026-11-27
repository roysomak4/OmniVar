#include <doctest.h>
#include <fstream>
#include <string>
#include <vector>

#include "util/hash_order.hpp"
#include "util/number_parse.hpp"
#include "util/regex_engine.hpp"

using namespace senri::util;

namespace {

std::string group1(std::string_view pattern, std::string_view text) {
    auto re = Regex::compile(pattern);
    REQUIRE_MESSAGE(re, (re ? "" : re.error()));
    auto m = re->find(text);
    if (!m) {
        return "<nomatch>";
    }
    auto g = m->group(text, 1);
    return g ? std::string(*g) : "<null>";
}

std::vector<std::string> split(std::string_view text, std::string_view regex) {
    auto r = regex_split(text, regex);
    REQUIRE(r);
    return {r->begin(), r->end()};
}

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

TEST_CASE("Upstream sample-name patterns (Patterns.SAMPLE_PATTERN / SAMPLE_PATTERN2)") {
    // Patterns SAMPLE_PATTERN / SAMPLE_PATTERN2 (string-literal escapes already removed).
    constexpr std::string_view kSample = "([^\\/\\._]+).sorted[^\\/]*.bam";
    constexpr std::string_view kSample2 = "([^\\/]+)[_\\.][^\\/]*bam";
    CHECK(group1(kSample, "/data/NA12878.sorted.bam") == "NA12878");
    CHECK(group1(kSample, "/data/tumor_1.sorted.dedup.bam") == "1");  // [^/._]+ stops at '_'
    CHECK(group1(kSample, "x.bam") == "<nomatch>");
    CHECK(group1(kSample2, "tests/regression/bams/Colo829-18_S3-sort.bam") == "Colo829-18_S3-sort");
    CHECK(group1(kSample2, "/a/b/Oncospan_DNA.realigned.bam") == "Oncospan_DNA.realigned");
    CHECK(group1(kSample2, "t.bam|n.bam") == "t.bam|n");  // raw -b value, greedy group
    CHECK(group1(kSample2, "nobam") == "<nomatch>");
    // Help text default for -n (lazy quantifier)
    CHECK(group1("([^\\/\\._]+?)_[^\\/]*.bam", "/x/sampleA_L001.bam") == "sampleA");
}

TEST_CASE("regex features") {
    CHECK(group1("a(b*)c", "xxabbbcx") == "bbb");
    CHECK(group1("a(b*?)b", "abbb") == "");
    CHECK(group1("(x|xy)z", "xyz") == "xy");
    CHECK(group1("(?:ab)+(c)", "ababc") == "c");
    CHECK(group1("^(\\d+)$", "12345") == "12345");
    CHECK(group1("^(\\d+)$", "123a") == "<nomatch>");
    CHECK(group1("(\\d{2,3})", "a12345") == "123");
    CHECK(group1("(\\d{2,})", "a12345") == "12345");
    CHECK(group1("(\\w+)\\s+(\\w+)", "  foo   bar") == "foo");
    CHECK(group1("(a)?b", "b") == "<null>");
    CHECK(group1("([a-c-])+", "zz-ab-q") == "-");  // last iteration's capture
    CHECK(group1("(.)$", "ab\n") == "b");          // $ before final line terminator
    CHECK(group1("([^]a])", "]]ab") == "b");
    CHECK(group1("(\\.)", "a.b") == ".");
    CHECK_FALSE(Regex::compile("a(?=b)"));
    CHECK_FALSE(Regex::compile("(a"));
    CHECK_FALSE(Regex::compile("a)"));
    CHECK_FALSE(Regex::compile("[ab"));
    CHECK_FALSE(Regex::compile("*a"));
    CHECK_FALSE(Regex::compile("a++"));
    CHECK_FALSE(Regex::compile("\\b"));
}

TEST_CASE("String.split semantics") {
    using V = std::vector<std::string>;
    CHECK(split("a\tb\tc", "\t") == V{"a", "b", "c"});
    CHECK(split("a\t\tc", "\t") == V{"a", "", "c"});
    CHECK(split("a\tb\t\t", "\t") == V{"a", "b"});    // trailing empties removed
    CHECK(split("\ta\tb", "\t") == V{"", "a", "b"});  // leading empty kept
    CHECK(split("\t", "\t").empty());                 // only delimiters -> length 0
    CHECK(split("", "\t") == V{""});                  // no match -> input itself
    CHECK(split("abc", ",") == V{"abc"});
    CHECK(split("a,b;c", "[,;]") == V{"a", "b", "c"});
    CHECK(split("a  b c", " +") == V{"a", "b", "c"});
    CHECK(split("a|b", "\\|") == V{"a", "b"});
    CHECK(split("a.b", ".").empty());             // '.' is a regex: everything matches
    CHECK(split("abc", "") == V{"a", "b", "c"});  // zero-width matches
    CHECK(split("1,2,3", ",") == V{"1", "2", "3"});
}

TEST_CASE("HashMap iteration order matches the JVM (tests/unit/data/hashmap_order_vectors.txt)") {
    std::ifstream in(SENRI_SOURCE_DIR "/tests/unit/data/hashmap_order_vectors.txt");
    REQUIRE(in);
    std::string line;
    int cases = 0;
    int skipped_tree = 0;
    while (std::getline(in, line)) {
        auto f = fields(line, '\t');
        REQUIRE(f.size() == 3);
        auto keys = fields(f[1], ',');
        std::vector<std::int32_t> hashes;
        for (const auto& k : keys) {
            hashes.push_back(f[0] == "S" ? string_hash_code(k) : integer_hash_code(parse_int(k).value()));
        }
        auto got = hashmap_order(hashes);
        if (!got) {
            ++skipped_tree;  // a bin was treeified; order not modelled
            continue;
        }
        std::vector<std::size_t> want;
        for (const auto& s : fields(f[2], ',')) {
            want.push_back(static_cast<std::size_t>(parse_int(s).value()));
        }
        INFO("case " << cases << " kind " << f[0] << " n=" << keys.size());
        CHECK(*got == want);
        ++cases;
    }
    CHECK(cases > 300);
    MESSAGE("HashMap vectors: " << cases << " checked, " << skipped_tree << " with tree bins skipped");
}

TEST_CASE("String.hashCode") {
    CHECK(string_hash_code("") == 0);
    // Values from the JVM
    CHECK(string_hash_code("chr1") == 3052836);
    CHECK(string_hash_code("chrX") == 3052875);
    CHECK(string_hash_code("chr1_KI270706v1_random") == 382476732);
    CHECK(string_hash_code("Oncospan_DNA") == -1765068243);
}
