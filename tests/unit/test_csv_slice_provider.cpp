#include <doctest.h>
#include <filesystem>

#include "io/csv_slice_provider.hpp"

using senri::io::CsvSliceProvider;
using senri::io::ReferenceError;

namespace {
CsvSliceProvider make(std::string_view text) {
    auto p = CsvSliceProvider::from_text(text, "test.csv");
    REQUIRE(p);
    return std::move(*p);
}
}  // namespace

TEST_CASE("serves a range inside one slice") {
    auto p = make("chr1,100,109,ACGTACGTAC\n");
    CHECK(p.sub_seq("chr1", 100, 109).value() == "ACGTACGTAC");
    CHECK(p.sub_seq("chr1", 102, 104).value() == "GTA");
    CHECK(p.sub_seq("chr1", 109, 109).value() == "C");
}

TEST_CASE("uses the slice with the greatest start <= query start") {
    auto p = make("chr1,100,109,AAAAAAAAAA\nchr1,105,124,CCCCCCCCCCCCCCCCCCCC\n");
    CHECK(p.sub_seq("chr1", 104, 106).value() == "AAA");  // only first slice starts <= 104
    CHECK(p.sub_seq("chr1", 105, 120).value() == std::string(16, 'C'));
    // floor slice for 104 is the first one, which ends at 109: not served even though the
    // union of slices covers it (a floorEntry lookup).
    CHECK_FALSE(p.sub_seq("chr1", 104, 115));
}

TEST_CASE("slice end comes from the sequence length, not the CSV end column") {
    auto p = make("chr1,100,999,ACGT\n");
    CHECK(p.sub_seq("chr1", 100, 103).value() == "ACGT");
    auto r = p.sub_seq("chr1", 100, 104);
    REQUIRE_FALSE(r);
    CHECK(r.error().kind == ReferenceError::Kind::no_slice);
    CHECK(r.error().message == "No region: >chr1:100-104 in file test.csv");
}

TEST_CASE("missing contig, query before first slice") {
    auto p = make("chr1,100,103,ACGT\n");
    CHECK_FALSE(p.sub_seq("chr2", 100, 101));
    CHECK_FALSE(p.sub_seq("chr1", 99, 101));
}

TEST_CASE("repeated chr,start keeps the last line; blank lines and CRLF are skipped") {
    auto p = make("chr1,100,103,AAAA\r\n\nchr1,100,103,GGGG\n");
    CHECK(p.sub_seq("chr1", 100, 103).value() == "GGGG");
}

TEST_CASE("case is preserved (caller upper-cases)") {
    auto p = make("chr1,1,4,acgN\n");
    CHECK(p.sub_seq("chr1", 1, 4).value() == "acgN");
}

TEST_CASE("malformed lines are errors") {
    CHECK_FALSE(CsvSliceProvider::from_text("chr1,100\n", "x.csv"));
    CHECK_FALSE(CsvSliceProvider::from_text("chr1,1x,4,ACGT\n", "x.csv"));
}

TEST_CASE("every regression reference CSV loads") {
    const std::filesystem::path dir = SENRI_SOURCE_DIR "/tests/regression/fastas";
    REQUIRE(std::filesystem::is_directory(dir));
    int n = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() != ".csv") {
            continue;
        }
        INFO(entry.path().string());
        CHECK(CsvSliceProvider::open(entry.path().string()));
        ++n;
    }
    CHECK(n == 9);
}
