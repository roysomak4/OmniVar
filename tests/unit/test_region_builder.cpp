#include <doctest.h>

#include "cli/config.hpp"
#include "region/region_builder.hpp"

using namespace omnivar;

namespace {

Config test_config(BedRowFormat format, int chr_column) {
    Config c;
    c.bed_row_format = format;
    c.zero_based = true;
    c.column_for_chromosome = chr_column;
    c.number_nucleotide_to_extend = 0;
    c.delimiter = "\t";
    return c;
}

const ChrLengths kNoChrs;

RegionBuilder spy(const Config& c) {
    auto b = RegionBuilder::create(kNoChrs, c);
    REQUIRE(b);
    b->set_chromosome_override("Y");  // Mockito.doReturn("Y").when(spy).correctChromosome(...)
    return std::move(*b);
}

}  // namespace

TEST_CASE("RegionBuilderTest.testNonAmpToRegions") {
    const Config default_config = test_config(BedRowFormat{0, 1, 2, 1, 2, 3}, 0);
    const std::vector<std::string> raws_1 = {
        "Y\t59023460\t59028460\tY_59023460_1", "Y\t59028310\t59033310\tY_59028310_2",
        "Y\t59033160\t59038160\tY_59033160_3", "Y\t59038010\t59043010\tY_59038010_4",
        "Y\t59042860\t59047860\tY_59042860_5",
    };
    const Segments result_1 = {
        {Region{"Y", 59023461, 59028460, "Y_59023460_1"}}, {Region{"Y", 59028311, 59033310, "Y_59028310_2"}},
        {Region{"Y", 59033161, 59038160, "Y_59033160_3"}}, {Region{"Y", 59038011, 59043010, "Y_59038010_4"}},
        {Region{"Y", 59042861, 59047860, "Y_59042860_5"}},
    };
    CHECK(spy(default_config).build_regions(raws_1, true).value() == result_1);

    const std::vector<std::string> raws_2 = {
        "Y\t59023460\t59028460\t59028460", "Y\t59028310\t59033310\t59033310", "Y\t59033160\t59038160\t59038160",
        "Y\t59038010\t59043010\t59043010", "Y\t59042860\t59047860\t59047860",
    };
    const Segments result_2 = {{}, {}, {}, {}, {}};
    const Config conf = test_config(BedRowFormat{0, 1, 2, 1, 2, 3}, -1);
    CHECK(spy(conf).build_regions(raws_2, true).value() == result_2);
}

TEST_CASE("RegionBuilderTest.testSingleRegion") {
    const std::pair<const char*, Region> cases[] = {
        {"Y:59023460-59028460:Y_59023460_1", Region{"Y", 59023461, 59028460, "Y_59023460_1"}},
        {"Y:59028310-59033310:Y_59028310_2", Region{"Y", 59028311, 59033310, "Y_59028310_2"}},
        {"Y:59033160-59038160:Y_59033160_3", Region{"Y", 59033161, 59038160, "Y_59033160_3"}},
        {"Y:59038010-59043010:Y_59038010_4", Region{"Y", 59038011, 59043010, "Y_59038010_4"}},
        {"Y:59042860-59047860:Y_59042860_5", Region{"Y", 59042861, 59047860, "Y_59042860_5"}},
    };
    for (const auto& [roi, expected] : cases) {
        Config c = test_config(BedRowFormat{0, 1, 2, 1, 2, 3}, 0);
        c.region_of_interest = roi;
        CHECK(spy(c).build_region_from_configuration().value() == Segments{{expected}});
    }
}

TEST_CASE("RegionBuilderTest.testAmp") {
    const Config default_config = test_config(BedRowFormat{0, 1, 2, 1, 2, 3}, 0);
    const std::vector<std::string> raws = {
        "Y\t59023460\t59028460\tY_59023460_1\t0\t.\t59023460\t59028460",
        "Y\t59028310\t59033310\tY_59028310_2\t0\t.\t59028310\t59033310",
        "Y\t59033160\t59038160\tY_59033160_3\t0\t.\t59033160\t59038160",
        "Y\t59038010\t59043010\tY_59038010_4\t0\t.\t59038010\t59043010",
        "Y\t59042860\t59047860\tY_59042860_5\t0\t.\t59042860\t59047860",
    };
    const Segments expected = {{
        Region{"Y", 59023461, 59028460, "Y_59023460_1", 59023461, 59028460},
        Region{"Y", 59028311, 59033310, "Y_59028310_2", 59028311, 59033310},
        Region{"Y", 59033161, 59038160, "Y_59033160_3", 59033161, 59038160},
        Region{"Y", 59038011, 59043010, "Y_59038010_4", 59038011, 59043010},
        Region{"Y", 59042861, 59047860, "Y_59042860_5", 59042861, 59047860},
    }};
    CHECK(spy(default_config).build_amp_regions(raws, true).value() == expected);
}

// ---- Beyond the upstream unit tests -------------------------------------------------------------

TEST_CASE("correctChromosome adds or strips chr when the name is not in the BAM header") {
    const ChrLengths chrs = {{"chr1", 1000}, {"2", 1000}};
    Config c;
    auto b = RegionBuilder::create(chrs, c).value();
    CHECK(b.correct_chromosome("chr1") == "chr1");
    CHECK(b.correct_chromosome("1") == "chr1");
    CHECK(b.correct_chromosome("chr2") == "2");
    CHECK(b.correct_chromosome("2") == "2");
    CHECK(b.correct_chromosome("chrZ") == "Z");  // corrected even if absent either way
}

TEST_CASE("regression harness BED with -c 1 -S 2 -E 3 -g 4 -z") {
    Config c;
    c.column_for_chromosome = 0;
    c.bed_row_format = BedRowFormat{0, 1, 2, 1, 2, 3};
    c.zero_based = true;
    const ChrLengths chrs = {{"chr1", 248956422}};
    auto b = RegionBuilder::create(chrs, c).value();
    auto segs = b.build_regions({"chr1\t100\t200\ttestbed"}, true).value();
    CHECK(segs == Segments{{Region{"chr1", 101, 200, "testbed"}}});
}

TEST_CASE("-x extension, start == end, single-position -R, commas in -R") {
    Config c;
    c.column_for_chromosome = 0;
    c.bed_row_format = BedRowFormat{0, 1, 2, 1, 2, 3};
    c.number_nucleotide_to_extend = 10;
    const ChrLengths chrs = {{"chr1", 1000}};
    auto b = RegionBuilder::create(chrs, c).value();
    CHECK(b.build_regions({"chr1\t100\t200\tg"}, std::nullopt).value() == Segments{{Region{"chr1", 90, 210, "g"}}});
    c.region_of_interest = "chr1:1,000";
    auto b2 = RegionBuilder::create(chrs, c).value();
    CHECK(b2.build_region_from_configuration().value() == Segments{{Region{"chr1", 990, 1010, "chr1"}}});
}

TEST_CASE("4-column BED without -c switches to the custom format for that line and every later one") {
    Config c;  // default refGene format, no -c
    const ChrLengths chrs = {{"chr1", 1000}};
    auto b = RegionBuilder::create(chrs, c).value();
    // Custom format reads thick start from column 4 and thick end from column 2 (as upstream wrote it):
    // here col4 = 150 > col2 = 100 -> cdsEnd(200) > thickEnd(100) -> break -> empty list.
    auto segs = b.build_regions({"chr1\t100\t200\t150"}, std::nullopt).value();
    CHECK(segs == Segments{{}});
    // A 4-column line whose column 4 is not a number fails like upstream (NumberFormatException).
    CHECK_FALSE(b.build_regions({"chr1\t100\t200\tGENE"}, std::nullopt));
}

TEST_CASE("readBedFile: comments, amplicon detection, zero-based default") {
    Config c;
    auto bed =
        read_bed_text(c, "#c\ntrack x\nbrowser y\nchr1\t100\t200\tA\t0\t.\t120\t180\r\nchr1\t300\t400\tB\n").value();
    CHECK(bed.lines.size() == 2);
    CHECK(bed.amplicon_parameters.value() == "10:0.95");
    CHECK(bed.zero_based.value());
    Config z;
    z.zero_based = false;
    auto bed2 = read_bed_text(z, "chr1\t100\t200\tA\t0\t.\t120\t180\n").value();
    CHECK_FALSE(bed2.zero_based.value());
    auto bed3 = read_bed_text(c, "chr1\t100\t200\tA\t0\t.\t90\t180\n").value();  // insert outside region
    CHECK_FALSE(bed3.amplicon_parameters);
}

TEST_CASE("buildAmpRegions groups by chromosome in JDK HashMap order") {
    Config c;
    const ChrLengths chrs = {{"chr1", 1000}, {"chr2", 1000}, {"chr10", 1000}};
    auto b = RegionBuilder::create(chrs, c).value();
    auto segs = b.build_amp_regions(
                     {
                         "chr2\t10\t90\ta\t0\t.\t20\t80",
                         "chr1\t10\t90\tb\t0\t.\t50\t80",
                         "chr1\t10\t90\tc\t0\t.\t20\t60",  // overlaps the next one after sorting -> same group
                         "chr10\t10\t90\td\t0\t.\t20\t80",
                     },
                     false)
                    .value();
    REQUIRE(segs.size() == 3);
    // Verified on the JVM: a HashMap with puts chr2, chr1, chr10 iterates [chr10, chr1, chr2].
    CHECK(segs[0][0].chr == "chr10");
    REQUIRE(segs[1].size() == 2);  // chr1 regions sorted by insert start; 20-60 and 50-80 overlap
    CHECK(segs[1][0].gene == "c");
    CHECK(segs[1][1].gene == "b");
    CHECK(segs[2][0].chr == "chr2");
}
