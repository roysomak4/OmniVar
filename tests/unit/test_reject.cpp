// -p (milestone M10): isGoodVar's first failing rule, and rejected VCF records that leave every
// other record as it is without them.

#include <doctest.h>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "annotate/var2vcf.hpp"
#include "cli/config.hpp"
#include "core/variant.hpp"

namespace {

using namespace omnivar;
using core::RejectReason;

// A variant that passes every isGoodVar rule with the default Config, at AF 0.1 (so the
// AF > 0.30 shortcut does not skip the later rules).
core::Variant good_snv() {
    core::Variant v;
    v.refallele = "A";
    v.varallele = "T";
    v.frequency = 0.1;
    v.hicnt = 10;
    v.mean_position = 20;
    v.mean_quality = 35;
    v.mean_mapping_quality = 60;
    v.high_quality_to_low_quality_ratio = 5;
    v.strand_bias_flag = "2;2";
    v.start_position = 100;
    v.end_position = 100;
    return v;
}

std::optional<RejectReason> reason(const core::Variant& v, const Config& conf = {}, const core::Variant* ref = nullptr,
                                   const std::unordered_set<std::string>& splice = {}) {
    return v.reject_reason(ref, "", splice, conf);
}

}  // namespace

TEST_CASE("reject_reason: one reason per isGoodVar rule") {
    const Config conf;
    core::Variant v = good_snv();
    CHECK_FALSE(reason(v).has_value());
    CHECK(v.is_good_var(nullptr, "", {}, conf));

    SUBCASE("empty reference allele") {
        v.refallele.clear();
        CHECK(reason(v) == RejectReason::no_ref);
    }
    SUBCASE("AF < -f") {
        v.frequency = 0.005;
        CHECK(reason(v) == RejectReason::af);
    }
    SUBCASE("hicnt < -r") {
        v.hicnt = 1;
        CHECK(reason(v) == RejectReason::min_reads);
    }
    SUBCASE("mean position < -P") {
        v.mean_position = 4;
        CHECK(reason(v) == RejectReason::pmean);
    }
    SUBCASE("mean base quality < -q") {
        v.mean_quality = 20;
        CHECK(reason(v) == RejectReason::qual);
    }
    SUBCASE("reference allele has much better MAPQ") {
        core::Variant ref = good_snv();
        ref.hicnt = 100;
        ref.mean_mapping_quality = 60;
        v.mean_mapping_quality = 2;  // d = 2 + 1 + 1 = 4: d - 2 < 5 and ref MQ > 20
        CHECK(reason(v, conf, &ref) == RejectReason::ref_mapq);
    }
    SUBCASE("deletion matching a splice junction") {
        v.refallele = "AC";
        v.varallele = "A";
        v.end_position = 101;
        CHECK(reason(v, conf, nullptr, {"100-101"}) == RejectReason::splice);
    }
    SUBCASE("quality ratio < -o") {
        v.high_quality_to_low_quality_ratio = 1;
        CHECK(reason(v) == RejectReason::qratio);
    }
    SUBCASE("MAPQ < -O") {
        Config c;
        c.mapq = 30;
        v.mean_mapping_quality = 20;
        CHECK(reason(v, c) == RejectReason::mapq);
        v.frequency = 0.4;  // AF > 0.30 is good before the MAPQ rule
        CHECK_FALSE(reason(v, c).has_value());
    }
    SUBCASE("monomer MSI") {
        v.msi = 15;
        v.msint = 1;
        CHECK(reason(v) == RejectReason::msi_monomer);
    }
    SUBCASE("non-monomer MSI") {
        v.msi = 12;
        v.msint = 2;
        v.frequency = 0.05;
        CHECK(reason(v) == RejectReason::msi);
    }
    SUBCASE("strand bias 2;1") {
        v.strand_bias_flag = "2;1";
        CHECK(reason(v) == RejectReason::strand_bias);
        v.frequency = 0.25;
        CHECK_FALSE(reason(v).has_value());
    }
}

TEST_CASE("reject_reason: the first failing rule wins, in the upstream order") {
    core::Variant v = good_snv();
    v.frequency = 0.005;  // AF
    v.hicnt = 1;          // and min reads
    v.mean_quality = 10;  // and quality
    CHECK(reason(v) == RejectReason::af);
    v.frequency = 0.1;
    CHECK(reason(v) == RejectReason::min_reads);
    v.hicnt = 10;
    CHECK(reason(v) == RejectReason::qual);
}

TEST_CASE("reject filter IDs and descriptions") {
    Config conf;
    conf.freq = 0.05;
    const auto d = core::reject_filter_descriptions(conf);
    CHECK(d.size() == 11);  // every reason except no_ref (never written)
    CHECK(d.front().first == "REJ_AF");
    CHECK(d.front().second.contains("0.05"));
    CHECK(core::reject_id(RejectReason::strand_bias) == "REJ_Bias");
}

namespace {

std::vector<std::string> lines_of(const std::string& path) {
    std::ifstream in(path);
    REQUIRE(in.good());
    std::vector<std::string> out;
    for (std::string l; std::getline(in, l);) {
        if (!l.empty()) {
            out.push_back(l);
        }
    }
    return out;
}

std::vector<std::string> records(const std::string& vcf, bool drop_rejected) {
    std::vector<std::string> out;
    std::istringstream in(vcf);
    for (std::string l; std::getline(in, l);) {
        if (!drop_rejected || l.find("\tREJ_") == std::string::npos) {
            out.push_back(l);
        }
    }
    return out;
}

}  // namespace

TEST_CASE("Var2Vcf: rejected records leave the other records unchanged") {
    // Calls (the jar without -p) and this port's rejected rows for a telomeric region with many
    // low-AF alleles: rejected rows sit at the same starts as calls, between clustered SNVs, and
    // above calls in AF order.
    const std::string name = "Simple;hg19.fa;Colo829-18_S3-sort.bam;chr1;9922-10122;-p;.txt";
    const auto calls = lines_of(std::string(OMNIVAR_SOURCE_DIR) + "/tests/golden/expected_nopileup/full/" + name);
    const auto rejected = lines_of(std::string(OMNIVAR_SOURCE_DIR) + "/tests/golden/expected_rejected/full/" + name);
    REQUIRE(calls.size() == 8);
    REQUIRE(rejected.size() > 100);

    for (int variant = 0; variant < 3; ++variant) {
        CAPTURE(variant);
        annotate::Var2VcfOptions opt;
        if (variant == 1) {
            opt.top_allele_only = true;
        } else if (variant == 2) {
            opt.cluster_bp = {200, "200"};
        }
        annotate::Var2Vcf plain(opt);
        annotate::Var2Vcf with_rejected(opt);
        for (const auto& l : calls) {
            plain.add_row(l);
        }
        // interleave as the caller does: rejected rows among the calls, by start
        // (one call after every 17 rejected rows; rows are regrouped by start inside Var2Vcf)
        std::size_t c = 0;
        for (std::size_t i = 0; i < rejected.size(); ++i) {
            const std::string& r = rejected[i];
            const auto tab = r.rfind('\t');
            if (c < calls.size() && i >= c * 17) {
                with_rejected.add_row(calls[c++]);
            }
            with_rejected.add_row(r.substr(0, tab), r.substr(tab + 1));
        }
        while (c < calls.size()) {
            with_rejected.add_row(calls[c++]);
        }
        const auto want = records(plain.records("s"), false);
        const auto got_all = records(with_rejected.records("s"), false);
        const auto got = records(with_rejected.records("s"), true);
        CHECK(got == want);
        CHECK(got_all.size() >= want.size() + rejected.size() - 5);  // a few duplicate keys may collapse
        for (const auto& l : got_all) {
            if (l.find("\tREJ_") != std::string::npos) {
                CHECK(l.find("Cluster") == std::string::npos);
            }
        }
    }
}
