#include <doctest.h>
#include <string_view>
#include <vector>

#include "cli/cli.hpp"
#include "version.hpp"

using senri::cli::Action;
using senri::cli::parse_command_line;

namespace {
auto parse(std::vector<std::string_view> args) {
    return parse_command_line(std::span<const std::string_view>(args), 8);
}
}  // namespace

// Parser behaviours that do not reach a Configuration (help, warnings, added options).

TEST_CASE("help: no options, -H, -?, and missing -b") {
    auto none = parse({});
    REQUIRE(none);
    CHECK(none->action == Action::help);
    CHECK(none->warnings == "Missing required option(s): b\n");
    CHECK(parse({"-b", "x.bam", "-H"})->action == Action::help);
    CHECK(parse({"-b", "x.bam", "-?"})->action == Action::help);
    CHECK(parse({"-b", "x.bam", "--?"})->action == Action::help);
    CHECK(parse({"a.bed"})->action == Action::help);  // no -b
    CHECK(parse({"--help"})->action == Action::help);
}

TEST_CASE("-h is the header flag, not help") {
    auto r = parse({"-b", "x.bam", "-h", "a.bed", "-G", "r.fa"});
    REQUIRE(r);
    CHECK(r->action == Action::run);
    CHECK(r->config.print_header);
    CHECK(r->config.bed.value() == "a.bed");
}

TEST_CASE("--version") {
    CHECK(parse({"--version"})->action == Action::version);
    auto v = senri::cli::version_string();
    CHECK(v.starts_with("senri " + std::string(senri::kVersion)));
    CHECK(v.contains("htslib "));
}

TEST_CASE("commons-cli errors") {
    CHECK(parse({"-b", "x.bam", "--bogus"}).error().message == "Unrecognized option: --bogus");
    CHECK(parse({"-b", "x.bam", "-f0.001"}).error().message == "Unrecognized option: -f0.001");
    CHECK(parse({"-b", "x.bam", "-f"}).error().message == "Missing argument for option: f");
    CHECK(parse({"-b", "x.bam", "-r", "-3"}).error().message == "Missing argument for option: r");
    CHECK(parse({"-b", "x.bam", "-f", "1e-3"}).error().message == "For input string: \"1e-3\"");
    CHECK(parse({"-b", "x.bam", "-z", "a.bed"}).error().message == "For input string: \"a.bed\"");
}

TEST_CASE("missing -G prints the upstream notice and uses the hg19 default") {
    auto r = parse({"-b", "x.bam", "a.bed"});
    REQUIRE(r);
    CHECK(r->warnings == "Reference file path wasn't set (option -G). Will be used the default fasta path.\n");
    CHECK(r->config.fasta == senri::Config::kHg19);
}

TEST_CASE("dev options are removed before the upstream rules") {
    auto r = parse({"--dev-tsv-to-vcf", "rows.tsv", "-b", "x.bam", "-G", "r.fa"});
    REQUIRE(r);
    CHECK(r->action == Action::tsv_to_vcf);
    CHECK(r->tsv_path == "rows.tsv");
    // after "--" they are positional
    auto p = parse({"-b", "x.bam", "-G", "r.fa", "--", "--version"});
    REQUIRE(p);
    CHECK(p->action == Action::run);
    CHECK(p->config.bed.value() == "--version");
}

TEST_CASE("-th without a value uses the processor count") {
    CHECK(parse({"-b", "x.bam", "-G", "r", "-th"})->config.threads == 8);
    CHECK(parse({"-b", "x.bam", "-G", "r", "-th", "0"})->config.threads == 1);
    CHECK(parse({"-b", "x.bam", "-G", "r"})->config.threads == 1);
}

TEST_CASE("output and VCF annotation options added by this port") {
    auto d = parse({"-b", "x.bam", "a.bed", "-G", "r.fa"});
    REQUIRE(d);
    CHECK_FALSE(d->config.output_tsv);  // VCF is the default output
    CHECK(d->config.vcf.min_af.text == "0.02");
    CHECK(d->config.vcf.min_qual.value == 22.5);
    CHECK(d->config.vcf.pstd);
    CHECK_FALSE(d->config.vcf.top_allele_only);

    auto r = parse({"-b", "x.bam", "--tsv", "--filter-min-af=0.050", "--filter-min-depth", "10", "a.bed",
                    "--no-filter-pstd", "--vcf-top-allele-only", "-G", "r.fa", "-f", "0.001"});
    REQUIRE(r);
    CHECK(r->action == Action::run);
    CHECK(r->config.output_tsv);
    CHECK(r->config.vcf.min_af.value == 0.05);
    CHECK(r->config.vcf.min_af.text == "0.050");  // tag name keeps the text, as Perl does: f0.050
    CHECK(r->config.vcf.min_depth.value == 10);
    CHECK_FALSE(r->config.vcf.pstd);
    CHECK(r->config.vcf.top_allele_only);
    // the upstream options around them are unaffected
    CHECK(r->config.freq == 0.001);
    CHECK(r->config.bed.value() == "a.bed");

    CHECK(parse({"-b", "x.bam", "a.bed", "--filter-min-qual"}).error().message ==
          "Missing argument for option: --filter-min-qual");
    CHECK(parse({"-b", "x.bam", "a.bed", "--filter-min-qual", "high"}).error().message ==
          "Invalid number for option --filter-min-qual: high");
}

TEST_CASE("--dev-tsv-to-vcf needs no caller options") {
    auto r = parse({"--dev-tsv-to-vcf", "rows.tsv", "--filter-cluster-bp", "10"});
    REQUIRE(r);
    CHECK(r->action == Action::tsv_to_vcf);
    CHECK(r->tsv_path == "rows.tsv");
    CHECK(r->config.vcf.cluster_bp.text == "10");
}

// var2vcf_paired.pl defaults apply to thresholds the user did not set.
TEST_CASE("paired VCF options: Perl defaults unless given") {
    const std::vector<std::string_view> plain = {"-b", "t.bam|n.bam", "-G", "ref.fa", "-R", "chr1:1-10"};
    auto p = senri::cli::parse_command_line(plain, 4);
    REQUIRE(p.has_value());
    const auto o = senri::annotate::paired_options(p->config.vcf);
    CHECK(o.min_depth.text == "5");
    CHECK(o.min_var_depth.text == "3");
    CHECK(o.min_mq.value == 0);
    CHECK(o.max_pvalue.text == "0.05");
    CHECK_FALSE(o.somatic_only);

    const std::vector<std::string_view> set = {"-b",
                                               "t.bam|n.bam",
                                               "-G",
                                               "ref.fa",
                                               "-R",
                                               "chr1:1-10",
                                               "--filter-min-depth",
                                               "7",
                                               "--filter-max-pvalue=0.01",
                                               "--vcf-somatic-only"};
    auto q = senri::cli::parse_command_line(set, 4);
    REQUIRE(q.has_value());
    const auto o2 = senri::annotate::paired_options(q->config.vcf);
    CHECK(o2.min_depth.text == "7");
    CHECK(o2.min_var_depth.text == "3");
    CHECK(o2.max_pvalue.value == 0.01);
    CHECK(o2.somatic_only);
}
