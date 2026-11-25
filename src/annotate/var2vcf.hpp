#pragma once

// Single-sample VCF from the Simple-mode rows: native replacement for
//   teststrandbias.R | var2vcf_valid.pl -A
// (tests/scripts/teststrandbias.R, tests/scripts/var2vcf_valid.pl, upstream 1.8.2 scripts).
//
// Input rows are the 36-column TSV rows (optionally with the CRISPR column), exactly as
// SimpleOutputVariant prints them without --fisher. The R step adds the strand-bias Fisher test
// on columns 10-13 (stats::FisherFlavor::r); the Perl step computes FILTER, QUAL, GT and AD.
//
// Values are written as the TSV formats them. The R + Perl chain prints the same numbers
// as R re-formats them (write.table, 15 significant digits), so the two agree numerically.
//
// Deliberate differences from the Perl script (all keep records the script would keep):
// - FILTER annotation never removes records: no -S. All alleles at a position are written by
//   default (-A); `top_allele_only` restores the script's default.
// - Thresholds of 0 are honoured (Perl's `$opt_x ? $opt_x : default` treats 0 as "unset").
// - Chromosome order: as reorder() (numeric, then X, Y, MT/chrM, then the rest); "the rest" is in
//   order of first appearance instead of Perl hash order, and contigs reorder() would drop
//   (chrX when X also exists, etc.) are kept.
// - The header lists contigs from the BAM header and declares CRISPR.

#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace omnivar::annotate {

// A threshold: the value, and its text as given (Perl writes the option text into the tag name,
// e.g. "f0.02", "q22.5").
struct Threshold {
    double value;
    std::string text;
    bool explicit_ = false;  // given on the command line (the paired script has other defaults)
};

// var2vcf_valid.pl options and their defaults (var2vcf_valid.pl:19-31).
struct Var2VcfOptions {
    Threshold min_depth{3, "3"};        // -d  d<x>
    Threshold min_var_depth{2, "2"};    // -v  v<x>
    Threshold min_af{0.02, "0.02"};     // -f  f<x>; also the 0/1 genotype limit
    Threshold min_pmean{8, "8"};        // -p  p<x>
    Threshold min_qual{22.5, "22.5"};   // -q  q<x>
    Threshold min_mq{10, "10"};         // -Q  Q<x>
    Threshold hom_freq{0.2, "0.2"};     // -F  1-AF below this -> 1/1
    Threshold min_sn{1.5, "1.5"};       // -o  SN<x>
    Threshold max_nm{5.25, "5.25"};     // -m  NM<x>
    Threshold msi{12, "12"};            // -I  MSI<x>
    Threshold cluster_bp{0, "0"};       // -c  Cluster<x>bp (0: off)
    bool pstd = true;                   // -P  pSTD
    Threshold min_split_reads{1, "1"};  // -T  SVs need this many split reads
    bool top_allele_only = false;       // not -A
    // var2vcf_paired.pl only
    Threshold max_pvalue{0.05, "0.05"};  // -P  P<x> (with -M)
    bool somatic_only = false;           // -M
    // -p (rejected variants): ##FILTER lines (ID, description) for the caller's REJ_* reasons.
    std::vector<std::pair<std::string, std::string>> reject_filters;
};

struct VcfHeaderInfo {
    std::string source;     // ##source=
    std::string sample;     // sample column name (paired: the tumour)
    std::string samplem;    // paired: the normal sample column (empty: "<sample>-match")
    std::string reference;  // ##reference= (empty: omitted)
    // ##contig lines, in this order; when empty, the contigs that have records (no length).
    std::vector<std::pair<std::string, long long>> contigs;
    std::string command_line;  // ##<name>CommandLine= (empty: omitted)
};

// var2vcf_paired.pl's defaults (var2vcf_paired.pl:19-31) for every threshold the user did not
// set: d5, v3, f0.02, p8, q22.5, Q0, F0.2, SN1.5, P0.05, MSI12, NM5.25, Cluster0bp.
[[nodiscard]] Var2VcfOptions paired_options(const Var2VcfOptions& user);

class Var2Vcf {
public:
    // `paired`: rows are the Somatic 55-column TSV, converted as
    // testsomatic.R | var2vcf_paired.pl -A (options from paired_options()).
    explicit Var2Vcf(Var2VcfOptions options, bool paired = false)
        : opt_(paired ? paired_options(options) : std::move(options)), paired_(paired) {}

    // One TSV row (no trailing newline). Rows are grouped by chromosome and start, keeping
    // the order they were added in within a start, as the Perl script reads them.
    // `reject` (-p): the caller rejected this variant; the record's FILTER starts with this ID.
    // Rejected records do not take part in the script's links between records (cluster tag,
    // top allele, duplicates), so every other record is written as without them. Single sample only.
    void add_row(std::string_view line, std::string_view reject = {});

    // VCF header (##fileformat ... #CHROM line).
    [[nodiscard]] std::string header(const VcfHeaderInfo& info) const;

    // All records, in var2vcf order. `sample_name` is the SAMPLE= INFO value (whitespace -> '_').
    [[nodiscard]] std::string records(const std::string& sample_name) const;

private:
    struct Row {
        std::vector<std::string> f;
        std::string reject;  // REJ_* ID, or empty
    };
    struct Chromosome {
        std::map<long long, std::vector<Row>> by_start;
    };
    [[nodiscard]] std::vector<std::string> chromosome_order() const;
    [[nodiscard]] std::string paired_header(const VcfHeaderInfo& info) const;
    [[nodiscard]] std::string paired_records(const std::string& sample_name) const;

    Var2VcfOptions opt_;
    bool paired_ = false;
    std::map<std::string, Chromosome> chrs_;
    std::vector<std::string> first_seen_;
};

}  // namespace omnivar::annotate
