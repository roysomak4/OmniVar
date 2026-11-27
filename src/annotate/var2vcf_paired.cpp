// Paired (tumour/normal) VCF: native replacement for
//   testsomatic.R | var2vcf_paired.pl -A
// (the upstream 1.8.2 R and Perl scripts).
//
// Input rows are the Somatic 55-column TSV (SomaticOutputVariant without --fisher). The R
// step adds three Fisher tests per row (tumour strands, normal strands, tumour vs normal, the
// last as the smaller of the one-sided p-values), stats::FisherFlavor::r; the Perl step computes
// the tumour and normal FILTER lists, STATUS, QUAL, GT and both FORMAT columns.
//
// Deliberate differences from the Perl script, as in the single-sample port (var2vcf.hpp):
// records are never dropped (no -S), all alleles at a position are written (-A), thresholds of 0
// are honoured, contigs come from the BAM header, and every reorder() contig is kept. The normal
// column is named after the second sample name (-N "T|N" or the BAM name) instead of
// "<tumour>-match" when one is known.

#include <algorithm>
#include <cmath>

#include "annotate/perl_fields.hpp"
#include "annotate/var2vcf.hpp"
#include "stats/fisher.hpp"
#include "util/num_format.hpp"

namespace senri::annotate {

using namespace detail;

namespace {

// Column indexes of the Somatic TSV row (SomaticOutputVariant, 55 columns). In the Perl
// script the normal columns are shifted by two (R inserts SBF/ODDRATIO after column 25) and the
// trailing ones by four.
enum Col : std::size_t {
    kSample = 0,
    kGene,
    kChr,
    kStart,
    kEnd,
    kRef,
    kAlt,
    kDp1,
    kVd1,
    kRfwd1,
    kRrev1,
    kVfwd1,
    kVrev1,
    kGt1,
    kAf1,
    kBias1,
    kPmean1,
    kPstd1,
    kQual1,
    kQstd1,
    kMapq1,
    kSn1,
    kHiaf1,
    kAdjaf1,
    kNm1,
    kDp2,
    kVd2,
    kRfwd2,
    kRrev2,
    kVfwd2,
    kVrev2,
    kGt2,
    kAf2,
    kBias2,
    kPmean2,
    kPstd2,
    kQual2,
    kQstd2,
    kMapq2,
    kSn2,
    kHiaf2,
    kAdjaf2,
    kNm2,
    kShift3,
    kMsi,
    kMsilen,
    kLseq,
    kRseq,
    kSeg,
    kStatus,
    kType,
    kDuprate1,
    kSv1,
    kDuprate2,
    kSv2
};

Threshold or_default(const Threshold& user, double v, const char* text) {
    return user.explicit_ ? user : Threshold{v, text, false};
}

// One sample's strand test (testsomatic.R) and the odds ratio as var2vcf_paired.pl rewrites it:
// Inf -> 0, values in (0, 1) -> sprintf("%.2f", 1/x).
struct Strand {
    double sbf;
    std::string sbf_text;
    double oddratio;
    std::string oddratio_text;
};

Strand strand_test(const std::vector<std::string>& f, std::size_t rfwd, std::size_t rrev, std::size_t vfwd,
                   std::size_t vrev) {
    const stats::FisherExact fisher(static_cast<int>(num(f, rfwd)), static_cast<int>(num(f, rrev)),
                                    static_cast<int>(num(f, vfwd)), static_cast<int>(num(f, vrev)),
                                    stats::FisherFlavor::r);
    Strand s;
    s.sbf = fisher.p_value();
    s.sbf_text = util::rounded_value_to_print(s.sbf, 5);
    s.oddratio = fisher.odds_ratio();
    if (std::isinf(s.oddratio)) {
        s.oddratio = 0;
        s.oddratio_text = "0";
    } else if (s.oddratio < 1 && s.oddratio > 0) {
        s.oddratio_text = util::printf_fixed(1 / s.oddratio, 2);
        s.oddratio = util::parse_double(s.oddratio_text).value_or(0);
    } else {
        s.oddratio_text = util::rounded_value_to_print(s.oddratio, 5);
    }
    return s;
}

// The genotype rule of var2vcf_paired.pl:220
std::string genotype(double af, const Var2VcfOptions& o) {
    return (1 - af < o.hom_freq.value) ? "1/1" : af >= 0.5 ? "1/0" : af >= o.min_af.value ? "0/1" : "0/0";
}

// $bias =~ s/;/,/; "0" -> "0,0"
std::string bias_text(std::string b) {
    if (auto p = b.find(';'); p != std::string::npos) {
        b[p] = ',';
    }
    return b == "0" ? "0,0" : b;
}

// Perl int() of a double, printed: integers, or Perl's "Inf" / "-Inf" / "NaN".
std::string perl_int(double v) {
    if (std::isnan(v) || std::isinf(v)) {
        return util::printf_fixed(v, 0);
    }
    return integer_text(std::trunc(v));
}

}  // namespace

Var2VcfOptions paired_options(const Var2VcfOptions& u) {
    Var2VcfOptions o = u;
    o.min_depth = or_default(u.min_depth, 5, "5");
    o.min_var_depth = or_default(u.min_var_depth, 3, "3");
    o.min_af = or_default(u.min_af, 0.02, "0.02");
    o.min_pmean = or_default(u.min_pmean, 8, "8");
    o.min_qual = or_default(u.min_qual, 22.5, "22.5");
    o.min_mq = or_default(u.min_mq, 0, "0");
    o.hom_freq = or_default(u.hom_freq, 0.2, "0.2");
    o.min_sn = or_default(u.min_sn, 1.5, "1.5");
    o.max_pvalue = or_default(u.max_pvalue, 0.05, "0.05");
    o.msi = or_default(u.msi, 12, "12");
    o.max_nm = or_default(u.max_nm, 5.25, "5.25");
    o.cluster_bp = or_default(u.cluster_bp, 0, "0");
    return o;
}

std::string Var2Vcf::paired_header(const VcfHeaderInfo& info) const {
    // var2vcf_paired.pl:51-115
    const Var2VcfOptions& o = opt_;
    std::string h = "##fileformat=VCFv4.2\n";
    h += "##source=" + info.source + "\n";
    if (!info.reference.empty()) {
        h += "##reference=" + info.reference + "\n";
    }
    if (info.contigs.empty()) {
        for (const std::string& name : chromosome_order()) {
            h += "##contig=<ID=" + name + ">\n";
        }
    }
    for (const auto& [name, length] : info.contigs) {
        h += "##contig=<ID=" + name + ",length=" + std::to_string(length) + ">\n";
    }
    if (!info.command_line.empty()) {
        h += "##senriCommandLine=\"" + info.command_line + "\"\n";
    }
    h +=
        R"VCF(##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MSI unit repeat length in bp">
##INFO=<ID=SSF,Number=1,Type=Float,Description="P-value">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Odds ratio">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=STATUS,Number=1,Type=String,Description="Somatic or germline status">
##INFO=<ID=P0.01Likely,Number=0,Type=Flag,Description="Likely candidate but p-value > 0.01/5**vd2 (means the evidence in tumor sample might be weak, e.g. small diff in AF)">
##INFO=<ID=InDelLikely,Number=0,Type=Flag,Description="Likely indels more than 2bp are not considered somatic (weak evidence of presence in normal samples)">
)VCF";
    h += "##FILTER=<ID=q" + o.min_qual.text + ",Description=\"Mean Base Quality Below " + o.min_qual.text + "\">\n";
    h += "##FILTER=<ID=Q" + o.min_mq.text + ",Description=\"Mean Mapping Quality Below " + o.min_mq.text + "\">\n";
    h += "##FILTER=<ID=p" + o.min_pmean.text + ",Description=\"Mean Position in Reads Less than " + o.min_pmean.text +
         "\">\n";
    h += "##FILTER=<ID=SN" + o.min_sn.text + ",Description=\"Signal to Noise Less than " + o.min_sn.text + "\">\n";
    h += "##FILTER=<ID=Bias,Description=\"Strand Bias\">\n";
    h += "##FILTER=<ID=pSTD,Description=\"Position in Reads has STD of 0\">\n";
    h += "##FILTER=<ID=MAF0.05,Description=\"Matched sample has AF > 0.05, thus not somatic\">\n";
    h += "##FILTER=<ID=d" + o.min_depth.text + ",Description=\"Total Depth < " + o.min_depth.text + "\">\n";
    h += "##FILTER=<ID=v" + o.min_var_depth.text + ",Description=\"Var Depth < " + o.min_var_depth.text + "\">\n";
    h += "##FILTER=<ID=f" + o.min_af.text + ",Description=\"Allele frequency < " + o.min_af.text + "\">\n";
    h += "##FILTER=<ID=P" + o.max_pvalue.text + ",Description=\"Not significant with p-value > " + o.max_pvalue.text +
         "\">\n";
    h += "##FILTER=<ID=DIFF0.2,Description=\"Non-somatic or LOH and allele frequency difference < 0.2\">\n";
    h += "##FILTER=<ID=MSI" + o.msi.text + ",Description=\"Variant in MSI region with " + o.msi.text +
         " non-monomer MSI or 12 monomer MSI\">\n";
    h += "##FILTER=<ID=NM" + o.max_nm.text + ",Description=\"Mean mismatches in reads >= " + o.max_nm.text +
         ", thus likely false positive\">\n";
    h += "##FILTER=<ID=InGap,Description=\"The somatic variant is in the deletion gap, thus likely false positive\">\n";
    h += "##FILTER=<ID=InIns,Description=\"The somatic variant is adjacent to an insertion variant\">\n";
    h += "##FILTER=<ID=Cluster" + o.cluster_bp.text + "bp,Description=\"Two somatic variants are within " +
         o.cluster_bp.text + " bp\">\n";
    h += "##FILTER=<ID=LongAT,Description=\"The somatic variant is flanked by long A/T (>=14)\">\n";
    for (const auto& [id, description] : o.reject_filters) {
        h += "##FILTER=<ID=" + id + ",Description=\"" + description + "\">\n";
    }
    h += R"VCF(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##FORMAT=<ID=BIAS,Number=2,Type=String,Description="Strand Bias Info">
##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">
##FORMAT=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##FORMAT=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##FORMAT=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##FORMAT=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Mean Mapping Quality">
##FORMAT=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##FORMAT=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##FORMAT=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
)VCF";
    const std::string samplem = info.samplem.empty() ? info.sample + "-match" : info.samplem;
    h += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + info.sample + "\t" + samplem + "\n";
    return h;
}

std::string Var2Vcf::paired_records(const std::string& sample_name) const {
    std::string sample_nowhitespace = sample_name;
    for (char& c : sample_nowhitespace) {
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v') {
            c = '_';
        }
    }
    const Var2VcfOptions& o = opt_;
    const std::string tag_d = "d" + o.min_depth.text;
    const std::string tag_v = "v" + o.min_var_depth.text;
    const std::string tag_f = "f" + o.min_af.text;
    const std::string tag_p = "p" + o.min_pmean.text;
    const std::string tag_q = "q" + o.min_qual.text;
    const std::string tag_Q = "Q" + o.min_mq.text;
    const std::string tag_sn = "SN" + o.min_sn.text;
    const std::string tag_nm = "NM" + o.max_nm.text;
    const std::string tag_msi = "MSI" + o.msi.text;
    const std::string tag_P = "P" + o.max_pvalue.text;
    const std::string tag_cluster = "Cluster" + o.cluster_bp.text + "bp";
    const auto has = [](const std::vector<std::string>& v, const std::string& t) {
        return std::find(v.begin(), v.end(), t) != v.end();
    };

    std::string out;
    for (const std::string& chr : chromosome_order()) {
        const Chromosome& chromosome = chrs_.at(chr);
        double pvs = 0;  // previous PASS SNV start
        std::string pinfo1;
        std::string pfilter;
        std::string pinfo2;
        // -p: rejected records that follow the pending record; written right after it.
        std::string prejected;
        const auto flush = [&] {
            if (!pinfo1.empty()) {
                out += pinfo1 + "\t" + pfilter + "\t" + pinfo2 + "\n";
            }
            out += prejected;
            prejected.clear();
        };
        for (const auto& [start_key, rows_at_start] : chromosome.by_start) {
            // sort { $b->[14] <=> $a->[14] }: tumour AF, descending (stable).
            std::vector<const Row*> tmp;
            for (const Row& r : rows_at_start) {
                tmp.push_back(&r);
            }
            std::stable_sort(tmp.begin(), tmp.end(),
                             [](const Row* a, const Row* b) { return num(a->f, kAf1) > num(b->f, kAf1); });
            std::vector<std::string> seen;
            std::vector<std::string> seen_rejected;
            std::size_t kept = 0;  // rows the script reads (-p rejected rows are not counted)
            for (const Row* row : tmp) {
                const bool rejected = !row->reject.empty();
                if (!rejected) {
                    if (o.top_allele_only && kept >= 1) {
                        continue;
                    }
                    ++kept;
                }
                const std::vector<std::string>& f = row->f;
                const std::string& ref = str(f, kRef);
                const std::string& alt = str(f, kAlt);
                if (!perl_true(ref)) {
                    continue;
                }
                const std::string key =
                    str(f, kChr) + "-" + str(f, kStart) + "-" + str(f, kEnd) + "-" + ref + "-" + alt;
                std::vector<std::string>& seen_here = rejected ? seen_rejected : seen;
                if (has(seen_here, key)) {
                    continue;
                }
                seen_here.push_back(key);
                std::string type = str(f, kType);
                if (!perl_true(type)) {
                    type = "REF";
                }
                const std::string& status = str(f, kStatus);
                const double start = num(f, kStart);
                const double end = num(f, kEnd);
                const double dp1 = num(f, kDp1), vd1 = num(f, kVd1), af1 = num(f, kAf1);
                const double dp2 = num(f, kDp2), vd2 = num(f, kVd2), af2 = num(f, kAf2);
                const double rd1 = num(f, kRfwd1) + num(f, kRrev1);
                const double rd2 = num(f, kRfwd2) + num(f, kRrev2);

                // testsomatic.R
                const Strand s1 = strand_test(f, kRfwd1, kRrev1, kVfwd1, kVrev1);
                const Strand s2 = strand_test(f, kRfwd2, kRrev2, kVfwd2, kVrev2);
                const int tref = std::max(0, static_cast<int>(dp1 - vd1));
                const int rref = std::max(0, static_cast<int>(dp2 - vd2));
                const stats::FisherExact tn(static_cast<int>(vd1), tref, static_cast<int>(vd2), rref,
                                            stats::FisherFlavor::r);
                const double pv_greater = tn.p_value_greater();
                const double pv_less = tn.p_value_less();
                // `if (h$p.value < pv)` takes "less"; both estimates are the same conditional MLE.
                const double pvalue = pv_less < pv_greater ? pv_less : pv_greater;
                const std::string pvalue_text = util::rounded_value_to_print(pvalue, 5);
                const double sor = tn.odds_ratio();
                const std::string sor_text = std::isinf(sor) ? "0" : util::rounded_value_to_print(sor, 5);

                std::vector<std::string> filters;
                std::vector<std::string> filters2;
                const bool strong_waiver = status == "StrongSomatic" && pvalue < 0.15 && af1 * vd1 >= 0.5;
                if (dp1 < o.min_depth.value && !strong_waiver) {
                    filters.push_back(tag_d);
                }
                if (vd1 < o.min_var_depth.value && !strong_waiver) {
                    filters.push_back(tag_v);
                }
                if (dp2 < o.min_depth.value) {
                    filters2.push_back(tag_d);
                }
                if (vd2 < o.min_var_depth.value) {
                    filters2.push_back(tag_v);
                }
                if (af1 < o.min_af.value) {
                    filters.push_back(tag_f);
                }
                if (num(f, kPmean1) < o.min_pmean.value) {
                    filters.push_back(tag_p);
                }
                if (o.pstd && num(f, kPstd1) == 0 && vd1 < o.min_depth.value) {
                    filters.push_back("pSTD");
                }
                if (num(f, kQual1) < o.min_qual.value) {
                    filters.push_back(tag_q);
                }
                const double mapq1 = num(f, kMapq1);
                if (mapq1 < o.min_mq.value) {
                    filters.push_back(tag_Q);
                }
                if (mapq1 < 10 && type == "SNV") {  // SNV somatic in low mapping quality regions
                    filters.push_back(tag_Q);
                }
                if (num(f, kSn1) < o.min_sn.value) {
                    filters.push_back(tag_sn);
                }
                if (num(f, kNm1) >= o.max_nm.value) {
                    filters.push_back(tag_nm);
                }
                const bool short_span = end - start < 100;
                if (str(f, kBias1) == "2;1" && s1.sbf < 0.01 && (s1.oddratio > 5 || s1.oddratio == 0) && short_span) {
                    filters.push_back("Bias");
                }
                if (af2 < o.min_af.value) {
                    filters2.push_back(tag_f);
                }
                if (num(f, kPmean2) < o.min_pmean.value) {
                    filters2.push_back(tag_p);
                }
                if (o.pstd && num(f, kPstd2) == 0 && vd2 < o.min_depth.value) {
                    filters2.push_back("pSTD");
                }
                if (num(f, kQual2) < o.min_qual.value) {
                    filters2.push_back(tag_q);
                }
                const double mapq2 = num(f, kMapq2);
                if (mapq2 < o.min_mq.value) {
                    filters2.push_back(tag_Q);
                }
                if (num(f, kSn2) < o.min_sn.value) {
                    filters2.push_back(tag_sn);
                }
                if (num(f, kNm2) >= o.max_nm.value) {
                    filters2.push_back(tag_nm);
                }
                if (!has(filters, "Bias") && str(f, kBias2) == "2;1" && s2.sbf < 0.01 &&
                    (s2.oddratio > 5 || s2.oddratio == 0) && short_span) {
                    filters.push_back("Bias");
                }
                // Require stringent statistics in regions with MSI
                const double msi = num(f, kMsi);
                const double msilen = num(f, kMsilen);
                if ((msi > o.msi.value && msilen > 1) || (msi > 12 && msilen == 1)) {
                    if (!(status == "StrongSomatic" && pvalue < 0.0005)) {
                        filters.push_back(tag_msi);
                    }
                }
                const double indel_len = std::abs(static_cast<double>(ref.size()) - static_cast<double>(alt.size()));
                // !grep(/^MSI<x>/, @filters)
                const bool msi_tagged = std::any_of(filters.begin(), filters.end(),
                                                    [&](const std::string& t) { return t.starts_with(tag_msi); });
                if (indel_len == msilen && !msi_tagged) {
                    if ((msi > o.msi.value && msilen > 1 && af1 < 0.35 && af2 < 0.35) ||
                        (msi > 12 && msilen == 1 && af1 < 0.35 && af2 < 0.35)) {
                        filters.push_back(tag_msi);
                    }
                }
                bool p_likely = false;
                bool indel_likely = false;
                if (o.somatic_only) {
                    if (pvalue > o.max_pvalue.value) {
                        if (!(status == "StrongSomatic" &&
                              ((pvalue < 0.25 && af1 > 0.1) || (pvalue < 0.5 && af1 > 0.20) ||
                               (pvalue < 0.15 && af1 > 0.05)))) {
                            filters.push_back(tag_P);
                        }
                    } else if (status.contains("LikelySomatic") && pvalue > 0.05 / std::pow(5.0, vd2)) {
                        p_likely = true;  // Increase the stringency for LikelySomatic
                    } else if (status.contains("Likely") && type != "SNV") {
                        indel_likely = !(ref.size() <= 2 && alt.size() <= 2);
                    }
                }
                std::string filter;
                for (const std::string& t : filters) {
                    filter += (filter.empty() ? "" : ";") + t;
                }
                if (filter.empty()) {
                    filter = "PASS";
                }
                // Unless somatic only (-M), good variants in the germline are reported as well.
                if (!o.somatic_only && filter != "PASS" && filters2.empty()) {
                    filter = "PASS";
                }
                if (rejected) {
                    filter = row->reject + (filter == "PASS" ? "" : ";" + filter);
                }
                const std::string gt = genotype(af1, o);
                const std::string gtm = genotype(af2, o);
                const std::string bias1 = bias_text(str(f, kBias1));
                const std::string bias2 = bias_text(str(f, kBias2));
                const std::string mapq1_text = util::printf_fixed(mapq1, 0);
                const std::string mapq2_text = util::printf_fixed(mapq2, 0);
                const std::string qual = vd1 > vd2 ? perl_int(std::log(vd1) / std::log(2.0) * num(f, kQual1))
                                                   : perl_int(std::log(vd2) / std::log(2.0) * num(f, kQual2));
                if (!rejected && pfilter == "PASS" && pinfo2.contains("Somatic") && pinfo2.contains("TYPE=SNV") &&
                    filter == "PASS" && status.contains("Somatic") && type == "SNV" &&
                    start - pvs < o.cluster_bp.value) {
                    pfilter = tag_cluster;
                    filter = tag_cluster;
                }
                if (!rejected) {
                    flush();
                }
                const auto sample_column =
                    [&](const std::string& g, std::size_t dp, std::size_t vd, std::size_t vfwd, std::size_t vrev,
                        std::size_t rfwd, std::size_t rrev, double rd, std::size_t af, const std::string& bias,
                        std::size_t pmean, std::size_t pstd, std::size_t qual_i, std::size_t qstd, const Strand& st,
                        const std::string& mapq, std::size_t sn, std::size_t hiaf, std::size_t adjaf, std::size_t nm) {
                        return g + ":" + str(f, dp) + ":" + str(f, vd) + ":" + str(f, vfwd) + "," + str(f, vrev) + ":" +
                               str(f, rfwd) + "," + str(f, rrev) + ":" + integer_text(rd) + "," + str(f, vd) + ":" +
                               str(f, af) + ":" + bias + ":" + str(f, pmean) + ":" + str(f, pstd) + ":" +
                               str(f, qual_i) + ":" + str(f, qstd) + ":" + st.sbf_text + ":" + st.oddratio_text + ":" +
                               mapq + ":" + str(f, sn) + ":" + str(f, hiaf) + ":" + str(f, adjaf) + ":" + str(f, nm);
                    };
                std::string info1 = chr + "\t" + str(f, kStart) + "\t.\t" + ref + "\t" + alt + "\t" + qual;
                std::string info2 =
                    "STATUS=" + status + ";SAMPLE=" + sample_nowhitespace + ";TYPE=" + type + ";DP=" + str(f, kDp1) +
                    ";VD=" + str(f, kVd1) + ";AF=" + str(f, kAf1) + ";SHIFT3=" + str(f, kShift3) +
                    ";MSI=" + str(f, kMsi) + ";MSILEN=" + str(f, kMsilen) + ";SSF=" + pvalue_text + ";SOR=" + sor_text +
                    ";LSEQ=" + str(f, kLseq) + ";RSEQ=" + str(f, kRseq) + (p_likely ? ";P0.01Likely" : "") +
                    (indel_likely ? ";InDelLikely" : "") +
                    "\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t" +
                    sample_column(gt, kDp1, kVd1, kVfwd1, kVrev1, kRfwd1, kRrev1, rd1, kAf1, bias1, kPmean1, kPstd1,
                                  kQual1, kQstd1, s1, mapq1_text, kSn1, kHiaf1, kAdjaf1, kNm1) +
                    "\t" +
                    sample_column(gtm, kDp2, kVd2, kVfwd2, kVrev2, kRfwd2, kRrev2, rd2, kAf2, bias2, kPmean2, kPstd2,
                                  kQual2, kQstd2, s2, mapq2_text, kSn2, kHiaf2, kAdjaf2, kNm2);
                if (rejected) {
                    std::string line = info1 + "\t" + filter + "\t" + info2 + "\n";
                    if (pinfo1.empty()) {
                        out += line;
                    } else {
                        prejected += line;
                    }
                    continue;
                }
                pinfo1 = std::move(info1);
                pfilter = filter;
                pinfo2 = std::move(info2);
                if (type == "SNV" && filter == "PASS") {
                    pvs = start;
                }
            }
        }
        flush();
    }
    return out;
}

}  // namespace senri::annotate
