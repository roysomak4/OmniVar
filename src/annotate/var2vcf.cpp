#include "annotate/var2vcf.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "annotate/perl_fields.hpp"
#include "stats/fisher.hpp"
#include "util/num_format.hpp"
#include "util/number_parse.hpp"

namespace omnivar::annotate {

using namespace detail;

namespace {

// Column indexes of the Simple TSV row (SimpleOutputVariant, 36 columns + optional CRISPR).
// In the Perl script these are shifted by two after column 20, where R inserts SBF and ODDRATIO.
enum Col : std::size_t {
    kSample = 0,
    kGene,
    kChr,
    kStart,
    kEnd,
    kRef,
    kAlt,
    kDp,
    kVd,
    kRfwd,
    kRrev,
    kVfwd,
    kVrev,
    kGenotype,
    kAf,
    kBias,
    kPmean,
    kPstd,
    kQual,
    kQstd,
    kMapq,
    kSn,
    kHiaf,
    kAdjaf,
    kShift3,
    kMsi,
    kMsilen,
    kNm,
    kHicnt,
    kHicov,
    kLseq,
    kRseq,
    kSeg,
    kType,
    kGamp /* Duprate */,
    kTamp /* SV_info */,
    kNcamp /* CRISPR */,
    kAmpflag
};

// reorder() (var2vcf_valid.pl:254): the digits of a name containing a digit and no '_'.
double chromosome_number(const std::string& c) {
    std::string digits;
    for (char ch : c) {
        if (ch >= '0' && ch <= '9') {
            digits += ch;
        }
    }
    return util::parse_double(digits).value_or(0);
}

bool has_digit(const std::string& c) {
    return std::any_of(c.begin(), c.end(), [](char ch) { return ch >= '0' && ch <= '9'; });
}

}  // namespace

void Var2Vcf::add_row(std::string_view line, std::string_view reject) {
    // var2vcf_valid.pl:35-42
    if (line.find("R_HOME") != std::string_view::npos) {
        return;
    }
    Row row{split_tabs(line), std::string(reject)};
    const std::string chr = str(row.f, kChr);
    const auto start = static_cast<long long>(num(row.f, kStart));
    auto [it, inserted] = chrs_.try_emplace(chr);
    if (inserted) {
        first_seen_.push_back(chr);
    }
    it->second.by_start[start].push_back(std::move(row));
}

std::vector<std::string> Var2Vcf::chromosome_order() const {
    std::vector<std::pair<double, std::string>> numeric;
    std::vector<std::string> rest;
    std::vector<std::string> special;
    for (const std::string& c : first_seen_) {
        if (has_digit(c) && !c.contains('_')) {
            numeric.emplace_back(chromosome_number(c), c);
        } else if (c == "X" || c == "chrX" || c == "Y" || c == "chrY" || c == "MT" || c == "chrM") {
            special.push_back(c);
        } else {
            rest.push_back(c);
        }
    }
    std::stable_sort(numeric.begin(), numeric.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::string> out;
    for (auto& [n, c] : numeric) {
        out.push_back(c);
    }
    for (const char* name : {"X", "chrX", "Y", "chrY", "MT", "chrM"}) {
        if (std::find(special.begin(), special.end(), name) != special.end()) {
            out.emplace_back(name);
        }
    }
    out.insert(out.end(), rest.begin(), rest.end());
    return out;
}

std::string Var2Vcf::header(const VcfHeaderInfo& info) const {
    if (paired_) {
        return paired_header(info);
    }
    // var2vcf_valid.pl:46-116
    std::string h = "##fileformat=VCFv4.2\n";
    h += "##source=" + info.source + "\n";
    if (!info.reference.empty()) {
        h += "##reference=" + info.reference + "\n";
    }
    if (info.contigs.empty()) {
        // No BAM header (e.g. --dev-tsv-to-vcf): declare the contigs that have records.
        for (const std::string& name : chromosome_order()) {
            h += "##contig=<ID=" + name + ">\n";
        }
    }
    for (const auto& [name, length] : info.contigs) {
        h += "##contig=<ID=" + name + ",length=" + std::to_string(length) + ">\n";
    }
    if (!info.command_line.empty()) {
        h += "##omnivarCommandLine=\"" + info.command_line + "\"\n";
    }
    h +=
        R"VCF(##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MicroSatellite unit length in bp">
##INFO=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=GDAMP,Number=1,Type=Integer,Description="No. of amplicons supporting variant">
##INFO=<ID=TLAMP,Number=1,Type=Integer,Description="Total of amplicons covering variant">
##INFO=<ID=NCAMP,Number=1,Type=Integer,Description="No. of amplicons don't work">
##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description="Top variant in amplicons don't match">
##INFO=<ID=HICNT,Number=1,Type=Integer,Description="High quality variant reads">
##INFO=<ID=HICOV,Number=1,Type=Integer,Description="High quality total reads">
##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="No. of split reads supporting SV">
##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="No. of pairs supporting SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type: INV DUP DEL INS FUS">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="The length of SV in bp">
##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="Duplication rate in fraction">
##INFO=<ID=CRISPR,Number=1,Type=Integer,Description="Variant is close to the CRISPR cutting site (-J)">
)VCF";
    const Var2VcfOptions& o = opt_;
    h += "##FILTER=<ID=q" + o.min_qual.text + ",Description=\"Mean Base Quality Below " + o.min_qual.text + "\">\n";
    h += "##FILTER=<ID=Q" + o.min_mq.text + ",Description=\"Mean Mapping Quality Below " + o.min_mq.text + "\">\n";
    h += "##FILTER=<ID=p" + o.min_pmean.text + ",Description=\"Mean Position in Reads Less than " + o.min_pmean.text +
         "\">\n";
    h += "##FILTER=<ID=SN" + o.min_sn.text + ",Description=\"Signal to Noise Less than " + o.min_sn.text + "\">\n";
    h += "##FILTER=<ID=Bias,Description=\"Strand Bias\">\n";
    h += "##FILTER=<ID=pSTD,Description=\"Position in Reads has STD of 0\">\n";
    h += "##FILTER=<ID=d" + o.min_depth.text + ",Description=\"Total Depth < " + o.min_depth.text + "\">\n";
    h += "##FILTER=<ID=v" + o.min_var_depth.text + ",Description=\"Var Depth < " + o.min_var_depth.text + "\">\n";
    h += "##FILTER=<ID=f" + o.min_af.text + ",Description=\"Allele frequency < " + o.min_af.text + "\">\n";
    h += "##FILTER=<ID=MSI" + o.msi.text + ",Description=\"Variant in MSI region with " + o.msi.text +
         " non-monomer MSI or 13 monomer MSI\">\n";
    h += "##FILTER=<ID=NM" + o.max_nm.text + ",Description=\"Mean mismatches in reads >= " + o.max_nm.text +
         ", thus likely false positive\">\n";
    h += "##FILTER=<ID=InGap,Description=\"The variant is in the deletion gap, thus likely false positive\">\n";
    h += "##FILTER=<ID=InIns,Description=\"The variant is adjacent to an insertion variant\">\n";
    h += "##FILTER=<ID=Cluster" + o.cluster_bp.text + "bp,Description=\"Two variants are within " + o.cluster_bp.text +
         " bp\">\n";
    h += "##FILTER=<ID=LongMSI,Description=\"The somatic variant is flanked by long A/T (>=14)\">\n";
    h += "##FILTER=<ID=AMPBIAS,Description=\"Indicate the variant has amplicon bias.\">\n";
    for (const auto& [id, description] : o.reject_filters) {
        h += "##FILTER=<ID=" + id + ",Description=\"" + description + "\">\n";
    }
    h += R"VCF(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
)VCF";
    h += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + info.sample + "\n";
    return h;
}

std::string Var2Vcf::records(const std::string& sample_name) const {
    if (paired_) {
        return paired_records(sample_name);
    }
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
    const std::string tag_cluster = "Cluster" + o.cluster_bp.text + "bp";

    std::string out;
    for (const std::string& chr : chromosome_order()) {
        const Chromosome& chromosome = chrs_.at(chr);
        double pvs = 0;  // previous PASS SNV start
        // The previous record is printed only when the next one is processed, so that a following
        // SNV can mark it as clustered (var2vcf_valid.pl:181-194).
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
            // sort { $b->[14] <=> $a->[14] }: AF, descending; Perl's merge sort is stable.
            std::vector<const Row*> tmp;
            for (const Row& r : rows_at_start) {
                tmp.push_back(&r);
            }
            std::stable_sort(tmp.begin(), tmp.end(),
                             [](const Row* a, const Row* b) { return num(a->f, kAf) > num(b->f, kAf); });
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
                if (!perl_true(ref)) {
                    continue;
                }
                std::string alt = str(f, kAlt);
                const std::string key =
                    str(f, kChr) + "-" + str(f, kStart) + "-" + str(f, kEnd) + "-" + ref + "-" + alt;
                std::vector<std::string>& seen_here = rejected ? seen_rejected : seen;
                if (std::find(seen_here.begin(), seen_here.end(), key) != seen_here.end()) {
                    continue;
                }
                seen_here.push_back(key);
                std::string type = str(f, kType);
                if (!perl_true(type)) {
                    type = "REF";
                }
                const bool isamp = f.size() > kAmpflag;
                const double start = num(f, kStart);
                const double end = num(f, kEnd);
                const double dp = num(f, kDp);
                const double rfwd = num(f, kRfwd);
                const double rrev = num(f, kRrev);
                const double af = num(f, kAf);
                const double pmean = num(f, kPmean);
                const double pstd = num(f, kPstd);
                const double qual = num(f, kQual);
                const double mapq = num(f, kMapq);
                const double sn = num(f, kSn);
                const double hiaf = num(f, kHiaf);
                const double msi = num(f, kMsi);
                const double msilen = num(f, kMsilen);
                const double nm = num(f, kNm);
                const double hicnt = num(f, kHicnt);
                const double rd = rfwd + rrev;

                // teststrandbias.R: fisher.test(matrix(c(RefFwd, RefRev, AltFwd, AltRev), nrow = 2))
                const stats::FisherExact fisher(static_cast<int>(rfwd), static_cast<int>(rrev),
                                                static_cast<int>(num(f, kVfwd)), static_cast<int>(num(f, kVrev)),
                                                stats::FisherFlavor::r);
                const double sbf = fisher.p_value();
                const std::string sbf_text = util::rounded_value_to_print(sbf, 5);
                double oddratio = fisher.odds_ratio();
                std::string oddratio_text;
                if (std::isinf(oddratio)) {
                    oddratio = 0;
                    oddratio_text = "0";
                } else if (oddratio < 1 && oddratio > 0) {
                    oddratio = 1 / oddratio;
                    oddratio_text =
                        oddratio == std::floor(oddratio) ? integer_text(oddratio) : util::double_to_string(oddratio);
                } else {
                    oddratio_text = fisher.odd_ratio();
                }

                std::vector<std::string> filters;
                if (dp < o.min_depth.value && !(hicnt * hiaf >= 0.5)) {
                    filters.push_back(tag_d);
                }
                if (hicnt < o.min_var_depth.value && !(hicnt * hiaf >= 0.5)) {
                    filters.push_back(tag_v);
                }
                if (af < o.min_af.value) {
                    filters.push_back(tag_f);
                }
                if (pmean < o.min_pmean.value) {
                    filters.push_back(tag_p);
                }
                if (o.pstd && pstd == 0 && !isamp && af < 0.35) {
                    filters.push_back("pSTD");
                }
                if (qual < o.min_qual.value) {
                    filters.push_back(tag_q);
                }
                if (mapq < o.min_mq.value && af < 0.8) {
                    filters.push_back(tag_Q);
                }
                if (sn < o.min_sn.value) {
                    filters.push_back(tag_sn);
                }
                if (nm > o.max_nm.value) {
                    filters.push_back(tag_nm);
                }
                const double indel_len = std::abs(static_cast<double>(ref.size()) - static_cast<double>(alt.size()));
                if ((msi > o.msi.value && msilen > 1 && af < 0.2 && indel_len == msilen) ||
                    (msi >= 13 && msilen == 1 && af <= 0.275 && indel_len == msilen)) {
                    filters.push_back(tag_msi);
                }
                if (hiaf < 0.25 && str(f, kBias) == "2;1" && sbf < 0.01 && (oddratio > 5 || oddratio == 0) &&
                    end - start < 100) {
                    filters.push_back("Bias");
                }
                if (indel_len == msilen) {
                    if ((hiaf <= 0.275 && msi >= 13) || (hiaf <= 0.2 && msi >= 8 && msilen > 1)) {
                        filters.push_back("LongMSI");
                    }
                }
                if (!rejected && type == "SNV" && filters.empty() && start - pvs < o.cluster_bp.value) {
                    pfilter = tag_cluster;
                    filters.push_back(tag_cluster);
                }
                if (isamp && (num(f, kGamp) < num(f, kTamp) - num(f, kNcamp) || perl_true(str(f, kAmpflag)))) {
                    filters.push_back("AMPBIAS");
                }
                std::string filter = row->reject;
                for (const std::string& t : filters) {
                    filter += (filter.empty() ? "" : ";") + t;
                }
                if (filter.empty()) {
                    filter = "PASS";
                }

                std::string bias = str(f, kBias);
                if (auto p = bias.find(';'); p != std::string::npos) {
                    bias[p] = ':';
                }
                // ($vd le 1) ? 0 : int(log($vd)/log(2) * $qual); `le` compares strings.
                const std::string& vd_text = str(f, kVd);
                const double vd = num(f, kVd);
                const std::string qual_out =
                    vd_text <= std::string("1") ? "0" : integer_text(std::trunc(std::log(vd) / std::log(2.0) * qual));
                if (!rejected) {
                    flush();
                }

                std::string svinfo;
                if (!isamp) {
                    std::string splitreads = "0";
                    std::string spanpairs = "0";
                    const bool has_tamp = f.size() > kTamp;
                    if (has_tamp && str(f, kTamp).contains('-')) {
                        const std::string& t = str(f, kTamp);
                        const auto a = t.find('-');
                        const auto b = t.find('-', a + 1);
                        splitreads = t.substr(0, a);
                        spanpairs = t.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1);
                    }
                    if (alt.contains('<')) {
                        auto sr = util::parse_double(splitreads);
                        if (!((sr ? *sr : 0) >= o.min_split_reads.value)) {
                            if (!rejected) {
                                pinfo1.clear();
                            }
                            continue;
                        }
                        double svlen = end - start;
                        if (alt.contains("INV")) {
                            svlen += 1;
                        }
                        std::string svtype;
                        for (char c : alt) {
                            if (c != '<' && c != '>') {
                                svtype += c;
                            }
                        }
                        svinfo = ";SVTYPE=" + svtype + ";SVLEN=" + integer_text(svlen);
                    }
                    if (has_tamp) {
                        svinfo += ";SPLITREAD=" + splitreads + ";SPANPAIR=" + spanpairs;
                    }
                }

                std::string gt;
                if (ref == alt) {
                    alt = ".";
                    gt = "0/0";
                } else {
                    gt = (1 - af < o.hom_freq.value) ? "1/1" : af >= 0.5 ? "1/0" : af >= o.min_af.value ? "0/1" : "0/0";
                }
                const std::string rd_text = integer_text(rd);
                const std::string ad = vd_text != "0" ? rd_text + "," + vd_text : rd_text;
                std::string ampinfo;
                std::string dupinfo;
                std::string crispr;
                if (isamp) {
                    ampinfo = ";GDAMP=" + str(f, kGamp) + ";TLAMP=" + str(f, kTamp) + ";NCAMP=" + str(f, kNcamp) +
                              ";AMPFLAG=" + str(f, kAmpflag);
                } else {
                    if (f.size() > kGamp) {
                        dupinfo = ";DUPRATE=" + str(f, kGamp);
                    }
                    if (f.size() > kNcamp) {
                        crispr = ";CRISPR=" + str(f, kNcamp);
                    }
                }
                std::string info1 = chr + "\t" + str(f, kStart) + "\t.\t" + ref + "\t" + alt + "\t" + qual_out;
                std::string info2 =
                    "SAMPLE=" + sample_nowhitespace + ";TYPE=" + type + ";DP=" + str(f, kDp) + ";END=" + str(f, kEnd) +
                    ";VD=" + vd_text + ";AF=" + str(f, kAf) + ";BIAS=" + bias + ";REFBIAS=" + str(f, kRfwd) + ":" +
                    str(f, kRrev) + ";VARBIAS=" + str(f, kVfwd) + ":" + str(f, kVrev) + ";PMEAN=" + str(f, kPmean) +
                    ";PSTD=" + str(f, kPstd) + ";QUAL=" + str(f, kQual) + ";QSTD=" + str(f, kQstd) +
                    ";SBF=" + sbf_text + ";ODDRATIO=" + oddratio_text + ";MQ=" + str(f, kMapq) + ";SN=" + str(f, kSn) +
                    ";HIAF=" + str(f, kHiaf) + ";ADJAF=" + str(f, kAdjaf) + ";SHIFT3=" + str(f, kShift3) +
                    ";MSI=" + str(f, kMsi) + ";MSILEN=" + str(f, kMsilen) + ";NM=" + str(f, kNm) +
                    ";HICNT=" + str(f, kHicnt) + ";HICOV=" + str(f, kHicov) + ";LSEQ=" + str(f, kLseq) +
                    ";RSEQ=" + str(f, kRseq) + ampinfo + dupinfo + crispr + svinfo + "\tGT:DP:VD:AD:AF:RD:ALD\t" + gt +
                    ":" + str(f, kDp) + ":" + vd_text + ":" + ad + ":" + str(f, kAf) + ":" + str(f, kRfwd) + "," +
                    str(f, kRrev) + ":" + str(f, kVfwd) + "," + str(f, kVrev);
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

}  // namespace omnivar::annotate
