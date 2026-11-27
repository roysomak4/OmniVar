#include "cli/cli.hpp"

#include <htslib/hts.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <vector>

#include "util/num_format.hpp"
#include "util/number_parse.hpp"
#include "util/regex_engine.hpp"
#include "version.hpp"

namespace senri {

BamNames BamNames::from(const std::string& value) {
    BamNames b;
    b.raw = value;
    for (auto part : util::regex_split(value, util::Regex::compile("\\|").value())) {
        b.names.emplace_back(part);
    }
    if (!b.names.empty()) {
        for (auto part : util::regex_split(b.names[0], util::Regex::compile(":").value())) {
            b.bams.emplace_back(part);
        }
    }
    return b;
}

}  // namespace senri

namespace senri::cli {

namespace {

enum class ArgKind { none, required, optional };
enum class ValueType { string, number };

// One commons-cli Option as declared in CmdParser.buildOptions.
struct OptionSpec {
    std::string_view opt;       // short name (may be several letters, e.g. "th")
    std::string_view long_opt;  // may be empty
    ArgKind arg;
    ValueType type;
    std::string_view arg_name;
    std::string_view description;
};

// Declaration order of CmdParser.buildOptions; used for lookup and for the help page.
constexpr OptionSpec kOptions[] = {
    {"H", "?", ArgKind::none, ValueType::string, "", "Print this help page"},
    {"h", "header", ArgKind::none, ValueType::string, "", "Print a header row describing columns"},
    {"v", "", ArgKind::none, ValueType::string, "", "VCF format output"},
    {"i", "splice", ArgKind::none, ValueType::string, "", "Output splicing read counts"},
    {"p", "", ArgKind::none, ValueType::string, "",
     "Also report the variants the caller rejects, with the first failing rule as FILTER (REJ_*; tumour/normal: "
     "REJ_*_T, REJ_*_N; TSV: last column). The calls are unchanged."},
    {"C", "", ArgKind::none, ValueType::string, "",
     "Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)"},
    {"D", "debug", ArgKind::none, ValueType::string, "",
     "Debug mode.  Will print some error messages and append full genotype at the end."},
    {"t", "dedup", ArgKind::none, ValueType::string, "",
     "Indicate to remove duplicated reads.  Only one pair with same start positions will be kept"},
    {"3", "", ArgKind::none, ValueType::string, "",
     "Indicate to move indels to 3-prime if alternative alignment can be achieved."},
    {"K", "", ArgKind::none, ValueType::string, "", "Include Ns in the total depth calculation"},
    {"u", "", ArgKind::none, ValueType::string, "",
     "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using "
     "forward read only."},
    {"UN", "", ArgKind::none, ValueType::string, "",
     "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using "
     "first read only."},
    {"chimeric", "", ArgKind::none, ValueType::string, "", "Indicate to turn off chimeric reads filtering."},
    {"deldupvar", "", ArgKind::none, ValueType::string, "",
     "Turn on deleting of duplicate variants. Variants in this mode are considered and outputted only if start "
     "position of variant is inside the region interest."},
    {"fisher", "", ArgKind::none, ValueType::string, "",
     "Experimental feature: Changes R script (teststrandbias.R and testsomatic.) to the built-in implementation of "
     "Fisher "
     "exact test."},
    {"U", "nosv", ArgKind::none, ValueType::string, "", "Turn off structural variant calling."},
    {"F", "", ArgKind::required, ValueType::string, "bit",
     "The hexical to filter reads using samtools. Default: 0x504 (filter 2nd alignments, unmapped reads and "
     "duplicates).  Use -F 0 to turn it off."},
    {"z", "", ArgKind::optional, ValueType::number, "0/1",
     "Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use "
     "0 to turn it off. When using the -R option, it's set to 0"},
    {"k", "", ArgKind::optional, ValueType::number, "0/1",
     "Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it.  For Ion or PacBio, 0 "
     "is recommended."},
    {"a", "amplicon", ArgKind::required, ValueType::number, "int:float",
     "Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is "
     "considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap "
     "fraction is at least float.  Default: 10:0.95"},
    {"c", "", ArgKind::required, ValueType::number, "INT", "The column for chromosome"},
    {"G", "", ArgKind::required, ValueType::string, "Genome fasta",
     "The reference fasta. Should be indexed (.fai). A thin-FASTA CSV slice file (chr,start,end,sequence) is "
     "also accepted (regression tests)."},
    {"R", "", ArgKind::required, ValueType::string, "Region",
     "The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED "
     "is needed."},
    {"d", "", ArgKind::required, ValueType::string, "delemiter",
     "The delimiter for split region_info, default to tab \"\\t\""},
    {"n", "", ArgKind::required, ValueType::string, "regular_expression",
     "The regular expression to extract sample name from BAM filenames.  Default to: /([^\\/\\._]+?)_[^\\/]*.bam/"},
    {"N", "", ArgKind::required, ValueType::string, "string",
     "The sample name to be used directly.  Will overwrite -n option"},
    {"b", "", ArgKind::required, ValueType::string, "string", "The indexed BAM file"},
    {"S", "", ArgKind::required, ValueType::number, "INT", "The column for region start, e.g. gene start"},
    {"E", "", ArgKind::required, ValueType::number, "INT", "The column for region end, e.g. gene end"},
    {"s", "", ArgKind::required, ValueType::number, "INT",
     "The column for segment starts in the region, e.g. exon starts"},
    {"e", "", ArgKind::required, ValueType::number, "INT", "The column for segment ends in the region, e.g. exon ends"},
    {"g", "", ArgKind::required, ValueType::number, "INT", "The column for gene name, or segment annotation"},
    {"x", "", ArgKind::required, ValueType::number, "INT",
     "The number of nucleotide to extend for each segment, default: 0"},
    {"f", "", ArgKind::required, ValueType::number, "double",
     "The threshold for allele frequency, default: 0.01 or 1%"},
    {"r", "", ArgKind::required, ValueType::number, "minimum reads", "The minimum # of variant reads, default 2"},
    {"B", "", ArgKind::required, ValueType::number, "INT",
     "The minimum # of reads to determine strand bias, default 2"},
    {"Q", "", ArgKind::required, ValueType::number, "INT",
     "If set, reads with mapping quality less than INT will be filtered and ignored"},
    {"q", "", ArgKind::required, ValueType::number, "INT",
     "The phred score for a base to be considered a good call.  Default: 22.5 (for Illumina) For PGM, set it to "
     "~15, as PGM tends to under estimate base quality."},
    {"m", "", ArgKind::required, ValueType::number, "INT",
     "If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as "
     "mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - "
     "Indels.  Default: 8, or reads with more than 8 mismatches will not be used."},
    {"T", "trim", ArgKind::required, ValueType::number, "INT", "Trim bases after [INT] bases in the reads"},
    {"X", "", ArgKind::required, ValueType::number, "INT",
     "Extension of bp to look for mismatches after insersion or deletion.  Default to 2 bp, or only calls when "
     "they're within 2 bp."},
    {"P", "", ArgKind::required, ValueType::number, "number",
     "The read position filter.  If the mean variants position is less that specified, it's considered false "
     "positive.  Default: 5"},
    {"Z", "downsample", ArgKind::required, ValueType::number, "double",
     "For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with "
     "caution.  The downsampling will be random and non-reproducible."},
    {"o", "", ArgKind::required, ValueType::number, "Qratio",
     "The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.  Default: "
     "1.5"},
    {"O", "", ArgKind::required, ValueType::number, "MapQ",
     "The reads should have at least mean MapQ to be considered a valid variant.  Default: no filtering"},
    {"V", "", ArgKind::required, ValueType::number, "freq",
     "The lowest frequency in the normal sample allowed for a putative somatic mutation.  Defaults to 0.05"},
    {"I", "", ArgKind::required, ValueType::number, "INT", "The indel size.  Default: 50bp"},
    {"y", "verbose", ArgKind::none, ValueType::string, "", ""},
    {"th", "", ArgKind::optional, ValueType::number, "INT", "Threads count."},
    {"M", "", ArgKind::required, ValueType::number, "INT",
     "The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, "
     "then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the "
     "matching is only the primers. Default: 0, or no filtering"},
    {"VS", "", ArgKind::required, ValueType::string, "STRICT | LENIENT | SILENT",
     "How strict to be when reading a SAM or BAM. STRICT - throw an exception if something looks wrong. LENIENT - "
     "Emit warnings but keep going if possible. SILENT - Like LENIENT, only don't emit warning messages. Default: "
     "LENIENT"},
    {"DP", "default-printer", ArgKind::required, ValueType::string, "string",
     "The printer type used for different outputs. Default: OUT (i.e. System.out)."},
    {"A", "", ArgKind::required, ValueType::number, "INT",
     "The number of STD. A pair will be considered for DEL if INSERT > INSERT_SIZE + INSERT_STD_AMT * INSERT_STD.  "
     "Default: 4"},
    {"W", "insert-std", ArgKind::required, ValueType::number, "INT",
     "The insert size STD.  Used for SV calling.  Default: 100"},
    {"w", "insert-size", ArgKind::required, ValueType::number, "INT",
     "The insert size.  Used for SV calling.  Default: 300"},
    {"L", "", ArgKind::required, ValueType::number, "INT",
     "The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc. Default: 1000. "
     "Any indel, complex variants less than this will be spelled out with exact nucleotides."},
    {"Y", "ref-extension", ArgKind::required, ValueType::number, "INT",
     "Extension of bp of reference to build lookup table. Default to 1200 bp. Increase the number will slowdown "
     "the program. The main purpose is to call large indels with 1000 bp that can be missed by discordant mate "
     "pairs."},
    {"adaptor", "", ArgKind::required, ValueType::string, "string",
     "Filter adaptor sequences so that they are not used in realignment. Multiple adaptors can be supplied by "
     "setting them with comma, like: --adaptor ACGTTGCTC,ACGGGGTCTC,ACGCGGCTAG ."},
    {"J", "crispr", ArgKind::required, ValueType::number, "CRISPR_cutting_site",
     "The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the "
     "guide.  For CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end sites to as close "
     "to this location as possible, if there are alternatives. The option should only be used for CRISPR mode."},
    {"j", "", ArgKind::required, ValueType::number, "CRISPR_filtering_bp",
     "In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not "
     "meet the criteria, it will not be used for variant calling, since it is likely just a partially amplified "
     "PCR.  Default: not set, or no filtering"},
    {"mfreq", "", ArgKind::required, ValueType::number, "double",
     "The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25"},
    {"nmfreq", "", ArgKind::required, ValueType::number, "double",
     "The variant frequency threshold to determine variant as good in case of non-monomer MSI. Default: 0.1"},
};

// Util.stripLeadingHyphens (commons-cli 1.2)
std::string_view strip_leading_hyphens(std::string_view s) {
    if (s.starts_with("--")) {
        return s.substr(2);
    }
    if (s.starts_with("-")) {
        return s.substr(1);
    }
    return s;
}

// Util.stripLeadingAndTrailingQuotes (commons-cli 1.2)
std::string strip_quotes(std::string_view s) {
    if (s.starts_with("\"")) {
        s.remove_prefix(1);
    }
    if (s.ends_with("\"")) {
        s.remove_suffix(1);
    }
    return std::string(s);
}

// Options.getOption / hasOption: short name first, then long name.
const OptionSpec* find_option(std::string_view token) {
    const std::string_view name = strip_leading_hyphens(token);
    for (const auto& o : kOptions) {
        if (o.opt == name) {
            return &o;
        }
    }
    for (const auto& o : kOptions) {
        if (!o.long_opt.empty() && o.long_opt == name) {
            return &o;
        }
    }
    return nullptr;
}

// The parsed commons-cli CommandLine: options in the order seen, with their values.
struct CommandLine {
    struct Entry {
        const OptionSpec* spec;
        std::vector<std::string> values;
    };
    std::vector<Entry> options;
    std::vector<std::string> args;

    // CommandLine.hasOption: by short or long name.
    bool has(std::string_view name) const { return first(name) != nullptr; }

    const Entry* first(std::string_view name) const {
        for (const auto& e : options) {
            if (e.spec->opt == name || (!e.spec->long_opt.empty() && e.spec->long_opt == name)) {
                return &e;
            }
        }
        return nullptr;
    }

    // CommandLine.getOptionValue: first value over all occurrences, in order.
    std::optional<std::string> value(std::string_view name) const {
        for (const auto& e : options) {
            if ((e.spec->opt == name || (!e.spec->long_opt.empty() && e.spec->long_opt == name)) && !e.values.empty()) {
                return e.values.front();
            }
        }
        return std::nullopt;
    }
};

std::expected<CommandLine, ParseError> commons_cli_parse(const std::vector<std::string_view>& tokens,
                                                         bool& missing_required_b) {
    CommandLine cmd;
    bool seen_b = false;
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string_view t = tokens[i];
        bool eat_the_rest = false;
        if (t == "--") {
            eat_the_rest = true;
        } else if (t != "-" && t.starts_with("-")) {  // a lone "-" is an argument
            const OptionSpec* spec = find_option(t);
            if (spec == nullptr) {
                return std::unexpected(ParseError{"Unrecognized option: " + std::string(t)});
            }
            CommandLine::Entry entry{spec, {}};
            if (spec->opt == "b") {
                seen_b = true;
            }
            if (spec->arg != ArgKind::none) {
                // processArgs: take following tokens while they are not options starting with '-';
                // each Option here accepts one value, so the second attempt ends the loop.
                while (i + 1 < tokens.size()) {
                    const std::string_view next = tokens[i + 1];
                    if (next.starts_with("-") && find_option(next) != nullptr) {
                        break;
                    }
                    if (!entry.values.empty()) {
                        break;  // "Cannot add value, list full."
                    }
                    entry.values.push_back(strip_quotes(next));
                    ++i;
                }
                if (entry.values.empty() && spec->arg == ArgKind::required) {
                    return std::unexpected(ParseError{"Missing argument for option: " + std::string(spec->opt)});
                }
            }
            cmd.options.push_back(std::move(entry));
        } else {
            cmd.args.emplace_back(t);
        }
        if (eat_the_rest) {
            for (++i; i < tokens.size(); ++i) {
                if (tokens[i] != "--") {
                    cmd.args.emplace_back(tokens[i]);
                }
            }
            break;
        }
    }
    missing_required_b = !seen_b;
    return cmd;
}

// CommandLine.getParsedOptionValue for a Number-typed option: null if absent or without value.
std::expected<std::optional<util::ParsedNumber>, ParseError> parsed_number(const CommandLine& cmd,
                                                                           std::string_view name) {
    auto v = cmd.value(name);
    if (!v) {
        return std::optional<util::ParsedNumber>{};
    }
    auto n = util::create_number(*v);
    if (!n) {
        return std::unexpected(ParseError{n.error()});
    }
    return std::optional<util::ParsedNumber>{*n};
}

// VCF annotation thresholds (var2vcf_valid.pl options, by long name). Returns the field an
// argument "--filter-<name>" or "--filter-<name>=VALUE" sets, or nullopt.
std::optional<annotate::Threshold*> vcf_threshold(std::string_view a, annotate::Var2VcfOptions& o) {
    const std::string_view name = a.substr(0, a.find('='));
    const std::pair<std::string_view, annotate::Threshold*> table[] = {
        {"--filter-min-depth", &o.min_depth},
        {"--filter-min-var-depth", &o.min_var_depth},
        {"--filter-min-af", &o.min_af},
        {"--filter-min-pmean", &o.min_pmean},
        {"--filter-min-qual", &o.min_qual},
        {"--filter-min-mq", &o.min_mq},
        {"--filter-min-sn", &o.min_sn},
        {"--filter-max-nm", &o.max_nm},
        {"--filter-msi", &o.msi},
        {"--filter-cluster-bp", &o.cluster_bp},
        {"--filter-hom-af", &o.hom_freq},
        {"--filter-min-split-reads", &o.min_split_reads},
        {"--filter-max-pvalue", &o.max_pvalue},
    };
    for (const auto& [n, field] : table) {
        if (name == n) {
            return field;
        }
    }
    return std::nullopt;
}

}  // namespace

std::expected<ParsedCommandLine, ParseError> parse_command_line(std::span<const std::string_view> all_args,
                                                                int available_processors) {
    ParsedCommandLine out;

    // Added options (long form only, never seen by the upstream rules).
    std::vector<std::string_view> tokens;
    bool want_version = false;
    bool want_help = false;
    bool after_double_dash = false;
    bool output_tsv = false;
    annotate::Var2VcfOptions vcf;
    for (std::size_t ai = 0; ai < all_args.size(); ++ai) {
        std::string_view a = all_args[ai];
        if (!after_double_dash) {
            if (auto t = vcf_threshold(a, vcf)) {
                // --filter-<name> VALUE or --filter-<name>=VALUE
                std::string_view value;
                if (auto eq = a.find('='); eq != std::string_view::npos) {
                    value = a.substr(eq + 1);
                } else if (ai + 1 < all_args.size()) {
                    value = all_args[++ai];
                } else {
                    return std::unexpected(ParseError{"Missing argument for option: " + std::string(a)});
                }
                auto v = util::parse_double(value);
                if (!v || std::isnan(*v)) {
                    return std::unexpected(
                        ParseError{"Invalid number for option " + std::string(a) + ": " + std::string(value)});
                }
                **t = annotate::Threshold{*v, std::string(value), true};
                continue;
            }
            if (a == "--") {
                after_double_dash = true;
            } else if (a == "--tsv") {
                output_tsv = true;
                continue;
            } else if (a == "--vcf-top-allele-only") {
                vcf.top_allele_only = true;
                continue;
            } else if (a == "--vcf-somatic-only") {
                vcf.somatic_only = true;
                continue;
            } else if (a == "--no-filter-pstd") {
                vcf.pstd = false;
                continue;
            } else if (a == "--version") {
                want_version = true;
                continue;
            } else if (a == "--help") {
                want_help = true;
                continue;
            } else if (a == "--dev-tsv-to-vcf") {
                if (ai + 1 >= all_args.size()) {
                    return std::unexpected(ParseError{"Missing argument for option: --dev-tsv-to-vcf"});
                }
                out.tsv_path = std::string(all_args[++ai]);
                continue;
            }
        }
        tokens.push_back(a);
    }
    if (want_version) {
        out.action = Action::version;
        return out;
    }
    if (want_help) {
        out.action = Action::help;
        return out;
    }
    if (!out.tsv_path.empty()) {
        out.action = Action::tsv_to_vcf;
        out.config.vcf = vcf;
        return out;
    }

    bool missing_b = false;
    auto parsed = commons_cli_parse(tokens, missing_b);
    if (!parsed) {
        return std::unexpected(parsed.error());
    }
    const CommandLine& cmd = *parsed;

    // CmdParser.parseParams: checkRequiredOptions runs inside parse, so a
    // missing -b prints its message and the help page before the "no options / -H" check.
    if (missing_b) {
        out.warnings += "Missing required option(s): b\n";
        out.action = Action::help;
        return out;
    }
    if (cmd.options.empty() || cmd.has("H")) {
        out.action = Action::help;
        return out;
    }

    Config& c = out.config;
    auto num = [&](std::string_view name) { return parsed_number(cmd, name); };
    auto int_value = [&](std::string_view name, int def) -> std::expected<int, ParseError> {
        auto v = num(name);
        if (!v) {
            return std::unexpected(v.error());
        }
        return *v ? (*v)->int_value() : def;
    };
    auto column_value = [&](std::string_view name, int def) -> std::expected<int, ParseError> {
        auto v = num(name);
        if (!v) {
            return std::unexpected(v.error());
        }
        return *v ? (*v)->int_value() - 1 : def;
    };
    auto double_value = [&](std::string_view name, double def) -> std::expected<double, ParseError> {
        auto v = num(name);
        if (!v) {
            return std::unexpected(v.error());
        }
        return *v ? (*v)->double_value() : def;
    };
#define VD_TRY(dst, expr)                         \
    do {                                          \
        auto vd_r = (expr);                       \
        if (!vd_r) {                              \
            return std::unexpected(vd_r.error()); \
        }                                         \
        (dst) = *vd_r;                            \
    } while (0)

    if (!cmd.args.empty()) {
        c.bed = cmd.args[0];
    }
    c.print_header = cmd.has("h");
    c.chromosome_name_is_number = cmd.has("C");
    c.debug = cmd.has("D");
    c.remove_duplicated_reads = cmd.has("t");
    c.move_indels_to_3 = cmd.has("3");
    c.samfilter = cmd.value("F").value_or("0x504");

    if (cmd.has("z")) {
        int z = 0;
        VD_TRY(z, int_value("z", 1));
        c.zero_based = z == 1;
    }
    if (auto a = cmd.value("a")) {
        c.amplicon_based_calling = *a;
    }
    {
        int k = 0;
        VD_TRY(k, int_value("k", 1));
        c.perform_local_realignment = k == 1;
    }

    if (auto g = cmd.value("G")) {
        c.fasta = *g;
        if (*g == "hg19") {
            c.fasta = Config::kHg19;
        } else if (*g == "hg38") {
            c.fasta = Config::kHg38;
        } else if (*g == "mm10") {
            c.fasta = Config::kMm10;
        }
    } else {
        out.warnings += "Reference file path wasn't set (option -G). Will be used the default fasta path.\n";
        c.fasta = Config::kHg19;
    }

    if (auto r = cmd.value("R")) {
        c.region_of_interest = *r;
    }
    c.delimiter = cmd.value("d").value_or("\t");
    if (auto n = cmd.value("N")) {
        c.sample_name = *n;
    }
    if (cmd.has("n")) {
        std::string regexp = cmd.value("n").value_or("");
        if (regexp.starts_with("/")) {
            regexp.erase(0, 1);
        }
        if (regexp.ends_with("/")) {
            regexp.pop_back();
        }
        c.sample_name_regexp = regexp;
    }
    c.bam = BamNames::from(cmd.value("b").value_or(""));
    if (c.bam.names.empty() || c.bam.bams.empty()) {
        // upstream: ArrayIndexOutOfBoundsException on bamNames[0]
        return std::unexpected(ParseError{"Invalid -b value: " + c.bam.raw});
    }

    int c_col = 0, s_col_upper = 0, e_col_upper = 0, s_col = 0, e_col = 0, g_col = 0;
    VD_TRY(c_col, column_value("c", kDefaultBedRowFormat.chr_column));
    VD_TRY(s_col_upper, column_value("S", kDefaultBedRowFormat.start_column));
    VD_TRY(e_col_upper, column_value("E", kDefaultBedRowFormat.end_column));
    VD_TRY(s_col, column_value("s", kDefaultBedRowFormat.thick_start_column));
    VD_TRY(e_col, column_value("e", kDefaultBedRowFormat.thick_end_column));
    VD_TRY(g_col, column_value("g", kDefaultBedRowFormat.gene_column));
    if (cmd.has("S") && !cmd.has("s")) {
        s_col = s_col_upper;
    }
    if (cmd.has("E") && !cmd.has("e")) {
        e_col = e_col_upper;
    }
    c.bed_row_format = BedRowFormat{c_col, s_col_upper, e_col_upper, s_col, e_col, g_col};
    VD_TRY(c.column_for_chromosome, column_value("c", -1));

    VD_TRY(c.number_nucleotide_to_extend, int_value("x", 0));
    VD_TRY(c.freq, double_value("f", 0.01));
    VD_TRY(c.minr, int_value("r", 2));
    VD_TRY(c.min_bias_reads, int_value("B", 2));
    if (cmd.has("Q")) {
        std::optional<util::ParsedNumber> q;
        VD_TRY(q, num("Q"));
        if (q) {
            c.mapping_quality = q->int_value();
        }
    }
    VD_TRY(c.goodq, double_value("q", 22.5));
    VD_TRY(c.mismatch, int_value("m", 8));
    VD_TRY(c.trim_bases_after, int_value("T", 0));
    VD_TRY(c.vext, int_value("X", 2));
    VD_TRY(c.read_pos_filter, int_value("P", 5));
    if (cmd.has("Z")) {
        double z = 0;
        VD_TRY(z, double_value("Z", 0));
        c.downsampling = z;
    }
    VD_TRY(c.qratio, double_value("o", 1.5));
    VD_TRY(c.mapq, double_value("O", 0));
    VD_TRY(c.lofreq, double_value("V", 0.05));
    VD_TRY(c.indelsize, int_value("I", 50));

    if (cmd.has("p")) {
        c.do_pileup = true;
    }
    c.y = cmd.has("y");
    c.output_splicing = cmd.has("i");

    if (cmd.has("M")) {
        std::optional<util::ParsedNumber> m;
        VD_TRY(m, num("M"));
        if (m) {
            c.minmatch = m->int_value();
        }
    }
    if (cmd.has("VS")) {
        std::string vs = cmd.value("VS").value_or("");
        std::transform(vs.begin(), vs.end(), vs.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
        if (vs == "STRICT") {
            c.validation_stringency = ValidationStringency::strict;
        } else if (vs == "LENIENT") {
            c.validation_stringency = ValidationStringency::lenient;
        } else if (vs == "SILENT") {
            c.validation_stringency = ValidationStringency::silent;
        } else {
            return std::unexpected(ParseError{"No enum constant htsjdk.samtools.ValidationStringency." + vs});
        }
    }

    c.include_n_in_total_depth = cmd.has("K");
    c.chimeric = cmd.has("chimeric");
    c.disable_sv = cmd.has("U");
    c.unique_mode_second_in_pair_enabled = cmd.has("UN");
    c.unique_mode_alignment_enabled = cmd.has("u");
    c.delete_duplicate_variants = cmd.has("deldupvar");

    VD_TRY(c.ins_size, int_value("w", 300));
    VD_TRY(c.ins_std, int_value("W", 100));
    VD_TRY(c.ins_std_amt, int_value("A", 4));
    VD_TRY(c.sv_min_len, int_value("L", 1000));

    {
        int threads = 0;
        if (cmd.has("th")) {
            std::optional<util::ParsedNumber> th;
            VD_TRY(th, num("th"));
            threads = th ? th->int_value() : available_processors;
        }
        c.threads = std::max(threads, 1);
    }
    VD_TRY(c.reference_extension, int_value("Y", 1200));

    if (cmd.has("adaptor")) {
        const std::string value = cmd.value("adaptor").value_or("");
        for (auto part : util::regex_split(value, util::Regex::compile(",").value())) {
            c.adaptor.emplace_back(part);
        }
    }
    c.fisher = cmd.has("fisher");
    if (cmd.has("DP")) {
        c.printer_type = cmd.value("DP").value_or("OUT") == "ERR" ? PrinterType::err : PrinterType::out;
    }
    VD_TRY(c.crispr_cutting_site, int_value("J", 0));
    VD_TRY(c.crispr_filtering_bp, int_value("j", 0));
    VD_TRY(c.monomer_msi_frequency, double_value("mfreq", 0.25));
    VD_TRY(c.non_monomer_msi_frequency, double_value("nmfreq", 0.1));
#undef VD_TRY

    out.config.output_tsv = output_tsv;
    out.config.vcf = vcf;
    return out;
}

std::string usage() {
    std::string s =
        "usage: senri [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] "
        "[-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info\n"
        "senri is a variant calling program for SNV, MNV, indels (<50 bp), complex and structural variants.\n"
        "Options:\n";
    for (const auto& o : kOptions) {
        std::string head = " -" + std::string(o.opt);
        if (!o.long_opt.empty()) {
            head += ",--" + std::string(o.long_opt);
        }
        if (o.arg == ArgKind::required || o.arg == ArgKind::optional) {
            head += " <" + std::string(o.arg_name) + ">";
        }
        s += head;
        if (!o.description.empty()) {
            s += head.size() < 28 ? std::string(28 - head.size(), ' ') : std::string("\n") + std::string(28, ' ');
            s += o.description;
        }
        s += '\n';
    }
    s += "\nOutput:\n"
         " --tsv                      Write the tab-separated variant table instead of VCF\n"
         " --vcf-top-allele-only      VCF: only the highest-AF allele per position (var2vcf_valid.pl without -A)\n"
         "\nVCF FILTER tags (var2vcf_valid.pl checks; they annotate, never remove records):\n"
         " --filter-min-depth <int>        d<x>: total depth below (waived if HICNT*HIAF >= 0.5). Default: 3\n"
         " --filter-min-var-depth <int>    v<x>: high-quality variant reads below (same waiver). Default: 2\n"
         " --filter-min-af <float>         f<x>: AF below; also the 0/1 genotype limit. Default: 0.02\n"
         " --filter-min-pmean <float>      p<x>: mean position in reads below. Default: 8\n"
         " --filter-min-qual <float>       q<x>: mean base quality below. Default: 22.5\n"
         " --filter-min-mq <float>         Q<x>: mean mapping quality below, when AF < 0.8. Default: 10\n"
         " --filter-min-sn <float>         SN<x>: signal to noise below. Default: 1.5\n"
         " --filter-max-nm <float>         NM<x>: mean mismatches above. Default: 5.25\n"
         " --filter-msi <float>            MSI<x>: non-monomer MSI above (monomer: 13), for low AF. Default: 12\n"
         " --filter-cluster-bp <int>       Cluster<x>bp: PASS SNVs closer than this. Default: 0 (off)\n"
         " --filter-hom-af <float>         GT 1/1 when 1-AF is below this. Default: 0.2\n"
         " --filter-min-split-reads <int>  SVs with fewer split reads are not written. Default: 1\n"
         " --no-filter-pstd                Do not tag pSTD (position STD of 0), e.g. for PCR panels\n"
         "\nTumour/normal VCF (two BAMs, -b \"T|N\"; var2vcf_paired.pl checks): the options above apply to\n"
         "each sample, with the paired script's defaults d5 v3 f0.02 p8 q22.5 Q0 SN1.5 NM5.25 MSI12. Also:\n"
         " --filter-max-pvalue <float>     P<x> with --vcf-somatic-only: tumour-vs-normal p above. Default: 0.05\n"
         " --vcf-somatic-only              var2vcf_paired.pl -M: P<x>, P0.01Likely, InDelLikely; no germline PASS\n"
         "\n --version                  Print version information\n";
    return s;
}

std::string version_string() { return "senri " + std::string(kVersion) + ", htslib " + hts_version(); }

}  // namespace senri::cli
