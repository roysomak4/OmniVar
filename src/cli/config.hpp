#pragma once

// Immutable run configuration. Field names follow the upstream Configuration (snake_case);
// defaults are upstream's. Nullable upstream fields (Boolean, Integer, Double) are std::optional. Filled by
// cli::parse_command_line (port of CmdParser.parseCmd).

#include <optional>
#include <string>
#include <vector>

#include "annotate/var2vcf.hpp"

namespace omnivar {

struct BedRowFormat {
    int chr_column;
    int start_column;
    int end_column;
    int thick_start_column;
    int thick_end_column;
    int gene_column;

    bool operator==(const BedRowFormat&) const = default;
};

// RegionBuilder.DEFAULT_BED_ROW_FORMAT (refGene.txt layout)
inline constexpr BedRowFormat kDefaultBedRowFormat{2, 6, 7, 9, 10, 12};
inline constexpr BedRowFormat kCustomBedRowFormat{0, 1, 2, 3, 1, 2};
inline constexpr BedRowFormat kAmpBedRowFormat{0, 1, 2, 6, 7, 3};

struct BamNames {
    std::string raw;                 // the -b value
    std::vector<std::string> names;  // raw split on '|' (String.split: trailing empties dropped)
    std::vector<std::string> bams;   // names[0] split on ':'

    [[nodiscard]] const std::string& bam1() const { return names.at(0); }
    [[nodiscard]] bool has_bam2() const { return names.size() > 1; }
    [[nodiscard]] const std::string& bam2() const { return names.at(1); }
    [[nodiscard]] const std::string& bam_x() const { return bams.at(0); }

    [[nodiscard]] static BamNames from(const std::string& value);
};

enum class ValidationStringency { strict, lenient, silent };
enum class PrinterType { out, err };

struct Config {
    // Fixed constants (Configuration static finals)
    static constexpr const char* kHg19 = "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
    static constexpr const char* kHg38 = "/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa";
    static constexpr const char* kMm10 = "/ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa";
    static constexpr int kLowQual = 10;
    static constexpr int kSeed1 = 17;
    static constexpr int kSeed2 = 12;
    static constexpr int kAdSeed = 6;
    static constexpr double kMinSvcDist = 1.5;
    static constexpr int kMinMapBase = 15;
    static constexpr int kMinSvPos = 25;
    static constexpr int kSvMaxLen = 150000;
    static constexpr int kSvFlank = 50;
    static constexpr int kDiscPairQual = 35;
    static constexpr int kExtension = 5000;
    static constexpr const char* kDefaultAmpliconParameters = "10:0.95";
    static constexpr int kMaxExceptionCount = 10;

    bool print_header = false;                          // -h
    std::string delimiter = "\t";                       // -d (a regex, util/regex_engine)
    std::optional<std::string> bed;                     // first positional argument
    int number_nucleotide_to_extend = 0;                // -x
    std::optional<bool> zero_based;                     // -z
    std::optional<std::string> amplicon_based_calling;  // -a
    int column_for_chromosome = -1;                     // -c (0-based), -1 if not set
    BedRowFormat bed_row_format = kDefaultBedRowFormat;
    std::optional<std::string> sample_name_regexp;  // -n
    std::optional<std::string> sample_name;         // -N
    std::string fasta;                              // -G
    BamNames bam;                                   // -b
    std::optional<double> downsampling;             // -Z
    bool chromosome_name_is_number = false;         // -C
    std::optional<int> mapping_quality;             // -Q
    bool remove_duplicated_reads = false;           // -t
    int mismatch = 8;                               // -m
    bool y = false;                                 // -y (verbose)
    double goodq = 22.5;                            // -q
    int vext = 2;                                   // -X
    int trim_bases_after = 0;                       // -T
    bool perform_local_realignment = true;          // -k
    int indelsize = 50;                             // -I
    double bias = 0.05;
    int min_bias_reads = 2;                                                      // -B
    int minr = 2;                                                                // -r
    bool debug = false;                                                          // -D
    double freq = 0.01;                                                          // -f
    bool move_indels_to_3 = false;                                               // -3
    std::string samfilter = "0x504";                                             // -F
    std::optional<std::string> region_of_interest;                               // -R
    int read_pos_filter = 5;                                                     // -P
    double qratio = 1.5;                                                         // -o
    double mapq = 0;                                                             // -O
    bool do_pileup = false;                                                      // -p
    double lofreq = 0.05;                                                        // -V
    int minmatch = 0;                                                            // -M
    bool output_splicing = false;                                                // -i
    ValidationStringency validation_stringency = ValidationStringency::lenient;  // -VS
    bool include_n_in_total_depth = false;                                       // -K
    bool unique_mode_alignment_enabled = false;                                  // -u
    bool unique_mode_second_in_pair_enabled = false;                             // -UN
    int threads = 1;                                                             // -th
    bool chimeric = false;                                                       // -chimeric
    bool disable_sv = false;                                                     // -U
    bool delete_duplicate_variants = false;                                      // -deldupvar
    bool fisher = false;                                                         // -fisher
    int ins_size = 300;                                                          // -w  (INSSIZE)
    int ins_std = 100;                                                           // -W  (INSSTD)
    int ins_std_amt = 4;                                                         // -A  (INSSTDAMT)
    int sv_min_len = 1000;                                                       // -L  (SVMINLEN)
    int reference_extension = 1200;                                              // -Y
    PrinterType printer_type = PrinterType::out;                                 // -DP
    std::vector<std::string> adaptor;                                            // -adaptor
    int crispr_filtering_bp = 0;                                                 // -j
    int crispr_cutting_site = 0;                                                 // -J
    double monomer_msi_frequency = 0.25;                                         // -mfreq
    double non_monomer_msi_frequency = 0.1;                                      // -nmfreq

    // Added options (long options only; not part of the upstream Configuration).
    bool output_tsv = false;       // --tsv: the TSV table instead of VCF
    annotate::Var2VcfOptions vcf;  // --filter-*, --vcf-top-allele-only

    [[nodiscard]] bool is_column_for_chromosome_set() const { return column_for_chromosome >= 0; }
    [[nodiscard]] bool is_downsampling() const { return downsampling.has_value(); }
    [[nodiscard]] bool has_mapping_quality() const { return mapping_quality.has_value(); }
    [[nodiscard]] bool is_zero_based_defined() const { return zero_based.has_value(); }
};

}  // namespace omnivar
