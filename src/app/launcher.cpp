#include "app/launcher.hpp"

#include <algorithm>
#include <charconv>
#include <condition_variable>
#include <cstdio>
#include <map>
#include <mutex>
#include <optional>
#include <print>
#include <set>
#include <stop_token>
#include <thread>
#include <unordered_set>
#include <vector>

#include "annotate/var2vcf.hpp"
#include "cli/cli.hpp"
#include "core/cigar_parser.hpp"
#include "core/realigner.hpp"
#include "core/vars_builder.hpp"
#include "io/bam.hpp"
#include "io/record_preprocessor.hpp"
#include "io/reference.hpp"
#include "modes/pipeline.hpp"
#include "modes/simple.hpp"
#include "modes/somatic.hpp"
#include "region/region_builder.hpp"
#include "util/regex_engine.hpp"
#include "util/sequence.hpp"

namespace omnivar::app {

namespace {

constexpr std::string_view kSamplePattern = "([^\\/\\._]+).sorted[^\\/]*.bam";
constexpr std::string_view kSamplePattern2 = "([^\\/]+)[_\\.][^\\/]*bam";

std::optional<std::string> find_group1(const util::Regex& re, std::string_view text) {
    auto m = re.find(text);
    if (!m) {
        return std::nullopt;
    }
    auto g = m->group(text, 1);
    // Matcher.group(1) of an unmatched group is null
    return g ? std::optional<std::string>(std::string(*g)) : std::nullopt;
}

std::vector<std::string> split_colon(const std::string& s) {
    static const util::Regex colon = util::Regex::compile(":").value();
    std::vector<std::string> out;
    for (auto p : util::regex_split(s, colon)) {
        out.emplace_back(p);
    }
    return out;
}

}  // namespace

std::expected<std::pair<std::string, std::string>, std::string> sample_names(const Config& conf) {
    std::optional<std::string> sample;
    std::string samplem;
    std::optional<util::Regex> user_re;
    if (conf.sample_name_regexp && !conf.sample_name_regexp->empty()) {
        auto re = util::Regex::compile(*conf.sample_name_regexp);
        if (!re) {
            return std::unexpected("unsupported -n regular expression: " + re.error());
        }
        user_re.emplace(std::move(*re));
    }
    static const util::Regex sample_re = util::Regex::compile(kSamplePattern).value();
    static const util::Regex sample_re2 = util::Regex::compile(kSamplePattern2).value();

    if (conf.sample_name) {
        sample = *conf.sample_name;
    } else {
        sample = find_group1(user_re ? *user_re : sample_re, conf.bam.raw);
    }
    if (!sample) {
        sample = find_group1(sample_re2, conf.bam.raw);
    }

    if (conf.region_of_interest && conf.bam.has_bam2()) {
        // getSampleNamesSomatic
        if (conf.sample_name_regexp) {
            auto re = util::Regex::compile(*conf.sample_name_regexp);
            if (!re) {
                return std::unexpected("unsupported -n regular expression: " + re.error());
            }
            if (auto m = re->find(conf.bam.bam1())) {
                auto g = m->group(conf.bam.bam1(), 1);
                sample = g ? std::optional<std::string>(std::string(*g)) : std::nullopt;
            }
            if (auto m = re->find(conf.bam.bam2())) {
                auto g = m->group(conf.bam.bam2(), 1);
                samplem = g ? std::string(*g) : std::string("null");
            }
        } else if (conf.sample_name) {
            static const util::Regex bar = util::Regex::compile("\\|").value();
            auto split = util::regex_split(*conf.sample_name, bar);
            if (!split.empty()) {
                sample = std::string(split[0]);
                samplem = split.size() > 1 ? std::string(split[1]) : std::string(split[0]) + "_match";
            }
        }
    }
    return std::pair{sample.value_or("null"), samplem};
}

std::expected<RunContext, std::string> init_resources(const Config& conf) {
    if (!conf.region_of_interest && !conf.bed) {
        return std::unexpected(
            std::string("The required BED file or region missed, please, set it with path to BED or with -R option."));
    }
    RunContext ctx;
    {
        auto bam = io::BamFile::open(conf.bam.bam_x());
        if (!bam) {
            return std::unexpected(bam.error().message);
        }
        ctx.chr_lengths = bam->chr_lengths();
        ctx.contigs = bam->contigs();
    }
    auto samples = sample_names(conf);
    if (!samples) {
        return std::unexpected(samples.error());
    }
    ctx.sample = samples->first;
    ctx.samplem = samples->second;

    auto builder = RegionBuilder::create(ctx.chr_lengths, conf);
    if (!builder) {
        return std::unexpected(builder.error());
    }
    if (conf.region_of_interest) {
        auto segs = builder->build_region_from_configuration();
        if (!segs) {
            return std::unexpected(segs.error());
        }
        ctx.segments = std::move(*segs);
    } else {
        auto bed = read_bed_file(conf);
        if (!bed) {
            return std::unexpected(bed.error());
        }
        ctx.amplicon_based_calling = bed->amplicon_parameters;
        auto segs = ctx.amplicon_based_calling ? builder->build_amp_regions(bed->lines, bed->zero_based.value_or(false))
                                               : builder->build_regions(bed->lines, bed->zero_based);
        if (!segs) {
            return std::unexpected(segs.error());
        }
        ctx.segments = std::move(*segs);
    }

    // Fill adaptor maps
    for (const std::string& sequence : conf.adaptor) {
        for (int i = 0; i <= 6 && i + Config::kAdSeed < static_cast<int>(sequence.size()); ++i) {
            const std::string forward = sequence.substr(static_cast<std::size_t>(i), Config::kAdSeed);
            const std::string reverse = util::complement(util::reverse(forward));
            ctx.adaptor_forward[forward] = i + 1;
            ctx.adaptor_reverse[reverse] = i + 1;
        }
    }
    return ctx;
}

int tsv_to_vcf(const Config& conf, const std::string& path) {
    std::FILE* in = path == "-" ? stdin : std::fopen(path.c_str(), "rb");
    if (in == nullptr) {
        std::println(stderr, "omnivar: cannot open {}", path);
        return 1;
    }
    std::string text;
    char buf[1 << 16];
    for (std::size_t n; (n = std::fread(buf, 1, sizeof buf, in)) > 0;) {
        text.append(buf, n);
    }
    if (in != stdin) {
        (void)std::fclose(in);
    }
    // Somatic rows (SomaticOutputVariant) have 55 columns, Simple ones 36 (+ CRISPR).
    const std::string_view first_line = std::string_view(text).substr(0, text.find('\n'));
    const bool paired = std::count(first_line.begin(), first_line.end(), '\t') + 1 >= 48;
    annotate::Var2Vcf vcf(conf.vcf, paired);
    std::string sample;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        if (nl == std::string::npos) {
            nl = text.size();
        }
        const std::string_view line = std::string_view(text).substr(start, nl - start);
        if (!line.empty()) {
            vcf.add_row(line);
            sample = std::string(line.substr(0, line.find('\t')));  // var2vcf: the last row's sample
        }
        start = nl + 1;
    }
    annotate::VcfHeaderInfo info;
    info.source = cli::version_string();
    info.sample = sample;
    std::print("{}", vcf.header(info));
    std::print("{}", vcf.records(sample));
    return 0;
}

namespace {

// One region through the Simple-mode pipeline (SimpleMode.processBamInPipeline / AbstractMode.pipeline).
void process_simple_region(const modes::PipelineEnv& env, const Region& region, std::string& out) {
    const Config& conf = env.conf;
    // AbstractMode.tryToGetReference: a failure stops the run
    auto reference = io::get_reference(env.provider, region, conf, env.run.chr_lengths);
    if (!reference) {
        throw core::FatalError(reference.error().message);
    }
    std::unordered_set<std::string> splice;
    if (conf.output_splicing) {
        core::RegionScope scope{conf,           env.run,       region,     *reference,     splice,
                                env.exceptions, &env.provider, &env.cache, conf.bam.bam1()};
        auto pre = io::RecordPreprocessor::create(split_colon(conf.bam.bam1()), region, conf, env.cache);
        if (!pre) {
            throw core::FatalError(pre.error().message);
        }
        core::CigarParser parser(scope, false, 0);
        out += parser.process(*pre).splicing_output;
        return;
    }
    modes::PipelineResult result = modes::run_pipeline(env, region, conf.bam.bam1(), 0, *reference, splice);
    core::RegionScope scope{conf,           env.run,       region,     *reference,     splice,
                            env.exceptions, &env.provider, &env.cache, conf.bam.bam1()};
    modes::simple_post_process(result.aligned, scope, out);
}

}  // namespace

int run(const Config& conf) {
    auto ctx = init_resources(conf);
    if (!ctx) {
        std::println(stderr, "omnivar: {}", ctx.error());
        return 1;
    }
    const bool splicing = conf.output_splicing;
    if (!splicing && !(conf.region_of_interest || !ctx->amplicon_based_calling)) {
        // Amplicon mode (an 8-column BED, or -a) is out of scope for this release.
        std::println(stderr,
                     "omnivar: amplicon mode (8-column amplicon BED or -a) is not supported in this "
                     "version; use a 3- or 4-column BED");
        return 2;
    }
    const bool somatic = !splicing && conf.bam.has_bam2();
    // VCF (default): the rows are formatted as the TSV prints them without --fisher or -D, which is
    // the input teststrandbias.R | var2vcf_valid.pl expects, and converted at the end (the
    // script orders chromosomes over the whole run).
    const bool vcf_output = !splicing && !conf.output_tsv;
    Config row_conf = conf;
    if (vcf_output) {
        row_conf.fisher = false;
        row_conf.debug = false;
    }
    // -p: also report the variants the caller rejects, with the first failing rule.
    const bool rejected_output = conf.do_pileup && !splicing;
    annotate::Var2VcfOptions vcf_options = conf.vcf;
    if (rejected_output && !somatic) {
        vcf_options.reject_filters = core::reject_filter_descriptions(conf);
    } else if (rejected_output) {
        // tumour/normal: each reason per sample, plus the somatic module's own drops
        for (const char* sample : {"T", "N"}) {
            const std::string name = sample[0] == 'T' ? "tumour" : "normal";
            for (const auto& [id, description] : core::reject_filter_descriptions(conf)) {
                const std::string pre = "Rejected by the caller (-p): ";
                vcf_options.reject_filters.emplace_back(id + "_" + sample,
                                                        pre + name + " " + description.substr(pre.size()));
            }
        }
        vcf_options.reject_filters.emplace_back(
            "REJ_Masked_T",
            "Rejected by the caller (-p): tumour allele dropped because a higher-ranked tumour allele "
            "at the position was rejected");
        vcf_options.reject_filters.emplace_back(
            "REJ_Combine_T",
            "Rejected by the caller (-p): tumour-only indel not supported when both BAMs are "
            "combined (combineAnalysis)");
        vcf_options.reject_filters.emplace_back(
            "REJ_Combine_N",
            "Rejected by the caller (-p): normal-only indel not supported when both BAMs are "
            "combined (combineAnalysis)");
    }
    annotate::Var2Vcf vcf(vcf_options, somatic);
    const auto emit = [&](const std::string& text) {
        if (!vcf_output && !rejected_output) {
            std::print("{}", text);
            return;
        }
        std::string tsv;
        std::size_t start = 0;
        while (start < text.size()) {
            auto nl = text.find('\n', start);
            if (nl == std::string::npos) {
                nl = text.size();
            }
            if (nl > start) {
                std::string_view line = std::string_view(text).substr(start, nl - start);
                std::string_view reason;
                if (line.front() == modes::kRejectedLineMarker) {
                    const auto tab = line.find('\t');
                    reason = line.substr(1, tab - 1);
                    line.remove_prefix(tab + 1);
                }
                if (vcf_output) {
                    vcf.add_row(line, reason);
                } else {
                    tsv.append(line).append("\t").append(reason.empty() ? std::string_view("PASS") : reason) += '\n';
                }
            }
            start = nl + 1;
        }
        std::print("{}", tsv);
    };
    const auto write_vcf = [&] {
        if (!vcf_output) {
            return;
        }
        std::string vcf_sample = ctx->sample;
        std::string vcf_samplem;
        if (somatic) {
            vcf_samplem = ctx->samplem;
            if (vcf_samplem.empty() && conf.sample_name && conf.sample_name->contains('|')) {
                const auto bar = conf.sample_name->find('|');
                vcf_sample = conf.sample_name->substr(0, bar);
                const auto bar2 = conf.sample_name->find('|', bar + 1);  // split(/\|/): field 2 only
                vcf_samplem =
                    conf.sample_name->substr(bar + 1, bar2 == std::string::npos ? std::string::npos : bar2 - bar - 1);
            }
        }
        annotate::VcfHeaderInfo info;
        info.source = cli::version_string();
        info.sample = vcf_sample;
        info.samplem = vcf_samplem;
        info.reference = conf.fasta;
        info.contigs = ctx->contigs;
        std::print("{}", vcf.header(info));
        std::print("{}", vcf.records(vcf_sample));
    };

    if (splicing) {
        if (conf.print_header) {
            std::print("Sample\tChr\tIntron\tIntron count\n");
        }
    } else if (!vcf_output) {
        std::string header = somatic ? modes::somatic_header(conf) : modes::simple_header(conf);
        if (rejected_output && !header.empty()) {
            header.insert(header.size() - 1, "\tFilter");  // -p: PASS or the REJ_* reason
        }
        std::print("{}", header);
    }
    (void)std::fflush(stdout);

    std::vector<const Region*> regions;
    for (const auto& list : ctx->segments) {
        for (const Region& region : list) {
            regions.push_back(&region);
        }
    }
    core::ExceptionCounter exceptions;

    // One region on one thread's resources (reference reader and open BAMs are per thread:
    // faidx and htslib handles are not thread-safe). A FatalError message, if any, is returned
    // with the output produced before it.
    struct RegionResult {
        std::string out;
        std::optional<std::string> fatal;
    };
    const auto process = [&](const Region& region, const io::ReferenceProvider* provider,
                             const std::string& provider_error, io::BamCache& cache) {
        RegionResult r;
        try {
            if (provider == nullptr) {
                throw core::FatalError(provider_error);
            }
            const modes::PipelineEnv env{row_conf, *ctx, *provider, cache, exceptions};
            if (somatic) {
                modes::process_somatic_region(env, region, r.out);
            } else {
                process_simple_region(env, region, r.out);
            }
        } catch (const core::FatalError& e) {
            r.fatal = e.what();
        }
        return r;
    };
    // Emits one result in region order; false when it was fatal (the run stops).
    const auto finish = [&](const Region& region, const RegionResult& r) {
        emit(r.out);
        if (!r.fatal) {
            return true;
        }
        write_vcf();
        (void)std::fflush(stdout);
        std::println(stderr, "Critical exception occurs on region: {}:{}-{}, program will be stopped.", region.chr,
                     region.start, region.end);
        std::println(stderr, "{}", *r.fatal);
        return false;
    };

    const int threads = std::max(1, conf.threads);
    if (threads == 1 || regions.size() <= 1) {
        auto provider = io::open_reference(conf.fasta);
        io::BamCache cache;
        for (const Region* region : regions) {
            const RegionResult r = process(*region, provider ? provider->get() : nullptr,
                                           provider ? std::string() : provider.error().message, cache);
            if (!finish(*region, r)) {
                return 1;
            }
        }
    } else {
        // AbstractParallelMode: a fixed pool; output in region order. Workers stay at most
        // `window` regions ahead of the output (upstream: a queue of CAPACITY = 10 futures).
        const std::size_t window = std::max<std::size_t>(10, 4 * static_cast<std::size_t>(threads));
        std::vector<std::optional<RegionResult>> results(regions.size());
        std::mutex mutex;
        std::condition_variable ready;         // a result was stored
        std::condition_variable_any advanced;  // the output moved on (or a stop was requested)
        std::size_t next = 0;
        std::size_t emitted = 0;
        // Declared after the state the workers use: destroying `workers` requests stop and joins
        // every thread before that state goes away. std::jthread needs libc++ 20 (Clang ≥ 20).
        std::vector<std::jthread> workers;
        workers.reserve(static_cast<std::size_t>(threads));
        for (int t = 0; t < threads; ++t) {
            workers.emplace_back([&](const std::stop_token& stop) {
                auto provider = io::open_reference(conf.fasta);
                const std::string provider_error = provider ? std::string() : provider.error().message;
                io::BamCache cache;
                while (true) {
                    std::size_t i = 0;
                    {
                        std::unique_lock lock(mutex);
                        advanced.wait(lock, stop, [&] { return next >= regions.size() || next < emitted + window; });
                        if (stop.stop_requested() || next >= regions.size()) {
                            return;
                        }
                        i = next++;
                    }
                    RegionResult r = process(*regions[i], provider ? provider->get() : nullptr, provider_error, cache);
                    {
                        std::lock_guard lock(mutex);
                        results[i] = std::move(r);
                    }
                    ready.notify_all();
                }
            });
        }
        int status = 0;
        for (std::size_t i = 0; i < regions.size(); ++i) {
            RegionResult r;
            {
                std::unique_lock lock(mutex);
                ready.wait(lock, [&] { return results[i].has_value(); });
                r = std::move(*results[i]);
                results[i].reset();
            }
            const bool ok = finish(*regions[i], r);
            if (!ok) {
                status = 1;
                break;
            }
            {
                std::lock_guard lock(mutex);
                emitted = i + 1;
            }
            advanced.notify_all();
        }
        // On a fatal region, stop the workers waiting for work; each finishes its current region.
        for (std::jthread& w : workers) {
            w.request_stop();
        }
        workers.clear();
        if (status != 0) {
            return status;
        }
    }
    write_vcf();
    (void)std::fflush(stdout);
    return 0;
}

}  // namespace omnivar::app
