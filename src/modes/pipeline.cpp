#include "modes/pipeline.hpp"

#include "core/cigar_parser.hpp"
#include "core/realigner.hpp"
#include "core/vars_builder.hpp"
#include "util/regex_engine.hpp"

namespace senri::modes {

std::vector<std::string> split_bams(const std::string& bam) {
    static const util::Regex colon = util::Regex::compile(":").value();
    std::vector<std::string> out;
    for (auto p : util::regex_split(bam, colon)) {
        out.emplace_back(p);
    }
    return out;
}

PipelineResult run_pipeline(const PipelineEnv& env, const Region& region, const std::string& bam, int max_read_length,
                            io::Reference& reference, std::unordered_set<std::string>& splice) {
    core::RegionScope scope{env.conf,       env.run,       region,     reference, splice,
                            env.exceptions, &env.provider, &env.cache, bam};
    auto pre = io::RecordPreprocessor::create(split_bams(bam), region, env.conf, env.cache);
    if (!pre) {
        throw core::FatalError(pre.error().message);
    }
    core::CigarParser parser(scope, false, max_read_length);
    core::CigarParserResult parsed = parser.process(*pre);
    core::realigner_process(parsed.data, scope);
    core::sv_process(parsed.data, scope);
    PipelineResult result;
    result.max_read_length = parsed.data.max_read_length;
    result.aligned = core::build_vars(parsed.data, scope);
    return result;
}

}  // namespace senri::modes
