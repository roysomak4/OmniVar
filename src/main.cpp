#include <cstdio>
#include <print>
#include <string_view>
#include <thread>
#include <vector>

#include "app/launcher.hpp"
#include "cli/cli.hpp"

int main(int argc, char** argv) {
    std::vector<std::string_view> args(argv + 1, argv + argc);

    const int cpus = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    auto parsed = senri::cli::parse_command_line(args, cpus);
    if (!parsed) {
        std::println(stderr, "senri: {}", parsed.error().message);
        return 1;
    }
    if (!parsed->warnings.empty()) {
        std::print(stderr, "{}", parsed->warnings);
    }

    switch (parsed->action) {
        case senri::cli::Action::help:
            // The help page, exit status 0.
            std::print("{}", senri::cli::usage());
            return 0;
        case senri::cli::Action::version: std::println("{}", senri::cli::version_string()); return 0;
        case senri::cli::Action::tsv_to_vcf: return senri::app::tsv_to_vcf(parsed->config, parsed->tsv_path);
        case senri::cli::Action::run: return senri::app::run(parsed->config);
    }
    return 2;
}
