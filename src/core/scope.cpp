#include "core/scope.hpp"

#include <cstdio>
#include <print>

namespace omnivar::core {

void print_exception_and_continue(const std::exception& e, const std::string& place, const std::string& place_def,
                                  const Region* region, ExceptionCounter& counter) {
    const std::string first = "There was Exception while processing " + place + " on " + place_def;
    const std::string second = ". The processing will be continued from the next " + place + ".";
    if (region != nullptr) {
        std::println(stderr, "{} on region {}{}", first, region->print_region(), second);
    } else {
        std::println(stderr, "{} but region is undefined{}", first, second);
    }
    std::println(stderr, "{}", e.what());
    if (const int current = ++counter.count; current > Config::kMaxExceptionCount) {
        throw FatalError("omnivar fails (there were " + std::to_string(counter.count.load()) +
                         " continued exceptions during the run).");
    }
}

}  // namespace omnivar::core
