#pragma once

// REF/ALT alleles and genotypes.

#include "core/scope.hpp"
#include "core/variant.hpp"
#include "core/variation.hpp"

namespace senri::core {

struct MsiResult {
    double msi;
    int shift3;
    std::string msint;
};

// ToVarsBuilder.findMSI; `left` empty means null.
[[nodiscard]] MsiResult find_msi(const std::string& tseq1, const std::string& tseq2, const std::string* left);

[[nodiscard]] AlignedVars build_vars(VariationData& data, const RegionScope& scope);

}  // namespace senri::core
