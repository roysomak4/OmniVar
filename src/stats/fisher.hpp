#pragma once

// R's fisher.test for 2x2 tables: p-values and the conditional-MLE odds ratio. Uses ports of
// commons-math3 HypergeometricDistribution / SaddlePointExpansion (Apache License 2.0).
// Verified against the JVM: tests/unit/data/fisher_vectors.tsv.
//
// Two flavours:
//   jvm   the upstream FisherExact, used for the `--fisher` TSV columns.
//   r     R 4.3's fisher.test as run by teststrandbias.R / testsomatic.R, used by the VCF
//         annotation layer: no 7-dp rounding of the log-densities (FisherExact.logdcDhyper),
//         long-double sums (R's sum()), and dhyper with R's numerical choices (see fisher.cpp).
//         Equal to R on every table in tests/unit/data/fisher_r_vectors.tsv

#include <string>
#include <vector>

namespace senri::stats {

enum class FisherFlavor { jvm, r };

class FisherExact {
public:
    FisherExact(int ref_fwd, int ref_rev, int alt_fwd, int alt_rev, FisherFlavor flavor = FisherFlavor::jvm);

    [[nodiscard]] double p_value() const;          // getPValue: two-sided, rounded like R (5 dp)
    [[nodiscard]] std::string odd_ratio() const;   // getOddRatio: "Inf", integer, or Double.toString
    [[nodiscard]] double odds_ratio() const;       // conditional MLE rounded to 5 dp (+Inf allowed)
    [[nodiscard]] double p_value_less() const;     // getPValueLess
    [[nodiscard]] double p_value_greater() const;  // getPValueGreater

private:
    [[nodiscard]] double mle(double x) const;
    [[nodiscard]] double mnhyper(double ncp) const;
    [[nodiscard]] std::vector<double> dnhyper(double ncp) const;
    [[nodiscard]] double pnhyper(int q, bool upper_tail) const;

    [[nodiscard]] double sum(const std::vector<double>& values) const;

    FisherFlavor flavor_;
    int m_;
    int n_;
    int k_;
    int x_;
    int lo_;
    int hi_;
    std::vector<int> support_;
    std::vector<double> logdc_;
    double p_less_ = 0;
    double p_greater_ = 0;
    double p_two_sided_ = 0;
};

// DoubleStream.sum() (JDK 21): Kahan summation (Collectors.sumWithCompensation / computeFinalSum).
[[nodiscard]] double stream_sum(const std::vector<double>& values);

// commons-math3 HypergeometricDistribution(populationSize, numberOfSuccesses, sampleSize)
[[nodiscard]] double hypergeometric_log_probability(int population, int successes, int sample, int x);
[[nodiscard]] double hypergeometric_cumulative(int population, int successes, int sample, int x);
[[nodiscard]] double hypergeometric_upper_cumulative(int population, int successes, int sample, int x);

}  // namespace senri::stats
