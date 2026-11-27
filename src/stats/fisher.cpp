#include "stats/fisher.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "util/num_format.hpp"

namespace senri::stats {

namespace {

// ---- commons-math3 3.6.1 SaddlePointExpansion (Apache License 2.0; see THIRD_PARTY_LICENSES) ----
constexpr double kHalfLog2Pi = 0x1.d67f1c864beb4p-1;
constexpr double kExactStirlingErrors[] = {
    0x0.0p0,
    0x1.3a37a020b8c22p-3,
    0x1.4c071bcda0a5bp-4,
    0x1.c1098b28dcf33p-5,
    0x1.52a9b923ea649p-5,
    0x1.0fab9626b44ffp-5,
    0x1.c579a268d80b3p-6,
    0x1.850ea113caf0ep-6,
    0x1.54a2662fd78a9p-6,
    0x1.2eea2e990f134p-6,
    0x1.10b4e513fcbedp-6,
    0x1.eff15b81c9cc5p-7,
    0x1.c6b167bebdf36p-7,
    0x1.a3c5f8a1e7d1dp-7,
    0x1.85d4d612e4a86p-7,
    0x1.6bdfcc7fbdb0ap-7,
    0x1.552805e7b3076p-7,
    0x1.411b75e41049cp-7,
    0x1.2f4871b12ab64p-7,
    0x1.1f553026fbce1p-7,
    0x1.10f9d4c0743a7p-7,
    0x1.03fc2d49c8fb9p-7,
    0x1.f0593088014f8p-8,
    0x1.dac773cc5b3cdp-8,
    0x1.c7018733aa9c6p-8,
    0x1.b4d04a067629fp-8,
    0x1.a40514700f36cp-8,
    0x1.947826547fe6cp-8,
    0x1.86076c002d4a7p-8,
    0x1.789583ffc5e8cp-8,
    0x1.6c08f6f194a1p-8,
};
constexpr double kS0 = 0.083333333333333333333;         // 1/12
constexpr double kS1 = 0.00277777777777777777778;       // 1/360
constexpr double kS2 = 0.00079365079365079365079365;    // 1/1260
constexpr double kS3 = 0.000595238095238095238095238;   // 1/1680
constexpr double kS4 = 0.0008417508417508417508417508;  // 1/1188
constexpr double kTwoPi = 6.283185307179586;            // MathUtils.TWO_PI

// SaddlePointExpansion.getStirlingError, for the integer arguments used here.
double stirling_error(double z) {
    if (z < 15.0) {
        const double z2 = 2.0 * z;
        if (std::floor(z2) == z2) {
            return kExactStirlingErrors[static_cast<int>(z2)];
        }
        // Gamma.logGamma path: never reached (z is an integer count)
        return std::lgamma(z + 1.0) - (z + 0.5) * std::log(z) + z - kHalfLog2Pi;
    }
    const double z2 = z * z;
    return (kS0 - (kS1 - (kS2 - (kS3 - kS4 / z2) / z2) / z2) / z2) / z;
}

// SaddlePointExpansion.getDeviancePart
double deviance_part(double x, double mu) {
    if (std::fabs(x - mu) < 0.1 * (x + mu)) {
        const double d = x - mu;
        double v = d / (x + mu);
        double s1 = v * d;
        double s = std::numeric_limits<double>::quiet_NaN();
        double ej = 2.0 * x * v;
        v = v * v;
        int j = 1;
        while (s1 != s) {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1);
            ++j;
        }
        return s1;
    }
    return x * std::log(x / mu) + mu - x;
}

// SaddlePointExpansion.logBinomialProbability
double log_binomial_probability(int x, int n, double p, double q) {
    if (x == 0) {
        if (p < 0.1) {
            return -deviance_part(n, n * q) - n * p;
        }
        return n * std::log(q);
    }
    if (x == n) {
        if (q < 0.1) {
            return -deviance_part(n, n * p) - n * q;
        }
        return n * std::log(p);
    }
    double ret = stirling_error(n) - stirling_error(x) - stirling_error(n - x) - deviance_part(x, n * p) -
                 deviance_part(n - x, n * q);
    const double f = (kTwoPi * x * (n - x)) / n;
    return -0.5 * std::log(f) + ret;
}

// ---- R flavour (R 4.3 dhyper via dbinom_raw): the same saddle-point method (Loader 2000) with
// R's numerical choices, which make results differ from commons-math in the last bits:
//   - the exact Stirling-error table is used for n <= 15 (commons: n < 15);
//   - the Stirling series is shortened for large n (> 500, > 80, > 35 terms);
//   - the normalising term is log(2 pi) + log(x) + log1p(-x/n) (commons: log(2 pi x (n-x) / n)).
// The deviance part and the x == 0 / x == n branches are the same as commons-math.
constexpr double kLn2Pi = 1.837877066409345483560659472811;  // M_LN_2PI

double r_stirling_error(double n) {
    if (n <= 15.0) {
        const double nn = n + n;
        // integer counts only, so nn is always integral here
        return kExactStirlingErrors[static_cast<int>(nn)];
    }
    const double nn = n * n;
    if (n > 500) {
        return (kS0 - kS1 / nn) / n;
    }
    if (n > 80) {
        return (kS0 - (kS1 - kS2 / nn) / nn) / n;
    }
    if (n > 35) {
        return (kS0 - (kS1 - (kS2 - kS3 / nn) / nn) / nn) / n;
    }
    return (kS0 - (kS1 - (kS2 - (kS3 - kS4 / nn) / nn) / nn) / nn) / n;
}

double r_log_binomial_raw(int x, int n, double p, double q) {
    if (p == 0) {
        return x == 0 ? 0 : -std::numeric_limits<double>::infinity();
    }
    if (q == 0) {
        return x == n ? 0 : -std::numeric_limits<double>::infinity();
    }
    if (x == 0) {
        if (n == 0) {
            return 0;
        }
        return p < 0.1 ? -deviance_part(n, n * q) - n * p : n * std::log(q);
    }
    if (x == n) {
        return q < 0.1 ? -deviance_part(n, n * p) - n * q : n * std::log(p);
    }
    if (x < 0 || x > n) {
        return -std::numeric_limits<double>::infinity();
    }
    const double lc = r_stirling_error(n) - r_stirling_error(x) - r_stirling_error(n - x) - deviance_part(x, n * p) -
                      deviance_part(n - x, n * q);
    const double lf = kLn2Pi + std::log(static_cast<double>(x)) + std::log1p(-static_cast<double>(x) / n);
    return lc - 0.5 * lf;
}

// dhyper(x, r, b, n, log = TRUE): r white, b black balls, n drawn.
double r_dhyper_log(int x, int r, int b, int n) {
    if (x < 0 || x > r || n - x < 0 || n - x > b) {
        return -std::numeric_limits<double>::infinity();
    }
    if (n == 0) {
        return x == 0 ? 0 : -std::numeric_limits<double>::infinity();
    }
    const double p = static_cast<double>(n) / static_cast<double>(r + b);
    const double q = static_cast<double>(r + b - n) / static_cast<double>(r + b);
    const double p1 = r_log_binomial_raw(x, r, p, q);
    const double p2 = r_log_binomial_raw(n - x, b, p, q);
    const double p3 = r_log_binomial_raw(n, r + b, p, q);
    return p1 + p2 - p3;
}

// dhyper(x, r, b, n, log = FALSE): p1 * p2 / p3 of the non-log dbinom_raw values (R 4.3 dhyper.c).
double r_dhyper(int x, int r, int b, int n) {
    if (x < 0 || x > r || n - x < 0 || n - x > b) {
        return 0;
    }
    if (n == 0) {
        return x == 0 ? 1 : 0;
    }
    const double p = static_cast<double>(n) / static_cast<double>(r + b);
    const double q = static_cast<double>(r + b - n) / static_cast<double>(r + b);
    const double p1 = std::exp(r_log_binomial_raw(x, r, p, q));
    const double p2 = std::exp(r_log_binomial_raw(n - x, b, p, q));
    const double p3 = std::exp(r_log_binomial_raw(n, r + b, p, q));
    return p1 * p2 / p3;
}

// phyper(x, NR, NB, n, lower.tail): R 4.3 phyper.c, the d(x) * pdhyper(x) form with the
// long-double series of pdhyper, after swapping tails when x is above the mean.
double r_phyper(double x, double nr, double nb, double n, bool lower_tail) {
    x = std::floor(x + 1e-7);
    if (x * (nr + nb) > n * nr) {  // Swap tails.
        std::swap(nr, nb);
        x = n - x - 1;
        lower_tail = !lower_tail;
    }
    if (x < 0 || x < n - nb) {
        return lower_tail ? 0 : 1;
    }
    if (x >= nr || x >= n) {
        return lower_tail ? 1 : 0;
    }
    const double d = r_dhyper(static_cast<int>(x), static_cast<int>(nr), static_cast<int>(nb), static_cast<int>(n));
    if (d == 0.) {
        return lower_tail ? 0 : 1;
    }
    // pdhyper
    long double sum = 0;
    long double term = 1;
    double xx = x;
    while (xx > 0 && term >= std::numeric_limits<double>::epsilon() * sum) {
        term *= xx * (nb - n + xx) / (n + 1 - xx) / (nr + 1 - xx);
        sum += term;
        xx--;
    }
    const double pd = 1 + static_cast<double>(sum);
    return lower_tail ? d * pd : (0.5 - d * pd + 0.5);
}

double hypergeometric_probability(int population, int successes, int sample, int x) {
    const double lp = hypergeometric_log_probability(population, successes, sample, x);
    return lp == -std::numeric_limits<double>::infinity() ? 0 : std::exp(lp);
}

// UnirootZeroIn.zeroinC
double zeroin(double ax, double bx, const std::function<double(double)>& f, double tol) {
    const double epsilon = std::ldexp(1.0, -52);  // Math.ulp(1.0)
    double a = ax;
    double b = bx;
    double fa = f(a);
    double fb = f(b);
    double c = a;
    double fc = fa;
    for (;;) {
        const double prev_step = b - a;
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const double tol_act = 2 * epsilon * std::fabs(b) + tol / 2.0;
        double new_step = (c - b) / 2.0;
        if (std::fabs(new_step) <= tol_act || fb == 0.0) {
            return b;
        }
        if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
            double p;
            double q;
            const double cb = c - b;
            if (a == c) {
                const double t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                q = fa / fc;
                const double t1 = fb / fc;
                const double t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if (p > 0.0) {
                q = -q;
            } else {
                p = -p;
            }
            if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2.0) && p < std::fabs(prev_step * q / 2.0)) {
                new_step = p / q;
            }
        }
        if (std::fabs(new_step) < tol_act) {
            new_step = new_step > 0.0 ? tol_act : -tol_act;
        }
        a = b;
        fa = fb;
        b += new_step;
        fb = f(b);
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a;
            fc = fa;
        }
    }
}

// FisherExact.round_as_r
double round_as_r(double value) {
    constexpr double kResultRoundR = 1E5;
    value = util::round_half_even(value * kResultRoundR, 0);
    value = value / kResultRoundR;
    return value == 0.0 ? 0 : (value == 1.0 ? 1 : value);
}

}  // namespace

double stream_sum(const std::vector<double>& values) {
    double s0 = 0;
    double s1 = 0;
    double simple = 0;
    for (double value : values) {
        const double tmp = value - s1;
        const double sum = s0;
        const double velvel = sum + tmp;
        s1 = (velvel - sum) - tmp;
        s0 = velvel;
        simple += value;
    }
    const double tmp = s0 - s1;
    if (std::isnan(tmp) && std::isinf(simple)) {
        return simple;
    }
    return tmp;
}

double hypergeometric_log_probability(int population, int successes, int sample, int x) {
    const int lower = std::max(0, successes - (population - sample));
    const int upper = std::min(sample, successes);
    if (x < lower || x > upper) {
        return -std::numeric_limits<double>::infinity();
    }
    const double p = static_cast<double>(sample) / static_cast<double>(population);
    const double q = static_cast<double>(population - sample) / static_cast<double>(population);
    const double p1 = log_binomial_probability(x, successes, p, q);
    const double p2 = log_binomial_probability(sample - x, population - successes, p, q);
    const double p3 = log_binomial_probability(sample, population, p, q);
    return p1 + p2 - p3;
}

double hypergeometric_cumulative(int population, int successes, int sample, int x) {
    const int lower = std::max(0, successes - (population - sample));
    const int upper = std::min(sample, successes);
    if (x < lower) {
        return 0.0;
    }
    if (x >= upper) {
        return 1.0;
    }
    // innerCumulativeProbability(lower, x, 1)
    int x0 = lower;
    double ret = hypergeometric_probability(population, successes, sample, x0);
    while (x0 != x) {
        x0 += 1;
        ret += hypergeometric_probability(population, successes, sample, x0);
    }
    return ret;
}

double hypergeometric_upper_cumulative(int population, int successes, int sample, int x) {
    const int lower = std::max(0, successes - (population - sample));
    const int upper = std::min(sample, successes);
    if (x <= lower) {
        return 1.0;
    }
    if (x > upper) {
        return 0.0;
    }
    int x0 = upper;
    double ret = hypergeometric_probability(population, successes, sample, x0);
    while (x0 != x) {
        x0 -= 1;
        ret += hypergeometric_probability(population, successes, sample, x0);
    }
    return ret;
}

// FisherFlavor::r skips the 7-dp rounding of the log-densities (see fisher.hpp).
FisherExact::FisherExact(int ref_fwd, int ref_rev, int alt_fwd, int alt_rev, FisherFlavor flavor)
    : flavor_(flavor), m_(ref_fwd + ref_rev), n_(alt_fwd + alt_rev), k_(ref_fwd + alt_fwd), x_(ref_fwd) {
    lo_ = std::max(0, k_ - n_);
    hi_ = std::min(k_, m_);
    for (int j = lo_; j <= hi_; ++j) {
        support_.push_back(j);
    }
    // logdcDhyper
    for (int element : support_) {
        if (m_ + n_ == 0) {
            logdc_.push_back(0.0);
            continue;
        }
        // R: logdc <- dhyper(support, m, n, k, log = TRUE)
        double value = flavor == FisherFlavor::jvm ? hypergeometric_log_probability(m_ + n_, m_, k_, element)
                                                   : r_dhyper_log(element, m_, n_, k_);
        if (std::isnan(value)) {
            value = 0.0;
        }
        logdc_.push_back(flavor == FisherFlavor::jvm ? util::round_half_even(value, 7) : value);
    }
    // calculatePValue
    p_less_ = pnhyper(x_, false);
    p_greater_ = pnhyper(x_, true);
    const double rel_err = 1 + 1E-7;
    const std::vector<double> d = dnhyper(1.0);
    // R: sum(d[d <= d[1 + x - lo] * relErr]), accumulated in long double like R's sum().
    std::vector<double> kept;
    for (double el : d) {
        if (el <= d.at(static_cast<std::size_t>(x_ - lo_)) * rel_err) {
            kept.push_back(el);
        }
    }
    if (flavor_ == FisherFlavor::jvm) {
        double total = 0.0;
        for (double el : kept) {
            total += el;
        }
        p_two_sided_ = total;
    } else {
        p_two_sided_ = sum(kept);
    }
}

double FisherExact::mle(double x) const {
    const double eps = std::ldexp(1.0, -52);
    if (x == lo_) {
        return 0.0;
    }
    if (x == hi_) {
        return std::numeric_limits<double>::infinity();
    }
    const double mu = mnhyper(1.0);
    const double tol = std::pow(eps, 0.25);
    if (mu > x) {
        return zeroin(0, 1, [&](double t) { return mnhyper(t) - x; }, tol);
    }
    if (mu < x) {
        return 1.0 / zeroin(eps, 1, [&](double t) { return mnhyper(1.0 / t) - x; }, tol);
    }
    return 1.0;
}

double FisherExact::mnhyper(double ncp) const {
    if (ncp == 0) {
        return lo_;
    }
    if (std::isinf(ncp)) {
        return hi_;
    }
    const std::vector<double> d = dnhyper(ncp);
    std::vector<double> multiply;
    multiply.reserve(support_.size());
    for (std::size_t i = 0; i < support_.size(); ++i) {
        multiply.push_back(support_[i] * d[i]);
    }
    return sum(multiply);
}

std::vector<double> FisherExact::dnhyper(double ncp) const {
    std::vector<double> result;
    result.reserve(support_.size());
    for (std::size_t i = 0; i < support_.size(); ++i) {
        result.push_back(logdc_[i] + std::log(ncp) * support_[i]);
    }
    // Collections.max (Double.compareTo order)
    double max_result = result.at(0);
    for (double r : result) {
        if (util::double_compare(r, max_result) > 0) {
            max_result = r;
        }
    }
    std::vector<double> exponent;
    exponent.reserve(result.size());
    for (double el : result) {
        exponent.push_back(std::exp(el - max_result));
    }
    const double total = sum(exponent);
    for (std::size_t i = 0; i < exponent.size(); ++i) {
        result[i] = exponent[i] / total;
    }
    return result;
}

double FisherExact::pnhyper(int q, bool upper_tail) const {
    if (m_ + n_ == 0) {
        return 1.0;
    }
    if (flavor_ == FisherFlavor::r) {
        // R: pnhyper(q, ncp = 1) is phyper(x - 1, m, n, k, lower.tail = FALSE) or phyper(x, m, n, k)
        return upper_tail ? r_phyper(q - 1, m_, n_, k_, false) : r_phyper(q, m_, n_, k_, true);
    }
    return upper_tail ? hypergeometric_upper_cumulative(m_ + n_, m_, k_, q)
                      : hypergeometric_cumulative(m_ + n_, m_, k_, q);
}

std::string FisherExact::odd_ratio() const {
    const double r = mle(x_);
    if (std::isinf(r)) {
        return "Inf";
    }
    if (r == static_cast<double>(util::math_round(r))) {
        return util::decimal_format(r, 0);
    }
    return util::double_to_string(round_as_r(r));
}

// DoubleStream.sum (Kahan). R: sum() accumulates in long double (rsum, LDOUBLE).
double FisherExact::sum(const std::vector<double>& values) const {
    if (flavor_ == FisherFlavor::jvm) {
        return stream_sum(values);
    }
    long double s = 0;
    for (double v : values) {
        s += v;
    }
    return static_cast<double>(s);
}

double FisherExact::odds_ratio() const {
    const double r = mle(x_);
    return std::isinf(r) ? r : round_as_r(r);
}

double FisherExact::p_value() const { return round_as_r(p_two_sided_); }
double FisherExact::p_value_less() const { return round_as_r(p_less_); }
double FisherExact::p_value_greater() const { return round_as_r(p_greater_); }

}  // namespace senri::stats
