// Local realignment (-k 1): CigarModifier against the jar.

#include <cstdio>
#include <doctest.h>
#include <fcntl.h>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "cli/config.hpp"
#include "core/cigar_modifier.hpp"
#include "core/scope.hpp"

namespace {

using namespace omnivar;

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> f;
    std::size_t start = 0;
    while (true) {
        const auto t = line.find('\t', start);
        f.push_back(line.substr(start, t == std::string::npos ? std::string::npos : t - start));
        if (t == std::string::npos) {
            return f;
        }
        start = t + 1;
    }
}

io::Reference make_reference(int start, const std::string& seq) {
    io::Reference ref;
    for (std::size_t i = 0; i < seq.size(); ++i) {
        ref.reference_sequences[start + static_cast<int>(i)] = seq[i];
    }
    for (std::size_t i = 0; i + 12 <= seq.size(); ++i) {
        ref.seed[seq.substr(i, 12)].push_back(start + static_cast<int>(i));
    }
    return ref;
}

// Silences stderr (printExceptionAndContinue messages) for the scope.
struct QuietStderr {
    QuietStderr() {
        std::fflush(stderr);
        saved = ::dup(2);
        const int null_fd = ::open("/dev/null", O_WRONLY);
        if (null_fd >= 0) {
            ::dup2(null_fd, 2);
            ::close(null_fd);
        }
    }
    ~QuietStderr() {
        std::fflush(stderr);
        if (saved >= 0) {
            ::dup2(saved, 2);
            ::close(saved);
        }
    }
    QuietStderr(const QuietStderr&) = delete;
    QuietStderr& operator=(const QuietStderr&) = delete;
    int saved = -1;
};

}  // namespace

TEST_CASE("CigarModifier: CigarModifierTest.threeIndels") {
    struct Case {
        const char* initial;
        const char* expected;
        int indel;
        int position;
        const char* query;
        const char* quality;
        const char* ref;
    };
    const Case cases[] = {
        {"15M1I1M1I2M1I27M", "27M3I18M", 3, 1, "GTGAAAACTGGAAACAAAAAAAAAAAAAAACAGCCTCTTCTTGGAAAG",
         "C@@FFFFFHHGDHIJIIJGJJJIIGJJJIHD6=3@A;;>;(;ACCCDCB", "GTGAAAACTGGAAACAAAAAAAAAAAACAGCCTCTTCTTGG"},
        {"14M3I2M2I2M2D26M1834N23M", "18M2D5I26M1834N23M", 7, 1,
         "GACTCTTTTGTAGGTGCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGCTACACTGTGTGTCATCTGC",
         "AAAAAEEEEEEEE6EEEEAEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEAE",
         "GACTCTTTTGTAGGTGCTNNTCTTAAGAAATTTAAATTTCT"},
        {"2S2M2I2M2D26M1834N10M", "6S28M1834N10M", 4, 1, "GCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGCTACACT",
         "AAAA6/EAEEEEEEEEEAEEAEE/6AEEE6EEA/EE/EEAA/EE", "GGTGCTGATCTTAAGAAATTTAAATTTCTAAAAGNNNNNNN"},
        {"34M3I2M2I2M2D26M1834N4M", "44M3I22M1834N4M", 7, 1,
         "GACTCCTTTTGTACCTCTGTGACTCTTTTGTAGGTGCTGATCTTCTTAAGAAATTTAAATTTCTAAAAGGTGC",
         "EEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAAAAA",
         "GACTCCTTTTGTACCTCTGTGACTCTTTTGTAGGTGCTGATCTTAAGAAATTTAAATTTCTAAAAGNNNNNNNNNNNNNNNN"},
        {"11M1I2M1I8M1I51M", "21M3I51M", 3, 1,
         "GTTTTTTGTTTCCTTTTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTT",
         "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEA<AAEEEAEEEEEE6EE/EEEEE6AEEEAEEEEEA<EEEEEEEEEA",
         "GTTTTTTGTTTCCTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTTNNNNNNNNNN"},
        {"10M2I2M1I8M1I51M", "21M3I50M", 4, 1,
         "GTTTTTTGTTTCCTTTTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTT",
         "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEA<AAEEEAEEEEEE6EE/EEEEE6AEEEAEEEEEA<EEEEEEEEEA",
         "GTTTTTTGTTTCCTTTTTTTTAACCATCTGATACTAAGAAGATGAATTTGCACAGATTTCTCTGCATAATTTNNNNNNNNNN"},
        {"16M1I1M1I2M1D3M1I50M", "23M2I50M", 4, 1,
         "CCTGTTCAAAGAAACCAAGAGGGAAAAAAAAAAGTGCAGAAGAAAAGTAATTTGCACTTAAAAATGAAAACAAAC",
         "EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/AEEEEEEEEEEAAAAA",
         "CCTGTTCAAAGAAACCAAGAGGGAAAAAAAAGTGCAGAAGAAAAGTAATTTGCACTTAAAAATGAAAACAAACCACTTTTCTC"},
    };
    const Config conf;
    for (const Case& c : cases) {
        io::Reference ref;  // the upstream test leaves Reference.seed empty
        const std::string rs = c.ref;
        for (std::size_t i = 0; i < rs.size(); ++i) {
            ref.reference_sequences[1 + static_cast<int>(i)] = rs[i];
        }
        core::ExceptionCounter exceptions;
        const Region region{"chr1", 1, 100, ""};
        auto r =
            core::modify_cigar(c.position, c.initial, c.query, c.quality, ref, c.indel, region, 150, conf, exceptions);
        CHECK_MESSAGE(r.cigar == c.expected, c.initial);
    }
}

TEST_CASE("CigarModifier matches the jar (tests/golden/cigar_modifier_vectors.tsv)") {
    std::ifstream in(OMNIVAR_SOURCE_DIR "/tests/golden/cigar_modifier_vectors.tsv");
    REQUIRE(in.good());
    const Config conf;
    const Region region{"chr1", 1, 500, ""};
    std::string line;
    int n = 0;
    int bad = 0;
    QuietStderr quiet;
    while (std::getline(in, line)) {
        const auto f = split_tabs(line);
        REQUIRE(f.size() == 13);
        const io::Reference ref = make_reference(std::stoi(f[4]), f[5]);
        core::ExceptionCounter exceptions;
        auto r = core::modify_cigar(std::stoi(f[0]), f[1], f[2], f[3], ref, std::stoi(f[6]), region, std::stoi(f[7]),
                                    conf, exceptions);
        const bool ok = r.position == std::stoi(f[8]) && r.cigar == f[9] && r.query_sequence == f[10] &&
                        r.query_quality == f[11] && exceptions.count == std::stoi(f[12]);
        if (!ok && ++bad <= 5) {
            std::ostringstream msg;
            msg << "pos " << f[0] << " " << f[1] << ": got " << r.position << " " << r.cigar << " exc "
                << exceptions.count << ", want " << f[8] << " " << f[9] << " exc " << f[12];
            FAIL_CHECK(msg.str());
        }
        ++n;
    }
    CHECK(bad == 0);
    CHECK(n > 1000);
}
