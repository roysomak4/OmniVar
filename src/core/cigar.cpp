#include "core/cigar.hpp"

namespace senri::core {

Cigar cigar_of(const io::BamRecord& record) {
    const bam1_t* b = record.raw();
    Cigar c;
    c.reserve(b->core.n_cigar);
    const uint32_t* cig = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        c.push_back(CigarElement{static_cast<int>(bam_cigar_oplen(cig[i])), bam_cigar_opchr(cig[i])});
    }
    return c;
}

std::string cigar_to_string(const Cigar& c) {
    if (c.empty()) {
        return "*";
    }
    std::string s;
    for (const auto& e : c) {
        s += std::to_string(e.length);
        s += e.op;
    }
    return s;
}

Cigar parse_cigar_text(std::string_view text) {
    Cigar c;
    if (text == "*") {
        return c;
    }
    int len = 0;
    for (char ch : text) {
        if (ch >= '0' && ch <= '9') {
            len = len * 10 + (ch - '0');
        } else {
            c.push_back(CigarElement{len, ch});
            len = 0;
        }
    }
    return c;
}

namespace {
template <class Pred>
int sum_if(const Cigar& c, Pred p) {
    int n = 0;
    for (const auto& e : c) {
        if (p(e.op)) {
            n += e.length;
        }
    }
    return n;
}
}  // namespace

int aligned_length(const Cigar& c) {
    return sum_if(c, [](char op) { return op == 'M' || op == 'D'; });
}
int soft_clipped_length(const Cigar& c) {
    return sum_if(c, [](char op) { return op == 'M' || op == 'I' || op == 'S'; });
}
int match_insertion_length(const Cigar& c) {
    return sum_if(c, [](char op) { return op == 'M' || op == 'I'; });
}
int insertion_deletion_length(const Cigar& c) {
    return sum_if(c, [](char op) { return op == 'I' || op == 'D'; });
}
int reference_length(const Cigar& c) {
    return sum_if(c, [](char op) { return consumes_reference(op); });
}

}  // namespace senri::core
