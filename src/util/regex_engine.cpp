#include "util/regex_engine.hpp"

#include <algorithm>
#include <bitset>
#include <limits>

namespace omnivar::util {

struct Regex::Node {
    enum class Kind { chars, any, start, end, group, alt, concat, repeat };
    Kind kind = Kind::concat;
    std::bitset<256> set;                     // chars
    std::vector<std::unique_ptr<Node>> kids;  // group (1), alt (n), concat (n), repeat (1)
    int group = -1;                           // capturing index, -1 for (?:...)
    int min = 0;
    int max = 0;  // -1 = unbounded
    bool greedy = true;
};

namespace {

using Node = Regex::Node;

std::bitset<256> class_of(char c) {
    std::bitset<256> s;
    auto add_range = [&](int a, int b) {
        for (int x = a; x <= b; ++x) {
            s.set(static_cast<std::size_t>(x));
        }
    };
    switch (c) {
        case 'd':
        case 'D': add_range('0', '9'); break;
        case 'w':
        case 'W':
            add_range('a', 'z');
            add_range('A', 'Z');
            add_range('0', '9');
            s.set('_');
            break;
        case 's':
        case 'S':
            // \s: [ \t\n\x0B\f\r]
            for (char x : {' ', '\t', '\n', '\x0B', '\f', '\r'}) {
                s.set(static_cast<unsigned char>(x));
            }
            break;
        default: break;
    }
    if (c == 'D' || c == 'W' || c == 'S') {
        s.flip();
    }
    return s;
}

class Parser {
public:
    explicit Parser(std::string_view p) : p_(p) {}

    std::expected<std::unique_ptr<Node>, std::string> parse() {
        auto n = alternation();
        if (!n) {
            return n;
        }
        if (i_ != p_.size()) {
            return std::unexpected(err(p_[i_] == ')' ? "Unmatched closing ')'" : "unexpected character"));
        }
        return n;
    }

    std::size_t groups() const { return groups_; }

private:
    std::string_view p_;
    std::size_t i_ = 0;
    std::size_t groups_ = 0;

    std::string err(std::string_view what) const {
        return std::string(what) + " near index " + std::to_string(i_) + " in regex \"" + std::string(p_) + "\"";
    }
    bool at_end() const { return i_ >= p_.size(); }
    char peek() const { return p_[i_]; }

    std::expected<std::unique_ptr<Node>, std::string> alternation() {
        auto alt = std::make_unique<Node>();
        alt->kind = Node::Kind::alt;
        while (true) {
            auto seq = sequence();
            if (!seq) {
                return seq;
            }
            alt->kids.push_back(std::move(*seq));
            if (!at_end() && peek() == '|') {
                ++i_;
                continue;
            }
            break;
        }
        if (alt->kids.size() == 1) {
            return std::move(alt->kids[0]);
        }
        return alt;
    }

    std::expected<std::unique_ptr<Node>, std::string> sequence() {
        auto seq = std::make_unique<Node>();
        seq->kind = Node::Kind::concat;
        while (!at_end() && peek() != '|' && peek() != ')') {
            auto atom_r = atom();
            if (!atom_r) {
                return atom_r;
            }
            auto q = quantifier(std::move(*atom_r));
            if (!q) {
                return q;
            }
            seq->kids.push_back(std::move(*q));
        }
        return seq;
    }

    std::expected<std::unique_ptr<Node>, std::string> quantifier(std::unique_ptr<Node> atom_node) {
        if (at_end()) {
            return atom_node;
        }
        int min = 0;
        int max = 0;
        const char c = peek();
        if (c == '*') {
            min = 0;
            max = -1;
            ++i_;
        } else if (c == '+') {
            min = 1;
            max = -1;
            ++i_;
        } else if (c == '?') {
            min = 0;
            max = 1;
            ++i_;
        } else if (c == '{') {
            std::size_t j = i_ + 1;
            auto number = [&](int& out) {
                std::size_t s = j;
                long v = 0;
                while (j < p_.size() && p_[j] >= '0' && p_[j] <= '9') {
                    v = v * 10 + (p_[j] - '0');
                    if (v > 1000000) {
                        return false;
                    }
                    ++j;
                }
                out = static_cast<int>(v);
                return j > s;
            };
            if (!number(min)) {
                return std::unexpected(err("Illegal repetition"));
            }
            max = min;
            if (j < p_.size() && p_[j] == ',') {
                ++j;
                if (j < p_.size() && p_[j] == '}') {
                    max = -1;
                } else if (!number(max)) {
                    return std::unexpected(err("Illegal repetition"));
                }
            }
            if (j >= p_.size() || p_[j] != '}') {
                return std::unexpected(err("Unclosed counted closure"));
            }
            if (max != -1 && max < min) {
                return std::unexpected(err("Illegal repetition range"));
            }
            i_ = j + 1;
        } else {
            return atom_node;
        }
        bool greedy = true;
        if (!at_end() && peek() == '?') {
            greedy = false;
            ++i_;
        } else if (!at_end() && peek() == '+') {
            return std::unexpected(err("possessive quantifiers are not supported"));
        }
        if (atom_node->kind == Node::Kind::start || atom_node->kind == Node::Kind::end) {
            return std::unexpected(err("quantified anchor is not supported"));
        }
        auto rep = std::make_unique<Node>();
        rep->kind = Node::Kind::repeat;
        rep->min = min;
        rep->max = max;
        rep->greedy = greedy;
        rep->kids.push_back(std::move(atom_node));
        // A quantifier directly after another is "Dangling meta character" upstream.
        if (!at_end() && (peek() == '*' || peek() == '+' || peek() == '?')) {
            return std::unexpected(err("Dangling meta character"));
        }
        return rep;
    }

    static std::unique_ptr<Node> chars_node(std::bitset<256> s) {
        auto n = std::make_unique<Node>();
        n->kind = Node::Kind::chars;
        n->set = s;
        return n;
    }

    // Escape after a backslash (i_ at the escaped char). Returns the char set.
    std::expected<std::bitset<256>, std::string> escape() {
        if (at_end()) {
            return std::unexpected(err("Unexpected internal error"));
        }
        const char c = p_[i_++];
        std::bitset<256> s;
        switch (c) {
            case 'd':
            case 'D':
            case 'w':
            case 'W':
            case 's':
            case 'S': return class_of(c);
            case 't': s.set('\t'); return s;
            case 'n': s.set('\n'); return s;
            case 'r': s.set('\r'); return s;
            case 'f': s.set('\f'); return s;
            case 'e': s.set(0x1B); return s;
            case 'a': s.set(0x07); return s;
            default: break;
        }
        const auto uc = static_cast<unsigned char>(c);
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')) {
            --i_;
            return std::unexpected(err(std::string("unsupported escape \\") + c));
        }
        s.set(uc);
        return s;
    }

    std::expected<std::unique_ptr<Node>, std::string> char_class() {
        // i_ is after '['
        bool negate = false;
        if (!at_end() && peek() == '^') {
            negate = true;
            ++i_;
        }
        std::bitset<256> s;
        bool first = true;
        while (true) {
            if (at_end()) {
                return std::unexpected(err("Unclosed character class"));
            }
            char c = peek();
            if (c == ']' && !first) {
                ++i_;
                break;
            }
            if (c == '[' || (c == '&' && i_ + 1 < p_.size() && p_[i_ + 1] == '&')) {
                return std::unexpected(err("nested classes and intersections are not supported"));
            }
            first = false;
            std::bitset<256> item;
            int lo = -1;
            if (c == '\\') {
                ++i_;
                auto e = escape();
                if (!e) {
                    return std::unexpected(e.error());
                }
                item = *e;
                if (item.count() == 1) {
                    for (int x = 0; x < 256; ++x) {
                        if (item.test(static_cast<std::size_t>(x))) {
                            lo = x;
                        }
                    }
                }
            } else {
                ++i_;
                lo = static_cast<unsigned char>(c);
                item.set(static_cast<std::size_t>(lo));
            }
            // Range a-b (a '-' before ']' is a literal).
            if (lo >= 0 && i_ + 1 < p_.size() && p_[i_] == '-' && p_[i_ + 1] != ']') {
                ++i_;
                int hi;
                if (p_[i_] == '\\') {
                    ++i_;
                    auto e = escape();
                    if (!e || e->count() != 1) {
                        return std::unexpected(err("Illegal character range"));
                    }
                    hi = 0;
                    for (int x = 0; x < 256; ++x) {
                        if (e->test(static_cast<std::size_t>(x))) {
                            hi = x;
                        }
                    }
                } else {
                    hi = static_cast<unsigned char>(p_[i_++]);
                }
                if (hi < lo) {
                    return std::unexpected(err("Illegal character range"));
                }
                for (int x = lo; x <= hi; ++x) {
                    item.set(static_cast<std::size_t>(x));
                }
            }
            s |= item;
        }
        if (negate) {
            s.flip();
        }
        return chars_node(s);
    }

    std::expected<std::unique_ptr<Node>, std::string> atom() {
        const char c = p_[i_];
        switch (c) {
            case '(': {
                ++i_;
                int index = -1;
                if (!at_end() && peek() == '?') {
                    if (i_ + 1 < p_.size() && p_[i_ + 1] == ':') {
                        i_ += 2;
                    } else {
                        return std::unexpected(err("only (?:...) groups are supported"));
                    }
                } else {
                    index = static_cast<int>(++groups_);
                }
                auto inner = alternation();
                if (!inner) {
                    return inner;
                }
                if (at_end() || peek() != ')') {
                    return std::unexpected(err("Unclosed group"));
                }
                ++i_;
                auto g = std::make_unique<Node>();
                g->kind = Node::Kind::group;
                g->group = index;
                g->kids.push_back(std::move(*inner));
                return g;
            }
            case '[': ++i_; return char_class();
            case '.': {
                ++i_;
                std::bitset<256> s;
                s.set();
                s.reset('\n');
                s.reset('\r');
                return chars_node(s);
            }
            case '^': {
                ++i_;
                auto n = std::make_unique<Node>();
                n->kind = Node::Kind::start;
                return n;
            }
            case '$': {
                ++i_;
                auto n = std::make_unique<Node>();
                n->kind = Node::Kind::end;
                return n;
            }
            case '\\': {
                ++i_;
                auto e = escape();
                if (!e) {
                    return std::unexpected(e.error());
                }
                return chars_node(*e);
            }
            case '*':
            case '+':
            case '?': return std::unexpected(err("Dangling meta character '" + std::string(1, c) + "'"));
            default: {
                ++i_;
                std::bitset<256> s;
                s.set(static_cast<unsigned char>(c));
                return chars_node(s);
            }
        }
    }
};

using Caps = std::vector<std::optional<std::pair<std::size_t, std::size_t>>>;

// A non-owning reference to a continuation `bool(std::size_t)`. The matcher's continuations
// are lambdas that live on the stack for the duration of the call that receives them, so a
// reference is enough, and unlike std::function it never allocates.
class Cont {
public:
    template <class F>
    Cont(const F& f) : obj_(&f), call_([](const void* o, std::size_t i) { return (*static_cast<const F*>(o))(i); }) {}
    bool operator()(std::size_t i) const { return call_(obj_, i); }

private:
    const void* obj_;
    bool (*call_)(const void*, std::size_t);
};

class Matcher {
public:
    Matcher(std::string_view t, Caps& caps) : t_(t), caps_(caps) {}

    bool m(const Node* n, std::size_t i, const Cont& k) {
        switch (n->kind) {
            case Node::Kind::chars:
                if (i < t_.size() && n->set.test(static_cast<unsigned char>(t_[i]))) {
                    return k(i + 1);
                }
                return false;
            case Node::Kind::any: return i < t_.size() && k(i + 1);
            case Node::Kind::start: return i == 0 && k(i);
            case Node::Kind::end: {
                // $ without MULTILINE: end of input, or before a final line terminator.
                const std::size_t len = t_.size();
                bool ok = i == len;
                if (!ok && i + 1 == len && (t_[i] == '\n' || t_[i] == '\r')) {
                    ok = true;
                }
                if (!ok && i + 2 == len && t_[i] == '\r' && t_[i + 1] == '\n') {
                    ok = true;
                }
                return ok && k(i);
            }
            case Node::Kind::group: {
                if (n->group < 0) {
                    return m(n->kids[0].get(), i, k);
                }
                const auto g = static_cast<std::size_t>(n->group);
                return m(n->kids[0].get(), i, [this, &k, i, g](std::size_t j) {
                    auto old = caps_[g];
                    caps_[g] = std::pair{i, j};
                    if (k(j)) {
                        return true;
                    }
                    caps_[g] = old;
                    return false;
                });
            }
            case Node::Kind::alt:
                for (const auto& kid : n->kids) {
                    if (m(kid.get(), i, k)) {
                        return true;
                    }
                }
                return false;
            case Node::Kind::concat: return seq(n, 0, i, k);
            case Node::Kind::repeat: return rep(n, 0, i, k);
        }
        return false;
    }

private:
    std::string_view t_;
    Caps& caps_;

    bool seq(const Node* n, std::size_t idx, std::size_t i, const Cont& k) {
        if (idx == n->kids.size()) {
            return k(i);
        }
        return m(n->kids[idx].get(), i, [this, n, idx, &k](std::size_t j) { return seq(n, idx + 1, j, k); });
    }

    bool rep(const Node* n, int count, std::size_t i, const Cont& k) {
        const bool can_more = n->max < 0 || count < n->max;
        auto more = [&]() {
            return can_more && m(n->kids[0].get(), i, [this, n, count, i, &k](std::size_t j) {
                       // As in the reference regex engine's Loop: an iteration that matched empty does not repeat.
                       if (j == i && count >= n->min) {
                           return false;
                       }
                       return rep(n, count + 1, j, k);
                   });
        };
        if (n->greedy) {
            if (more()) {
                return true;
            }
            return count >= n->min && k(i);
        }
        if (count >= n->min && k(i)) {
            return true;
        }
        return more();
    }
};

}  // namespace

Regex::Regex(Regex&&) noexcept = default;
Regex& Regex::operator=(Regex&&) noexcept = default;
Regex::~Regex() = default;

std::optional<std::string_view> Regex::Match::group(std::string_view text, std::size_t i) const {
    if (i >= groups.size() || !groups[i]) {
        return std::nullopt;
    }
    return text.substr(groups[i]->first, groups[i]->second - groups[i]->first);
}

namespace {

// True when every match must begin at index 0 (the pattern starts with ^ outside any alternation).
bool anchored_at_start(const Node* n) {
    switch (n->kind) {
        case Node::Kind::start: return true;
        case Node::Kind::concat: return !n->kids.empty() && anchored_at_start(n->kids[0].get());
        case Node::Kind::group: return anchored_at_start(n->kids[0].get());
        default: return false;
    }
}

// Characters every match must contain (single-character classes in mandatory positions).
void required_chars(const Node* n, std::bitset<256>& out) {
    switch (n->kind) {
        case Node::Kind::chars:
            if (n->set.count() == 1) {
                for (std::size_t c = 0; c < 256; ++c) {
                    if (n->set.test(c)) {
                        out.set(c);
                    }
                }
            }
            return;
        case Node::Kind::group: required_chars(n->kids[0].get(), out); return;
        case Node::Kind::concat:
            for (const auto& kid : n->kids) {
                required_chars(kid.get(), out);
            }
            return;
        case Node::Kind::repeat:
            if (n->min >= 1) {
                required_chars(n->kids[0].get(), out);
            }
            return;
        default:  // any, anchors, alternation: nothing certain
            return;
    }
}

}  // namespace

std::expected<Regex, std::string> Regex::compile(std::string_view pattern) {
    Parser p(pattern);
    auto root = p.parse();
    if (!root) {
        return std::unexpected(root.error());
    }
    Regex re;
    re.pattern_ = std::string(pattern);
    re.groups_ = p.groups();
    re.root_ = std::move(*root);
    re.anchored_ = anchored_at_start(re.root_.get());
    required_chars(re.root_.get(), re.required_);
    return re;
}

std::optional<Regex::Match> Regex::find(std::string_view text, std::size_t from) const {
    // Fast rejections that cannot change the result: an anchored pattern can only match at 0,
    // and a match needs every required character somewhere after `from`.
    const bool anchored = anchored_;
    if (anchored && from > 0) {
        return std::nullopt;
    }
    const std::bitset<256>& required = required_;
    if (required.any()) {
        std::bitset<256> present;
        for (std::size_t i = from; i < text.size(); ++i) {
            present.set(static_cast<unsigned char>(text[i]));
        }
        if ((required & ~present).any()) {
            return std::nullopt;
        }
    }
    const std::size_t last_start = anchored ? from : text.size();
    Caps caps(groups_ + 1);
    for (std::size_t start = from; start <= last_start; ++start) {
        std::fill(caps.begin(), caps.end(), std::nullopt);
        Matcher matcher(text, caps);
        std::size_t end = 0;
        if (matcher.m(root_.get(), start, [&](std::size_t j) {
                end = j;
                return true;
            })) {
            Match mt;
            mt.begin = start;
            mt.end = end;
            caps[0] = std::pair{start, end};
            mt.groups = caps;
            return mt;
        }
    }
    return std::nullopt;
}

std::vector<std::string_view> regex_split(std::string_view text, const Regex& re) {
    std::vector<std::string_view> out;
    std::size_t index = 0;
    std::size_t search = 0;
    bool matched = false;
    while (search <= text.size()) {
        auto mt = re.find(text, search);
        if (!mt) {
            break;
        }
        if (mt->end == 0 && mt->begin == 0) {
            // No empty leading substring for a zero-width match at the beginning (JDK 8+).
            search = 1;
            continue;
        }
        matched = true;
        out.push_back(text.substr(index, mt->begin - index));
        index = mt->end;
        // Matcher.find after a zero-length match advances one position.
        search = mt->end == mt->begin ? mt->end + 1 : mt->end;
    }
    if (!matched) {
        return {text};
    }
    out.push_back(text.substr(index));
    while (!out.empty() && out.back().empty()) {
        out.pop_back();
    }
    return out;
}

std::expected<std::vector<std::string_view>, std::string> regex_split(std::string_view text, std::string_view regex) {
    auto re = Regex::compile(regex);
    if (!re) {
        return std::unexpected(re.error());
    }
    return regex_split(text, *re);
}

}  // namespace omnivar::util
