# Contributing to omnivar

omnivar's variant-calling core reproduces VarDictJava 1.8.3 exactly, quirks included (integer
division, rounding, tie-break order). Every change must keep that parity unless it adds a new
feature that is off by default. This guide covers building, testing and the checks a change
must pass.

## Requirements

- CMake ≥ 3.25, Ninja
- One of: GCC ≥ 14; Clang ≥ 18 with libstdc++ 14; Clang ≥ 20 with libc++ 20
- zlib, libbz2, liblzma (libdeflate optional). htslib 1.22.1 is downloaded (pinned by SHA-256) and
  built automatically; `-DOMNIVAR_HTSLIB_TARBALL=<file>` builds offline.
- For the parity tools only: a Java 21 runtime (the reference jar runs as the oracle), R
  (`Rscript`) and Perl (the reference strand-bias and VCF scripts), bcftools. For truth-set
  benchmarks: [aardvark](https://github.com/PacificBiosciences/aardvark).
- For lint: `clang-format-20`, `clang-tidy-20`.

CMake options: `OMNIVAR_STATIC`, `OMNIVAR_BUILD_TESTS` (on), `OMNIVAR_WERROR`,
`OMNIVAR_SANITIZE` (ASan + UBSan), `OMNIVAR_LIBDEFLATE` (on), `OMNIVAR_HTSLIB_TARBALL`,
`OMNIVAR_DOWNLOAD_CACHE` (default `~/.cache/omnivar/downloads`).

## Build and test

```
# debug build with warnings as errors
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 -DOMNIVAR_WERROR=ON
cmake --build build
# unit tests, CLI tests and the golden parity gates (no Java needed)
ctest --test-dir build --output-on-failure -j 8

# Clang + libc++ with ASan and UBSan (also catches standard-library portability issues)
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=clang-20 -DCMAKE_CXX_COMPILER=clang++-20 \
      -DCMAKE_CXX_FLAGS=-stdlib=libc++ -DOMNIVAR_SANITIZE=ON
cmake --build build-asan && ctest --test-dir build-asan --output-on-failure -j 8

# static Linux binary, as released (Debian trixie, glibc): writes out/omnivar
docker build -f docker/debian-static.Dockerfile --target export --output type=local,dest=out .
```

Unit tests use doctest (`tests/third_party/doctest.h`) in `tests/unit/`; add new test files to
`tests/unit/CMakeLists.txt`. Every new behaviour needs a unit or golden test.

## Parity and validation tools

All test tooling is C++: `tests/testkit/`, built as `build/tests/testkit/omnivar-testkit`, with one
subcommand per harness (`omnivar-testkit --help`). The repository has no Python.

```
TK=build/tests/testkit/omnivar-testkit
# golden cases (docs/golden, from the VarDictJava test suite): the same set of output rows as the
# reference, compared after sorting. --stage full | U (-U) | k0_U (-k 0 -U); --mode Simple|Somatic;
# --diff N shows differing rows. The 2 amplicon cases are skipped (amplicon mode is out of scope).
$TK golden --binary build/omnivar --stage full --diff 5
# VCF output against the reference R + Perl chain (checked-in expectations)
$TK vcf --binary build/omnivar -j 16 --stage full [--mode Somatic]
# -th N must give the same output as -th 1
$TK threads --binary build/omnivar -th 4
# -p must never change the calls
$TK reject-invariant --binary build/omnivar -j 16
# the reference jar on the golden cases (needs Java)
$TK oracle --check -j 16
# a whole panel against the jar and the R + Perl chain (your own data)
$TK panel --binary build/omnivar --out DIR --reference ref.fa --bed panel.bed -b sample.bam -N sample -- -f 0.01
# truth-set benchmark with aardvark: recall per variant type, optional call snapshot
$TK bench --binary build/omnivar --reference ref.fa --bed targets.bed -b sample.bam -N sample \
    --truth truth.vcf.gz [--snapshot FILE] [--min-recall 1] -th 16 -- -f 0.01
# tumour/normal: score an existing VCF (PASS, STATUS Somatic, AF >= X)
$TK bench --vcf calls.vcf --reference ref.fa --bed eval.bed --truth truth.vcf.gz --somatic --min-af 0.05
```

Never edit `vardictjavasrc/` (the reference source), `docs/golden/` (expected output) or
`tests/scripts/` (the reference R and Perl scripts): they are read-only references. Generated
expectations live under `tests/golden/`. When output differs from the reference, find the first
stage where it differs; do not tune thresholds to make a test pass.

## Code rules

- C++23 without compiler extensions, building with every compiler listed above. Allowed: all of
  C++20, plus `std::expected`, `std::print`/`std::println`, `std::ranges::to`,
  `string::contains`, `std::to_underlying`, `std::unreachable`, `std::byteswap`. Not yet:
  `std::mdspan`, `std::flat_map`, `std::generator`, modules, deducing `this`, `std::stacktrace`.
- No `std::regex` (use `util/regex_engine`), no floating-point `std::from_chars` (use
  `util::parse_double`). Numbers that appear in output go through `util/num_format`, so rounding
  matches the reference.
- Runtime dependencies: htslib and its compression libraries only.
- No global mutable state: an immutable `Config` is passed by const reference.
- The code does not name the reference implementation: comments say "upstream".
- Formatting: `.clang-format`. Static analysis: `.clang-tidy`, every enabled check clean (checks
  that would change upstream arithmetic or control flow are off; see the file).

```
git ls-files 'src/*.[ch]pp' 'tests/unit/*.cpp' 'tests/testkit/*.[ch]pp' | xargs clang-format-20 --dry-run --Werror   # fix: -i
cmake -S . -B build-lint -G Ninja -DCMAKE_C_COMPILER=clang-19 -DCMAKE_CXX_COMPILER=clang++-19 && cmake --build build-lint --target htslib_ext
run-clang-tidy-20 -p build-lint -quiet "$PWD/(src|tests/testkit)/"
```

## Branches, CI and releases

Changes go from a feature branch to `test`, then to `main`. CI
(`.github/workflows/ci.yml`) runs on the maintainer's self-hosted runner when `test` is pushed,
which needs write access. Pull requests do not run CI. Jobs run in this order:

1. **lint**: no Python, no upstream names in the code, clang-format, clang-tidy.
2. **static**: the Docker static build, which runs the unit and golden tests inside, then the
   accuracy gate. On the Oncospan reference material every truth variant must be recalled, and
   the call set must equal `tests/validation/oncospan_snapshot.tsv`.
3. **sanitizers**: ASan + UBSan, and TSan, on the tests and a subset of the panel. Their output
   must match the static binary's byte for byte.

`main` only receives commits that passed on `test`.

Releases: bump `project(VERSION)` in `CMakeLists.txt` and push a tag `v<VERSION>`. The release
workflow runs the whole CI with the Clang static analyzer, packages the static binary, and
drafts a GitHub release. Its notes include `docs/release-notes/<tag>.md` when that file exists.
