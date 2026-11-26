# Contributing to omnivar

omnivar started as a C++ rewrite of VarDictJava 1.8.3. Version 0.1.0 reproduces its output
exactly (see `docs/release-notes/v0.1.0.md`). Later versions improve on it, so differences from
VarDictJava are expected. What must not happen is an *unintended* change, and the regression
tests below exist to catch those. A deliberate change is recorded by updating the expected
output in the same commit.

## Requirements

- CMake ≥ 3.25, Ninja
- One of: GCC ≥ 14; Clang ≥ 18 with libstdc++ 14; Clang ≥ 20 with libc++ 20
- zlib, libbz2, liblzma (libdeflate optional). htslib 1.22.1 is downloaded (pinned by SHA-256) and
  built automatically; `-DOMNIVAR_HTSLIB_TARBALL=<file>` builds offline.
- Optional: bcftools (the regression tests also check that every VCF parses);
  [aardvark](https://github.com/PacificBiosciences/aardvark) for truth-set benchmarks;
  `clang-format-20` and `clang-tidy-20` for lint.

CMake options: `OMNIVAR_STATIC`, `OMNIVAR_BUILD_TESTS` (on), `OMNIVAR_WERROR`,
`OMNIVAR_SANITIZE` (ASan + UBSan), `OMNIVAR_LIBDEFLATE` (on), `OMNIVAR_HTSLIB_TARBALL`,
`OMNIVAR_DOWNLOAD_CACHE` (default `~/.cache/omnivar/downloads`).

## Build and test

```
# debug build with warnings as errors
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 -DOMNIVAR_WERROR=ON
cmake --build build
# unit, CLI, regression and invariant tests
ctest --test-dir build --output-on-failure -j 8

# Clang + libc++ with ASan and UBSan (also catches standard-library portability issues)
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=clang-20 -DCMAKE_CXX_COMPILER=clang++-20 \
      -DCMAKE_CXX_FLAGS=-stdlib=libc++ -DOMNIVAR_SANITIZE=ON
cmake --build build-asan && ctest --test-dir build-asan --output-on-failure -j 8

# static Linux binary, as released (Debian trixie, glibc): writes out/omnivar
docker build -f docker/debian-static.Dockerfile --target export --output type=local,dest=out .
```

Unit tests use doctest (`tests/third_party/doctest.h`) in `tests/unit/`, with their data in
`tests/unit/data/`; add new test files to `tests/unit/CMakeLists.txt`. Every new behaviour needs
a unit or regression test.

## Regression tests

All test tooling is C++: `tests/testkit/`, built as `build/tests/testkit/omnivar-testkit`
(`omnivar-testkit --help`). The repository has no Python.

`tests/regression/` holds 90 cases (single sample and tumour/normal: SNVs, indels, complex
variants, SVs, soft-clips, splicing, `-p`), with their BAMs and reference slices, which come
from the VarDictJava test suite:

- `cases/<case>.txt`: line 1 is the case (mode, reference, BAM, region, options); the rest is
  omnivar's expected `--tsv` output.
- `vcf/<case>.vcf`: the expected VCF of the same case.
- `synthetic/`: crafted TSV rows (many contigs, SVs, clusters, edge values) and the VCF that
  `omnivar --dev-tsv-to-vcf` writes for them with each set of VCF options.

```
TK=build/tests/testkit/omnivar-testkit
# compare with the expected output; --diff N shows differing rows, -k selects cases by name
$TK regress --binary build/omnivar [--part tsv,vcf,synthetic] [--mode Simple|Somatic] [-k substr] --diff 5
# after a deliberate change of output: rewrite the expectations, review, commit them with the change
$TK regress --binary build/omnivar --update && git diff --stat tests/regression
# invariants: -th N writes exactly what -th 1 writes; -p never changes the calls
$TK threads --binary build/omnivar -th 4
$TK reject-invariant --binary build/omnivar
```

## Benchmarks

```
# truth-set benchmark with aardvark: recall per variant type, optional call snapshot
$TK bench --binary build/omnivar --reference ref.fa --bed targets.bed -b sample.bam -N sample \
    --truth truth.vcf.gz [--snapshot FILE [--update-snapshot]] [--min-recall 1] -th 16 -- -f 0.01
# tumour/normal: score an existing VCF (PASS, STATUS Somatic, AF >= X)
$TK bench --vcf calls.vcf --reference ref.fa --bed eval.bed --truth truth.vcf.gz --somatic --min-af 0.05
```

CI runs the Oncospan benchmark: every truth variant must be recalled, and the call set must
equal `tests/validation/oncospan_snapshot.tsv` (update it with `--update-snapshot` when the
calls change on purpose).

## Code rules

- C++23 without compiler extensions, building with every compiler listed above. Allowed: all of
  C++20, plus `std::expected`, `std::print`/`std::println`, `std::ranges::to`,
  `string::contains`, `std::to_underlying`, `std::unreachable`, `std::byteswap`. Not yet:
  `std::mdspan`, `std::flat_map`, `std::generator`, modules, deducing `this`, `std::stacktrace`.
- No `std::regex` (use `util/regex_engine`), no floating-point `std::from_chars` (use
  `util::parse_double`). Numbers that appear in output go through `util/num_format`, which
  defines how they are written.
- Runtime dependencies: htslib and its compression libraries only.
- No global mutable state: an immutable `Config` is passed by const reference.
- Formatting: `.clang-format`. Static analysis: `.clang-tidy`, every enabled check clean.

```
git ls-files 'src/*.[ch]pp' 'tests/unit/*.cpp' 'tests/testkit/*.[ch]pp' | xargs clang-format-20 --dry-run --Werror   # fix: -i
cmake -S . -B build-lint -G Ninja -DCMAKE_C_COMPILER=clang-19 -DCMAKE_CXX_COMPILER=clang++-19 && cmake --build build-lint --target htslib_ext
run-clang-tidy-20 -p build-lint -quiet "$PWD/(src|tests/testkit)/"
```

## Branches, CI and releases

Changes go from a feature branch to `test`, then to `main`. CI
(`.github/workflows/ci.yml`) runs on the maintainer's self-hosted runner when `test` is pushed,
which needs write access. Pull requests do not run CI. Jobs run in this order:

1. **lint**: no Python, clang-format, clang-tidy.
2. **static**: the Docker static build, which runs all the tests inside, then the Oncospan
   accuracy gate.
3. **sanitizers**: ASan + UBSan, and TSan, on the tests and a subset of the panel. Their output
   must match the static binary's byte for byte.

`main` only receives commits that passed on `test`.

Releases: bump `project(VERSION)` in `CMakeLists.txt` and push a tag `v<VERSION>`. The release
workflow runs the whole CI with the Clang static analyzer, packages the static binary, and
drafts a GitHub release. Its notes include `docs/release-notes/<tag>.md` when that file exists.
