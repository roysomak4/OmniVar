# OmniVar

It is an ultra-fast variant caller for short-read DNA sequencing: SNVs, MNVs, indels, complex and
structural variants in single samples and tumour/normal pairs. It is built and validated for
hybrid-capture, paired-end targeted sequencing.

Its core is a rewrite of the robust and versatile **[VarDictJava v1.8.3](https://github.com/AstraZeneca-NGS/VarDictJava)**
variant caller (not actively maintained anymore) in C++23. The rewrite is focused on performance improvements, minimal external dependencies for easy deployment (no Java, R or Perl needed). Additional work is underway to extend the variant calling features, such as ITD and large deletion detection, improving signal to noise ratio, and UMI-aware variant calling. Future plans include adding the capability for variant calling on long-read sequencing data. Variant calling on amplicon data has not been included in this release. The app builds into one C++ static binary that can be run on any linux system. MacOS systems will require a compilation for the target architecture (AMD64 or ARM64/AppleSilicon).

[Claude Code](https://claude.com/product/claude-code) was used to automate tests for performance and accuracy benchmarks, sanitize builds (ASan/UBSan/TSan), security scans, CI, and generate documentation. All generated content was reviewed by the maintainer prior to commit and release.

## Highlights

- OmniVar v0.1.0 initial release is meant to be a near drop in replacement for **version 1.8.3 of VarDictJava**.

- **Validation**: On the [OncoSpan](https://horizondiscovery.com/en/reference-standards/products/oncospan-gdna) reference material, OmniVar finds all 187 expected variants between 2% and 48% VAF (ROI corresponding to the [CinCSeq](https://www.cincinnatichildrens.org/clinical-labs/our-labs/pathology/molecular-genomic-pathology-services/-/media/C99B1C995FC8428DA7080059BBF3AD70.ashx) NGS panel; benchmarked with [aardvark](https://github.com/PacificBiosciences/aardvark)).

- **Faster.** On the Oncospan RM (1200x average coverage). Runtime is 44s (OmniVar) versus 3m 16s (VarDictJava, not including R and Perl scripts steps) with 16 threads (`-th 16`).

- **One file to install.** A fully static x86-64 Linux binary (glibc, built on Debian trixie) that needs no libraries or dependencies on the target machine.

- **VCF by default.** VCF 4.2 output with the upstream FILTER tags. FILTER annotation never removes records. `

  --tsv` writes the upstream table.

- **New: `-p` reports rejected variants.** It also writes candidates the caller rejected, with the first rule they failed as FILTER (

  `REJ_*`). The calls themselves don't change.

## Install

Download `omnivar-v0.1.0-linux-x86_64.tar.gz` from the
[Releases](https://github.com/roysomak4/OmniVar/releases) page. It holds a fully static binary
that runs on any x86-64 Linux, with nothing else to install (the archive also holds `LICENSE`,
`THIRD_PARTY_LICENSES` and this README):

```
sha256sum -c SHA256SUMS
tar xzf omnivar-v0.1.0-linux-x86_64.tar.gz
sudo cp omnivar-v0.1.0-linux-x86_64/omnivar /usr/local/bin/
omnivar --version
```

## Usage

```
# single sample: BED with chrom, start, end, gene in columns 1-4
omnivar -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -f 0.01 -th 16 targets.bed > calls.vcf

# tumour/normal pair
omnivar -G ref.fa -b "tumour.bam|normal.bam" -N "tumour|normal" -c 1 -S 2 -E 3 -g 4 -f 0.01 -th 16 targets.bed > calls.vcf

# one region, no BED; -p also reports the candidates the caller rejected and why (FILTER REJ_*)
omnivar -G ref.fa -b sample.bam -N sample -p -R chr7:55181300-55181500 > region.vcf
```

The reference must be indexed (`.fai`) and the BAMs indexed. Output is VCF 4.2 with the upstream
FILTER tags; annotation never removes records. `--tsv` writes the upstream tab-separated table
instead. Options, letters and defaults are VarDictJava's; `omnivar -H` lists them all. Output is
the same for any thread count (`-th`).

To build from source, run the tests or benchmark omnivar against a truth set, see
[CONTRIBUTING.md](CONTRIBUTING.md).

## Building

Requirements: CMake ≥ 3.25, Ninja, and GCC ≥ 14, Clang ≥ 18 with libstdc++ 14, or Clang ≥ 20
with libc++ 20; zlib, libbz2, liblzma (libdeflate optional). htslib 1.22.1 is downloaded and
built automatically (pass `-DOMNIVAR_HTSLIB_TARBALL=…` to build offline).

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14
cmake --build build
ctest --test-dir build --output-on-failure
```

Static Linux binary (glibc, Debian trixie), as released:

```
docker build -f docker/debian-static.Dockerfile --target export --output type=local,dest=out .
```

## Licence

omnivar is released under the MIT License (`LICENSE`). The release binary is statically linked
with htslib, libdeflate, zlib, libbz2, liblzma, the GNU C Library and the GCC runtime; their
licences are in `THIRD_PARTY_LICENSES`, which ships with every release.

## Acknowledgements

OmniVar's variant-calling core is derived from **VarDictJava 1.8.3**
(https://github.com/AstraZeneca-NGS/VarDictJava), the Java implementation of VarDict by
Zhongwu Lai and colleagues, and reproduces its output near identically to be used as a drop in replacement. Its strand-bias and VCF
annotation reproduce VarDict's R and Perl scripts (`teststrandbias.R`, `testsomatic.R`,
`var2vcf_valid.pl`, `var2vcf_paired.pl`; the Perl scripts were written by Zhongwu Lai,
AstraZeneca, and are free software that may be changed and redistributed).
VarDictJava is distributed under the MIT License (See VarDictJava repository).
