# OmniVar
Fast and comprehensive DNA variant caller for Short Read sequencing data

`omnivar` is a single C++23 binary for SNV, MNV, indel, complex and structural variant
calling in single samples and tumour/normal pairs, with the strand-bias tests and VCF writing
built in. See `docs/design.md` for the design and `docs/validation.md` for checking it on your
own data.

## Building

Requirements: CMake ≥ 3.25, Ninja, GCC ≥ 14 / Clang ≥ 18 / AppleClang ≥ 16,
and zlib, libbz2, liblzma (libdeflate optional). htslib 1.22.1 is downloaded
and built automatically (pass `-DOMNIVAR_HTSLIB_TARBALL=…` to build offline).

```
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=g++-14
cmake --build build
ctest --test-dir build --output-on-failure
```

Static Linux binary (glibc, Debian trixie):

```
docker build -f docker/debian-static.Dockerfile --target export --output type=local,dest=out .
```

## Acknowledgements

omnivar's variant-calling core is derived from **VarDictJava 1.8.3**
(https://github.com/AstraZeneca-NGS/VarDictJava), the Java implementation of VarDict by
Zhongwu Lai and colleagues, and reproduces its output exactly. Its strand-bias and VCF
annotation reproduce VarDict's R and Perl scripts (`teststrandbias.R`, `testsomatic.R`,
`var2vcf_valid.pl`, `var2vcf_paired.pl`; the Perl scripts were written by Zhongwu Lai,
AstraZeneca, and are free software that may be changed and redistributed).
`docs/provenance.tsv` maps the C++ sources to their upstream origins. VarDictJava is
distributed under the MIT License:

```
The MIT License (MIT)

Copyright (c) 2019 AstraZeneca - NGS Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
