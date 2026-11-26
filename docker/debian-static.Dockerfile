# Fully static Linux release build of omnivar on glibc (Debian trixie, GCC 14).
#
#   docker build -f docker/debian-static.Dockerfile -t omnivar-static .
#   docker run --rm omnivar-static --version
#
# Extract the binary without running a container:
#   docker build -f docker/debian-static.Dockerfile --target export --output type=local,dest=out .
#   -> out/omnivar
#
# Downloads are cached across builds with BuildKit cache mounts (Docker 23+ uses BuildKit by
# default): apt packages, and the htslib tarball (OMNIVAR_DOWNLOAD_CACHE, see cmake/htslib.cmake).
# `docker build --no-cache` also starts these caches afresh (BuildKit 0.34), so prefer plain
# rebuilds. Clear them with `docker builder prune --filter type=exec.cachemount`.

FROM debian:trixie AS build

# Keep downloaded .debs (the image deletes them by default) so the cache mount is useful.
RUN rm -f /etc/apt/apt.conf.d/docker-clean \
 && echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache
RUN --mount=type=cache,id=omnivar-apt-cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,id=omnivar-apt-lists,target=/var/lib/apt/lists,sharing=locked \
    apt-get update \
 && apt-get install -y --no-install-recommends \
        g++ gcc make cmake ninja-build ca-certificates curl bzip2 file \
        zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev

WORKDIR /src
COPY CMakeLists.txt ./
COPY cmake ./cmake
COPY src ./src
COPY tests ./tests

ENV OMNIVAR_DOWNLOAD_CACHE=/var/cache/omnivar-downloads
RUN --mount=type=cache,id=omnivar-downloads,target=/var/cache/omnivar-downloads \
    cmake -S . -B /build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
        -DOMNIVAR_STATIC=ON -DOMNIVAR_WERROR=ON \
 && cmake --build /build \
 && ctest --test-dir /build --output-on-failure \
 && strip /build/omnivar \
 && file /build/omnivar \
 && (ldd /build/omnivar 2>&1 || true) | grep -q "not a dynamic executable" \
 && /build/omnivar --version

# Only the binary.
FROM scratch AS export
COPY --from=build /build/omnivar /omnivar

# Default target: minimal runnable image.
FROM scratch
COPY --from=build /build/omnivar /omnivar
ENTRYPOINT ["/omnivar"]
