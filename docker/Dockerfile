FROM --platform=linux/arm/v7 debian:buster

# Update sources to use archived Debian Buster repos (if needed for older packages)
RUN sed -i 's/deb.debian.org/archive.debian.org/g' /etc/apt/sources.list && \
    sed -i 's/security.debian.org/archive.debian.org/g' /etc/apt/sources.list

# Install basic build tools and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libboost-all-dev \
    libgmp-dev \
    libmpfr-dev \
    libeigen3-dev \
    wget \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Build and install CMake 3.28 from source (required for CGAL 6.1)
RUN wget https://github.com/Kitware/CMake/releases/download/v3.28.0/cmake-3.28.0.tar.gz && \
    tar -xzf cmake-3.28.0.tar.gz && \
    cd cmake-3.28.0 && \
    ./bootstrap -- -DOPENSSL_ROOT_DIR=/usr && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf cmake-3.28.0*

# Build and install Boost 1.84 from source (required for CGAL)
RUN wget https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.bz2 && \
    tar -xjf boost_1_84_0.tar.bz2 && \
    cd boost_1_84_0 && \
    ./bootstrap.sh --prefix=/usr/local && \
    ./b2 install && \
    cd .. && rm -rf boost_1_84_0*

# Ensure Boost is found
ENV BOOST_ROOT=/usr/local

# Optional: Build GMP and MPFR statically from source for static linking
# Uncomment if needed
# RUN wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz \
#     && tar -xf gmp-6.3.0.tar.xz && cd gmp-6.3.0 \
#     && ./configure --enable-static --disable-shared --prefix=/usr/local \
#     && make -j$(nproc) && make install && cd .. && rm -rf gmp-6.3.0* \
#     && wget https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz \
#     && tar -xf mpfr-4.2.1.tar.xz && cd mpfr-4.2.1 \
#     && ./configure --enable-static --disable-shared --prefix=/usr/local --with-gmp=/usr/local \
#     && make -j$(nproc) && make install && cd .. && rm -rf mpfr-4.2.1*

# Set working directory
WORKDIR /build

# Default command: bash for interactive building
CMD ["bash"]
