FROM ubuntu:24.04

WORKDIR /workspace

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    autoconf \
    automake \
    libgmp-dev \
    libspdlog-dev \
    libtool \   
    libssl-dev \
    libmpfr-dev \
    nasm \
    python3 \
    python3-pip \
    vim \
    git \
    iproute2 \
    net-tools \
    curl \
    jq && \
    rm -rf /var/lib/apt/lists/*

# Install tcconfig for network interface configuration
RUN curl -fsSL \
    https://raw.githubusercontent.com/thombashi/tcconfig/master/scripts/installer.sh \
    | bash

# Install thirdparty dependencies
COPY --chmod=755 shell_install_dependencies.sh /workspace

RUN chmod +x shell_install_dependencies.sh && \
    ./shell_install_dependencies.sh && \
    rm -rf /workspace/thirdparty

# Copy source code and build files
COPY shell_build_cmd.sh \
    CMakeLists.txt \
    /workspace/
COPY ./fuzzyPSI/ /workspace/fuzzyPSI/

# Build executable file
RUN chmod +x shell_build_cmd.sh && \
    /workspace/shell_build_cmd.sh

# Copy runtime and benchmark files
COPY README.md \
    shell_run_bench_fpsi.sh \
    /workspace/

RUN chmod +x shell_run_bench_fpsi.sh