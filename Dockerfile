FROM ubuntu:24.04

WORKDIR /home

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    vim \
    git \
    python3 \
    python3-pip \
    cmake \
    libgmp-dev \
    libspdlog-dev \
    libtool \
    nasm \
    libssl-dev \
    libmpfr-dev \
    iproute2 \
    net-tools \
    curl\
    jq  && \
    # install tcconfig for network interface configuration
    curl -sSL https://raw.githubusercontent.com/thombashi/tcconfig/master/scripts/installer.sh | bash

# Install thirdparty dependencies
COPY ./shell_install_dependencies.sh ./

RUN chmod +x ./*.sh && \
    ./shell_install_dependencies.sh

RUN rm -rf ./thirdparty

#  Copy sourcode files and build executable file
COPY ./shell_build_cmd.sh \
    ./CMakeLists.txt \
    ./
COPY ./fuzzyPSI/ ./fuzzyPSI/

# Build executable file
RUN chmod +x ./*.sh && \
    ./shell_build_cmd.sh

# # copy other files
COPY ./README.md \
    # ./shell_run_bench_fmap.sh \
    ./shell_run_bench_fpsi.sh \
    ./

RUN chmod +x ./*.sh

