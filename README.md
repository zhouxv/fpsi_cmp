# Constrain Data Length: Enhanced Fuzzy PSI Through Secure Interval Membership Testing

Fuzzy PSI extends the conventional Private Set Intersection (PSI) by allowing intersections to be computed based on approximate matches. Concretely, a receiver and a sender each hold a set of d-dimensional points, the receiver learns which points from the sender’s set lie within a distance δ of its own points under a chosen distance metric.

We propose a fuzzy Private Set Intersection (PSI) protocol for the Lp distance with p ∈ [1,+∞], built upon efficient symmetric-key primitives. Our protocol achieves linear complexity in both the dimensionality d and logδ. The key building blocks are: (1) a novel fuzzy mapping protocol whose complexity is independent of δ, and (2) a multi-point secret-shared fuzzy matching scheme (mp-ssFMat). The latter can be viewed as a special case of fuzzy PSI in which the sender holds only a single item; we present an efficient construction of mp-ssFMat using symmetric-key techniques.

## 1. Project Structure

The following is the structure description of this project. We have provided some scripts to assist with building and running tests. Please grant execute permissions to these scripts, all of which are designed to be run from the project root directory.

```bash
├── fuzzyPSI/                          # Core FPSI protocol implementation directory
│   ├── andpair/                       # Implementation for AND pair operations
│   ├── main/                          # Main entry points and test files for various components
│   ├── mPeqt/                         # Multi-point equality test implementation
│   ├── ole/                           # Oblivious Linear Evaluation (OLE) protocol
│   ├── oprf_so/                       # Oblivious Pseudorandom Function
│   ├── psi/                           # Private Set Intersection (PSI) implementations for different metrics (L1, L2, Linfty)
│   └── utils/                         # Utility functions and helpers
├── CMakeLists.txt                     # CMake configuration file for building the project
├── Dockerfile                         # Docker image definition for containerized builds and runs
├── README.md                          # Project documentation and usage guide
├── shell_build_cmd.sh                 # Shell script to build the project executables
├── shell_install_dependencies.sh      # Shell script to install required system dependencies
├── shell_run_bench_fmap.sh            # Benchmark script to reproduce FMAP test results
└── shell_run_bench_fpsi.sh            # Benchmark script to reproduce FPSI test results
```

## 2. Prerequisites

**Supported OS :** `Ubuntu 20.04+` meet all project runtime specifications.

**Memory :** `100GB or above recommended`. Peak memory usage during image build ranges from 66GB to 80GB, which can cause out-of-memory issues on machines with less than 64 GB RAM. Please ensure sufficient memory is available before proceeding.

**Docker Version :** `Docker 28.3.3+`  (Earlier versions have not been validated).

**Dependencies :** local build must install the following dependencies :

```bash
gcc13           # required to ensure full *C++20* support
cmake >= 3.15   # cmake_minimum_required
git
python3
python3-pip
libgmp-dev
libspdlog-dev
libssl-dev
libmpfr-dev
libfmt-dev
libtool
nasm
```

## 3. Docker Build and Run benchmarks

### 3.1 Obtain Docker Image

*Option 1 : Build the image locally.*

```bash
# build docker image and run the docker container with the necessary capabilities
docker build -t fpsi_cmp:latest .
docker run -dit --name fpsi_cmp --cap-add=NET_ADMIN fpsi_cmp:latest
```

*Option 2 : Pull the latest image directly from Docker Hub public repository using the following command:*

```bash
# pull image from Docker Hub
docker pull blueobsidian/fpsi_cmp:latest
docker run -dit --name fpsi_cmp --cap-add=NET_ADMIN blueobsidian/fpsi_cmp:latest
```

### 3.2 Run Benchmark Scripts

```bash
# Reproduce fmap test results
./shell_run_bench_fmap.sh
# Reproduce fpsi test results
./shell_run_bench_fpsi.sh
```

## 4. Local Build and Run benchmarks

### 4.1 Install dependencies and build executables

```bash
# Install dependencies
./shell_install_dependencies.sh
# build the executables
./shell_build_cmd.sh
```

### 4.2 Run Benchmark Scripts in project root

```bash
# Reproduce fmap test results
./shell_run_bench_fmap.sh
# Reproduce fpsi test results
./shell_run_bench_fpsi.sh
```

## 5. Usage Guide for Executables

This section describes the usage of the executables located at `./build/fpsi` and `./build/fmap`.

### 5.1 Command-Line Options for `./build/fpsi`

| Flag | Meaning | Optional Values | Description |
|:----:|:--------|:----------------|:------------|
| **n** | Set Size (direct) | Positive integer, default: `4096` | Direct input set size |
| **nn** | Set Size (logarithm) | Positive integer, default: `12` | Input set size = 2^nn (overrides -n) |
| **dim** | Dimension | Positive integer, default: `6` | Dimension of the points |
| **delta** | Distance Threshold | Positive integer, default: `60` | Distance threshold δ for fuzzy matching |
| **metric** | Distance Metric | `0`: L∞ (default)<br/>`1`: L₁<br/>`2`: L₂ | Which p for L_p distance, 0 is infinite norm |
| **ip** | Server IP | IP address string, default: `"localhost"` | IP address for network communication |
| **port** | Server Port | Port number, default: `1212` | Port number for connections |
| **trait** | Number of Trials | Positive integer, default: `5` | Number of test runs for averaging results |
| **h/help** | Help Message | Flag (no value) | Print help message and usage |

### 5.2 Command-Line Options for `./build/fmap`

| Flag | Meaning | Optional Values | Description |
|:----:|:--------|:----------------|:------------|
| **n** | Set Size (direct) | Positive integer, default: `4096` | Direct input set size |
| **nn** | Set Size (logarithm) | Positive integer, default: `12` | Input set size = 2^nn (overrides -n) |
| **dim** | Dimension | Positive integer, default: `6` | Dimension of the points |
| **delta** | Distance Threshold | Positive integer, default: `60` | Distance threshold δ for fuzzy matching |
| **trait** | Number of Trials | Positive integer, default: `5` | Number of test runs for averaging results |
| **ip** | Server IP | IP address string, default: `"localhost"` | IP address for network communication |
| **port** | Server Port | Port number, default: `1213` | Port number for connections |
| **h/help** | Help Message | Flag (no value) | Print help message and usage |
