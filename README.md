# Modulo Interval Membership Testing: Towards Practical Fuzzy Private Set Intersection

This artifact implements the fuzzy PSI protocols evaluated in Section 8 of
the paper, with online costs in Table 2 and overall costs (offline + online)
in Table 3 and Appendix D. A receiver learns which sender points are within a
distance threshold of its own points under the Linf, L1, or L2 metric.

**Badges requested:** Artifact Available, Artifact Functional, and Artifact
Results Reproduced. These are evaluation targets, not awarded badges.
The final evaluated artifact is intended for permanent public archival with a
DOI; the archive link will be added when available.

## 1. Project Structure

The protocol implementation is located in `fuzzyPSI/`, with a benchmark entry
point, metric-specific fuzzy PSI protocols, and supporting cryptographic
components. The scripts in the project root handle dependency installation,
building, network configuration, and benchmark execution, as outlined below.

```text
.
├── fuzzyPSI/                      Protocol implementation
│   ├── main/
│   │   └── main.cpp               Benchmark entry: inputs, both parties, measurements, and CSV output
│   ├── psi/
│   │   ├── psi.h                  Sender/receiver interfaces and protocol parameters
│   │   ├── Defines.h              Shared types and definitions
│   │   ├── psiLinfty.cpp          Linf fuzzy PSI: offline and online phases
│   │   ├── psiL1.cpp              L1 fuzzy PSI: offline and online phases
│   │   ├── psiL2.cpp              L2 fuzzy PSI: offline and online phases
│   │   ├── fmap.h / fmap.cpp      Proximity-collision mapping for candidate matches
│   │   ├── mImt.h / mImt.cpp      Modulo interval membership testing
│   │   └── l2_offline_cache.h     One-time L2 offline-material storage and consumption
│   ├── mPeqt/                    Private equality testing and preprocessing/OT support
│   ├── andpair/                  Boolean triples and AND-pair correlations
│   └── ole/                      OLE and modular VOLE support for arithmetic operations
├── CMakeLists.txt                 Build configuration
├── Dockerfile                     Container build
├── shell_install_dependencies.sh  Dependency installation
├── shell_build_cmd.sh             Project build
├── shell_config_network.sh        Network configuration
└── shell_run_bench_fpsi.sh        Experiments and CSV output
```

## 2. Requirements

### 2.1 Hardware

An x86-64 CPU is required; no GPU is used. The paper's experiments used an
Intel Xeon Gold 6338 and 128 GB RAM. This describes the paper's machine,
not a measured minimum memory requirement. Peak memory requirements for the
minimal, quick, and full experiments remain to be documented.

### 2.2 Software

The Dockerfile uses Ubuntu 24.04. Local builds require GCC 13 or a compatible
C++20 compiler, CMake 3.15 or later, Python 3.9 or later, and the system packages
below. Network configuration uses tcconfig, iproute2, and the Linux
`sch_netem` module.

On Ubuntu 24.04:

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake autoconf automake ca-certificates curl git iproute2 jq \
  libfmt-dev libgmp-dev libmpfr-dev libspdlog-dev libssl-dev libtool nasm \
  python3 python3-pip python3-venv
```

Network changes require root/sudo on a local machine or the `NET_ADMIN`
capability inside Docker.

## 3. Installation

### 3.1 Build Docker Image
Build the image and open a shell inside it:

```bash
docker build -t fpsi_cmp_artifact:latest .
docker run -d --cap-add=NET_ADMIN \
  --name fpsi_cmp_ours \
  fpsi_cmp_artifact:latest \
  sleep infinity
```

**Both baseline images can also be built from source; prebuilt images are
not required.** For Docker build and container startup instructions, see
the **Build and Run with Docker** section in the corresponding repository's
README:

- [Exp11: so-OPPRF-based fuzzy PSI — Docker build instructions](https://github.com/zhouxv/fpsi_ssoprf/blob/fpsi-cmp_artifact_20260919/README.md#build-and-run-with-docker)
- [Exp12: da-ROT-based fuzzy PSI — Docker build instructions](https://github.com/zhouxv/fpsi_daOT/blob/fpsi-cmp_artifact_20260916/README.md#build-and-run-with-docker)

Allow approximately **25–35 minutes per image** to build Ours, Exp11, or
Exp12 from source, including dependency downloads and compilation. This is a
rough planning estimate based on our build logs, not a measured minimum or
upper bound. Actual build time depends on hardware, download speed, and
Docker's build cache; slower downloads may take longer.

### 3.2 Prebuilt Docker images

Using the prebuilt Docker images is the recommended way to reproduce the
experiments. The source repositories and corresponding Docker images are:

| Implementation | Source repository | Docker image |
|---|---|---|
| Ours | `https://github.com/zhouxv/fpsi_cmp` | `blueobsidian/fpsi_cmp_artifact:latest` |
| Exp11: so-OPPRF-based fuzzy PSI [\[1\]](#ref-1) | `https://github.com/zhouxv/fpsi_ssoprf/tree/fpsi-cmp_artifact_20260919` | `blueobsidian/fpsi_cmp_artifact_exp11:latest` |
| Exp12: da-ROT-based fuzzy PSI [\[2\]](#ref-2) | `https://github.com/zhouxv/fpsi_daOT/tree/fpsi-cmp_artifact_20260916` | `blueobsidian/fpsi_cmp_artifact_exp12:latest` |

Pull the three images:

```bash
docker pull blueobsidian/fpsi_cmp_artifact:latest
docker pull blueobsidian/fpsi_cmp_artifact_exp11:latest
docker pull blueobsidian/fpsi_cmp_artifact_exp12:latest
```

Start the three containers:

```bash
docker run -d --cap-add=NET_ADMIN \
  --name fpsi_cmp_ours \
  blueobsidian/fpsi_cmp_artifact:latest \
  sleep infinity

docker run -d --cap-add=NET_ADMIN \
  --name fpsi_cmp_exp11 \
  blueobsidian/fpsi_cmp_artifact_exp11:latest \
  sleep infinity

docker run -d --cap-add=NET_ADMIN \
  --name fpsi_cmp_exp12 \
  blueobsidian/fpsi_cmp_artifact_exp12:latest \
  sleep infinity
```

The containers are kept running so that the experiments in Section 4 can be
executed interactively with `docker exec`. The `NET_ADMIN` capability is
required to configure the LAN and WAN network profiles.

See [Quick Start](#41-quick-start) in Section 4
for commands to enter the containers and run experiments.

### 3.3 Standalone Build

Install the project dependencies and compile:

```bash
./shell_install_dependencies.sh
./shell_build_cmd.sh
```

The executable is `./build/fpsi`. Use `./shell_build_cmd.sh --clean` when a
fresh project build is needed.

## 4. Running Experiments

The executable runs both parties in one process using protocol threads and
TCP over `127.0.0.1`. Network configuration and benchmark execution use
separate scripts.

### 4.1 Quick Start

This section provides a first end-to-end run of our fuzzy PSI implementation
and the two baselines used in the paper. It helps you check that the programs
run successfully and obtain initial runtime and communication measurements
for comparison, without running the full experiment matrix. Start with Ours;
then run Exp11 and Exp12 to collect the corresponding baseline results.
The instructions assume that the containers in Section 3 are already running.

The commands below use LAN; replace `lan` with `wan` for a WAN run, and
replace `--preset quick` with `--preset full` for the full matrix.

Estimated running times for the **quick** preset on our evaluation machine
are listed below. Both offline and online
phases are included; actual times depend on hardware and system load.

| Implementation | LAN | WAN |
| --- | ---: | ---: |
| Ours | About 7 minutes | About 17 minutes |
| Exp11: so-OPPRF-based fuzzy PSI | About 7 minutes | About 19 minutes |
| Exp12: da-ROT-based fuzzy PSI | About 15 minutes | About 22 minutes |


#### Ours

Enter the container for our implementation:

```bash
docker exec -it fpsi_cmp_ours bash
```

Inside the container (or from the project root after a standalone build):

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset quick
```

#### Exp11: so-OPPRF-based fuzzy PSI [\[1\]](#ref-1)

Enter the container for the so-OPPRF-based fuzzy PSI implementation:

```bash
docker exec -it fpsi_cmp_exp11 bash
```

Inside the container:

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset quick
```

#### Exp12: da-ROT-based fuzzy PSI [\[2\]](#ref-2)

Enter the container for the da-ROT-based fuzzy PSI implementation:

```bash
docker exec -it fpsi_cmp_exp12 bash
```

Inside the container:

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset quick
```

The comparison scripts write `fpsi_ssoprf_results_<network>_<timestamp>.csv`
and `fpsi_daot_results_<network>_<timestamp>.csv`, respectively, using the same
11 CSV columns as Ours (offline, online, then total measurements).


#### Minimal Run

To check that the executable runs, use one trial of the smallest paper
configuration for each metric:

```bash
./shell_run_bench_fpsi.sh \
  --metric 0 1 2 --nn 8 --dim 2 --delta 60 --trials 1
```

This produces three result rows and is a basic execution check, not a
reproduction of the paper's performance conclusions. An unconfigured network
is labeled `unknown`; select a paper network profile before comparing timings.

### 4.2 Network Configuration

Use `shell_config_network.sh` to reproduce the paper's LAN or WAN conditions
before measuring protocol performance. It limits bandwidth and adds network
delay; it does not run an experiment. After selecting a profile, run
`shell_run_bench_fpsi.sh` separately to collect results under that profile.
Apply the same profile in each implementation's container for comparisons.

| Profile | Bandwidth | Target RTT |
|---|---|---|
| LAN | 10 Gbps | 0 ms added delay |
| WAN | 100 Mbps | 80 ms |
| Custom | User-defined | User-defined |

Run the commands below inside the container being evaluated, or from the
project root for a standalone build. The default interface is `lo` (loopback),
because both parties communicate through `127.0.0.1` on the same machine.
Containers must be started with `NET_ADMIN`; standalone configuration requires
root or sudo.

**Select LAN or WAN.** Choose one profile, then start the benchmark. For LAN:

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset quick
```

To run the same test under WAN instead:

```bash
./shell_config_network.sh wan
./shell_run_bench_fpsi.sh --preset quick
```

Selecting a profile replaces the previous settings on the interface. To
collect both LAN and WAN results, finish the first run before switching
profiles and running again. Do not change the profile during a benchmark.

**Use a custom configuration.** For example, set bandwidth to 1 Gbps and
target round-trip time (RTT) to 20 ms:

```bash
./shell_config_network.sh custom --rate 1Gbps --rtt 20ms
```

`--rate` specifies bandwidth and is required for `custom`; `--rtt` specifies
the target RTT and defaults to `0ms` when omitted. Run the benchmark after
applying the configuration, just as with LAN or WAN.

**Inspect or clear settings.** Check the active profile without changing it:

```bash
./shell_config_network.sh show
```

The benchmark also prints the detected settings at startup. When finished,
remove the traffic-control rules from the selected interface:

```bash
./shell_config_network.sh clear
```

Settings remain active between runs until replaced, cleared, or the container
is removed. Clearing the rules does not select LAN; it removes the bandwidth
and delay constraints. Use `./shell_config_network.sh --help` for all options,
including `--interface IFACE` to select another interface. Running the script
without arguments also displays usage.

### 4.3 Benchmark Script Usage

Use `shell_run_bench_fpsi.sh` to collect runtime and communication measurements
across a set of experiment configurations without launching each case
manually. The script runs the selected parameter combinations, averages the
requested trials for each case, and saves the results to a CSV file. Choose
a predefined quick/full test set or supply your own parameters as described
below. All three projects provide the same script interface; run it from the
corresponding project root after configuring the network in Section 4.2.

#### 4.3.1 Presets: Quick and Full

Quick is the default: `./shell_run_bench_fpsi.sh` is equivalent to
`./shell_run_bench_fpsi.sh --preset quick`. Full averages three trials for
each parameter combination.

| Parameter | Quick | Full |
|---|---|---|
| Metrics | Linf, L1, L2 | Linf, L1, L2 |
| Set size N | 2^12 | 2^8, 2^12, 2^16 |
| Dimension d | 2, 6, 10 | 2, 6, 10 |
| Threshold delta | 60, 250 | 60, 250 |
| Trials per combination | 1 | 3 |
| Combinations: Ours and Exp11 | 18 | 54 |
| Supported combinations: Exp12 | 14 | 42 |

Exp11 is so-OPPRF-based fuzzy PSI [\[1\]](#ref-1). Exp12 is da-ROT-based
fuzzy PSI [\[2\]](#ref-2), which supports L2 only at d=2: its script skips
L2 cases with d=6,10 and produces no CSV rows for them. Linf and L1 use all
listed dimensions in every project.

The presets correspond to Section 8.2 and Appendix D of the camera-ready
paper as follows:

| Paper results | Quick coverage | Full coverage | CSV columns to compare |
|---|---|---|---|
| Table 2: online costs | N=2^12 subset | Complete supported parameter matrix | `Online(s)`, `Online_Com.(MB)` |
| Table 3: overall costs (offline + online) | N=2^12 subset | Select the N=2^12,2^16 rows | `Total(s)`, `Total_Com.(MB)` |

Quick illustrates representative results; full also checks small and large
sets, including cases where a baseline has an advantage. For comparisons,
run Ours, Exp11, and Exp12 with matching metrics, N, d, delta, and network
settings. Evaluate LAN and WAN separately. Section 5 describes the expected
observations and limitations.

After selecting a network profile, choose either preset:

```bash
./shell_run_bench_fpsi.sh --preset quick
# Or run the full matrix:
./shell_run_bench_fpsi.sh --preset full
```

#### 4.3.2 Custom Parameters and Help

Explicit options override preset defaults regardless of argument order.
For example:

```bash
./shell_run_bench_fpsi.sh --preset quick --dim 6 --delta 60 --trials 5
./shell_run_bench_fpsi.sh --preset full --dry-run
./shell_run_bench_fpsi.sh --help
```

`--dry-run` prints the commands without executing the protocols.
Use `--output-dir DIR` to choose where results are saved.

## 5. Results and Paper Claims

### 5.1 Output files and measurements

CSV files are written to the project root unless `--output-dir` is supplied:

```text
fpsi_cmp_results_lan_YYYYMMDD_HHMMSS.csv
fpsi_cmp_results_wan_YYYYMMDD_HHMMSS.csv
fpsi_cmp_results_unknown_YYYYMMDD_HHMMSS.csv
```

The label is derived from the detected network settings: 10 Gbps/0 ms is
`lan`, 100 Mbps/80 ms is `wan`, and other settings are `unknown`.
Network details and trial counts are not included as CSV columns.

| Column | Meaning |
|---|---|
| Protocol | fpsi_cmp |
| Metric | Linf, L1, or L2 |
| Dim, Delta, Size | Dimension, threshold, and set size per party |
| Online_Com.(MB), Online(s) | Average online communication and runtime |
| Offline_Com.(MB), Offline(s) | Average preprocessing communication and runtime |
| Total_Com.(MB), Total(s) | Sum of online and offline measurements |

Runtime is measured at the receiver in seconds. Communication is the sum of
bytes sent and received at the receiver, counting both directions once.
The columns labeled MB use bytes divided by 1024^2 (MiB).
Trials are averaged arithmetically; there is no separate warm-up run.
Total time is the sum of the measured phases, not the entire process wall time.

Offline columns identify the preprocessing overhead. Absolute times can vary
with hardware and system load.

### 5.2 Relation to the paper

The following observations summarize the paper's conclusions for the
experiments defined in [Presets: Quick and Full](#presets-quick-and-full).

| Paper conclusion | Cases to compare | Expected observation |
|---|---|---|
| Online advantages over so-OPPRF-based fuzzy PSI [\[1\]](#ref-1) (Exp11) | All matching Table 2 cases; use N=2^12 for a quick comparison | Ours has lower online communication and lower LAN/WAN runtime across the tested configurations. At N=2^12, reported LAN speedups are 13–27x (Linf), 9.4–18x (L1), and 18–33x (L2). |
| Dimension-dependent advantages over da-ROT-based fuzzy PSI [\[2\]](#ref-2) (Exp12) | Linf/L1 at d=2,6,10, keeping N and delta fixed | Ours has lower online LAN runtime in every tested configuration and lower communication when d > 2. At N=2^12, reported LAN speedups are 11–137x (Linf) and 6.7–81.5x (L1). |
| L2 comparison with da-ROT-based fuzzy PSI [\[2\]](#ref-2) | Quick L2 cases | Reported online improvements are 10–15x for LAN runtime, 1.3–2.1x for WAN runtime, and 3.6–11x for communication. |
| Limits of the online advantage | Include N=2^8, d=2 cases from full, and inspect both network profiles | da-ROT-based fuzzy PSI can use less communication and have lower WAN runtime at small N and d. The paper claims the fastest online LAN performance across tested cases, but the fastest online WAN performance only in most cases. |
| Overall cost, including preprocessing (Appendix D) | Full subset corresponding to Table 3 | Relative to so-OPPRF-based fuzzy PSI, Ours has lower total LAN time at d=2,6, but not always at d=10. Offline communication and latency-sensitive rounds can outweigh online gains in WAN. Total communication can still be lower, especially for Linf/L1 at larger dimensions and thresholds. |

For example, at L1, N=2^16, d=10, delta=250, Ours takes 116.91 s online in
WAN versus 903.56 s for so-OPPRF-based fuzzy PSI, but 1274.03 s in total
versus 921.15 s. This illustrates how preprocessing can reverse an online
runtime advantage.

The dimension experiments cover only d=2,6,10. They do not establish
asymptotic scaling or validate the paper's model-based crossover estimates
against so-OPPRF-based fuzzy PSI. Our asymptotic dependence on d is less
favorable than that baseline's; the observed advantage should not be
extended to arbitrarily high dimensions.

## 6. Executable Usage

Run `./build/fpsi` directly to evaluate a single parameter configuration.
The command-line options below control the set size, dimension, distance
metric, threshold, trial count, and optional CSV output. Display the executable's help message with:

```bash
./build/fpsi --help
```

| Flag | Meaning | Default |
|---|---|---|
| -n N | Direct set size per party | 4096 |
| -nn E | Set size 2^E; overrides -n | Not set |
| -dim D | Point dimension | 6 |
| -delta T | Distance threshold | 60 |
| -metric M | 0=Linf, 1=L1, 2=L2 | 0 |
| -ip ADDR | Socket address | 127.0.0.1 |
| -port P | Base port | 1212 |
| -trait N | Number of trials to average | 5 |
| -out FILE | Append a CSV row | No CSV output |
| -h, -help, --help | Display help | — |

Use the benchmark scripts for presets, network detection, and automatic result
filenames. Direct executable invocations do not detect the network profile.

## 7. License

The original code, scripts, and documentation contributed to `fpsi_cmp` are
licensed under the [MIT License](LICENSE).

Third-party code (including code adapted from other projects) and dependencies
remain subject to their respective licenses and copyright notices. This MIT
license does not relicense the Exp11 or Exp12 baseline implementations.

## References

<a id="ref-1"></a>

\[1\] X. Yang, M. Hao, C. Weng, R. H. Deng, Y. Wen, and T. Zhang.
“Efficient Fuzzy Private Set Intersection from Secret-Shared OPRF.” 2026.

<a id="ref-2"></a>

\[2\] L. Piske, J. Singh, N. Trieu, V. Kolesnikov, and V. Zikas.
“Distance-Aware OT with Application to Fuzzy PSI.”
Cryptology ePrint Archive, Paper 2025/996, 2025.
