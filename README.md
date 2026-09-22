# Modulo Interval Membership Testing: Towards Practical Fuzzy Private Set Intersection

This artifact implements the fuzzy PSI protocols evaluated in Section 8 and
Table 2 of the paper. A receiver learns which sender points are within a
distance threshold of its own points under the Linf, L1, or L2 metric.

**Badges requested:** Artifact Available, Artifact Functional, and Artifact
Results Reproduced. These are evaluation targets, not awarded badges.
The final evaluated artifact is intended for permanent public archival with a
DOI; the archive link will be added when available.

## 1. Project Structure

```text
.
├── fuzzyPSI/                      Protocol implementation
├── CMakeLists.txt                 Build configuration
├── Dockerfile                     Container build
├── shell_install_dependencies.sh  Dependency installation
├── shell_build_cmd.sh             Project build
├── shell_config_network.sh        Network configuration
└── shell_run_bench_fpsi.sh        Experiments and CSV output
```

## 2. Requirements

### Hardware

An x86-64 CPU is required; no GPU is used. The paper's experiments used an
Intel Xeon Gold 6338 and 128 GB RAM. This describes the paper's machine,
not a measured minimum memory requirement. Peak memory requirements for the
minimal, quick, and full experiments remain to be documented.

### Software

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

### Build Docker Image
Build the image and open a shell inside it:

```bash
docker build -t fpsi_cmp_artifact:latest .
docker run -d --cap-add=NET_ADMIN \
  --name fpsi_cmp_ours \
  fpsi_cmp_artifact:latest \
  sleep infinity
```
### Prebuilt Docker images

Using the prebuilt Docker images is the recommended way to reproduce the
experiments. The source repositories and corresponding Docker images are:

| Implementation | Source repository | Docker image |
|---|---|---|
| Ours | `https://github.com/zhouxv/fpsi_cmp` | `blueobsidian/fpsi_cmp_artifact:latest` |
| [11] so-OPPRF-based fuzzy PSI | `https://github.com/zhouxv/fpsi_ssoprf/tree/fpsi-cmp_artifact_20260919` | `blueobsidian/fpsi_cmp_artifact_exp11:latest` |
| [12] da-ROT-based fuzzy PSI | `https://github.com/zhouxv/fpsi_daOT/tree/fpsi-cmp_artifact_20260916` | `blueobsidian/fpsi_cmp_artifact_exp12:latest` |

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

To enter the container for our implementation:

```bash
docker exec -it fpsi_cmp_ours bash
```

The comparison containers are used in the comparison reproduction subsection
of Section 4.

### Standalone Build

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

### Minimal example

Run one trial of the smallest paper configuration for each metric:

```bash
./shell_run_bench_fpsi.sh \
  --metric 0 1 2 --nn 8 --dim 2 --delta 10 --trials 1
```

This produces three result rows and checks that the program runs. Select the
LAN profile below first if you want to compare the measurements with Table 2.
An unconfigured network is labeled `unknown`.

### Network configuration

```bash
./shell_config_network.sh lan
./shell_config_network.sh wan
./shell_config_network.sh custom --rate 1Gbps --rtt 20ms
```

Choose one configuration before running an experiment:

| Profile | Bandwidth | Target RTT |
|---|---|---|
| LAN | 10 Gbps | 0 ms added delay |
| WAN | 100 Mbps | 80 ms |
| Custom | User-defined | User-defined |

The default interface is `lo`. The benchmark prints the detected settings at
startup. LAN corresponds to the paper's environment; WAN is a supplementary
experiment.

View settings, display help, or clear the rules when finished:

```bash
./shell_config_network.sh show
./shell_config_network.sh --help
./shell_config_network.sh clear
```

Running the network script without arguments displays usage. Network settings
persist until changed, cleared, or the container is removed.

### Quick reproduction

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset quick
```

Quick is the default: `./shell_run_bench_fpsi.sh` is equivalent to the command
above. It selects 18 combinations and runs one trial per combination.
It uses the moderate set size emphasized in Section 8.2 and retains every
metric and dimension, with the two threshold endpoints.

### Full reproduction

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset full
```

Quick is the default when no preset is specified. Full runs all 81 combinations
and averages three trials per combination.

| Parameter | Quick | Full |
|---|---|---|
| Metrics | Linf, L1, L2 | Linf, L1, L2 |
| Set size N | 2^12 | 2^8, 2^12, 2^16 |
| Dimension d | 2, 6, 10 | 2, 6, 10 |
| Threshold delta | 10, 250 | 10, 60, 250 |
| Trials per combination | 1 | 3 |

To collect WAN results, select `wan` and run either preset:

```bash
./shell_config_network.sh wan
./shell_run_bench_fpsi.sh --preset quick
```

### Custom parameters and help

Explicit options override preset defaults regardless of argument order.
For example:

```bash
./shell_run_bench_fpsi.sh --preset quick --dim 6 --delta 60 --trials 5
./shell_run_bench_fpsi.sh --preset full --dry-run
./shell_run_bench_fpsi.sh --help
```

`--dry-run` prints the commands without executing the protocols.
Use `--output-dir DIR` to choose where results are saved.

### Reproducing comparison experiments

The comparison implementations and their Docker images are listed in
Section 3. The following commands reproduce the experiments used for
comparison with [11] and [12]. All three projects use
`./shell_run_bench_fpsi.sh --preset quick|full`, with quick as the default.
The options `--metric`, `--nn`, `--dim`, `--delta`, `--trials`, `--interface`,
`--output-dir`, and `--dry-run` are shared.

The same network profile should be used for our implementation and the
comparison implementation when comparing their results. The commands below
use the LAN profile corresponding to the main experimental setting.

#### Comparison with [11]

Enter the container for the so-OPPRF-based fuzzy PSI implementation:

```bash
docker exec -it fpsi_cmp_exp11 bash
```

Inside the container, configure the LAN network profile:

```bash
cd /workspace
./shell_config_network.sh lan
```

For a quick reproduction, run:

```bash
./shell_run_bench_fpsi.sh --preset quick
```

The quick mode evaluates:

```text
metric = Linf, L1, L2
N      = 2^12
d      = 2, 6, 10
delta  = 10, 250
trials = 1
```

It contains 18 parameter combinations.

For the complete reproduction, run:

```bash
./shell_run_bench_fpsi.sh --preset full
```

The full mode evaluates:

```text
metric = Linf, L1, L2
N      = 2^8, 2^12, 2^16
d      = 2, 6, 10
delta  = 10, 60, 250
trials = 3
```

It contains 81 parameter combinations.

#### Comparison with [12]

Enter the container for the da-ROT-based fuzzy PSI implementation:

```bash
docker exec -it fpsi_cmp_exp12 bash
```

Inside the container, configure the LAN network profile:

```bash
cd /workspace
./shell_config_network.sh lan
```

For a quick reproduction, run:

```bash
./shell_run_bench_fpsi.sh --preset quick
```

For Linf and L1, the quick mode evaluates:

```text
N      = 2^12
d      = 2, 6, 10
delta  = 10, 250
trials = 1
```

For L2, the implementation of [12] supports only `d = 2`. Therefore, the
quick mode contains 14 parameter combinations in total.

For the complete reproduction, run:

```bash
./shell_run_bench_fpsi.sh --preset full
```

For Linf and L1, the full mode evaluates:

```text
N      = 2^8, 2^12, 2^16
d      = 2, 6, 10
delta  = 10, 60, 250
trials = 3
```

For L2, only `d = 2` is evaluated. The full mode therefore contains 63
parameter combinations in total.

The comparison scripts write `fpsi_ssoprf_results_<network>_<timestamp>.csv`
and `fpsi_daot_results_<network>_<timestamp>.csv`, respectively, using the same
11 CSV columns as Ours (offline, online, then total measurements). daOT skips
L2 dimensions other than 2; those combinations have no CSV rows.

When reproducing cross-protocol comparisons, compare results only under the
same metric, set size, dimension, threshold, and network profile.

## 5. Results and Paper Claims

### Output files and measurements

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

Compare the online columns under LAN with the Ours entries in Table 2.
Absolute times can vary with hardware and system load.

### Relation to the paper

| Paper result | Experiment subset | What to inspect |
|---|---|---|
| Section 8.2: advantages over [11] at N=2^12 | All quick cases | Online runtime and communication for the same metric, dimension, and threshold. |
| Section 8.2: high-dimensional advantages over [12] | Quick Linf/L1 at d=2,6,10 with delta fixed | How costs change with dimension, retaining d=2 as a reference. |
| Section 8.2: L2 comparison with [12] | Quick L2 at d=2 | Compare only the supported dimension; [12] does not support higher-dimensional L2. |
| Mild threshold dependence | Quick delta=10 versus 250 at fixed metric and dimension | How much cost changes when the threshold increases 25-fold. |
| Scaling with set size and large-set results | Full | Compare N=2^8,2^12,2^16 at fixed metric, dimension, and threshold. |

For example, at N=2^12, d=10, delta=10, Table 2 reports approximately
178x/102x runtime speedups over [12] for Linf/L1. At the same N and d, raising
delta from 10 to 250 increases Ours communication by approximately
1.20x, 1.28x, and 1.25x for Linf, L1, and L2 respectively. These are published
reference results, not measurements generated by this README.

Quick omits the middle threshold and other set sizes to reduce evaluation
work. It illustrates representative results but does not cover every reported
speedup extremum or establish an asymptotic bound. Full covers the complete
parameter matrix for Ours. Independent cross-protocol speedup validation also
requires matching baseline runs as described in Section 4.

## 6. Executable Usage

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

Use the benchmark script for presets, network detection, and automatic result
filenames. Direct executable invocations do not detect the network profile.
