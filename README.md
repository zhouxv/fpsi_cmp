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

Run the following commands from the repository root.

### Docker

Build the image and open a shell inside it:

```bash
docker build -t fpsi_cmp:latest .
docker run --rm -it --cap-add=NET_ADMIN fpsi_cmp:latest bash
```

Run the experiment commands in Section 4 from that container shell.
To keep CSV files after the container exits, mount a host output directory:

```bash
mkdir -p results
docker run --rm -it --cap-add=NET_ADMIN \
  -v "$PWD/results:/results" fpsi_cmp:latest bash
```

Then pass `--output-dir /results` to the benchmark script.

The project also provides a Docker Hub image:

```bash
docker pull blueobsidian/fpsi_cmp:latest
```

### Local

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

Quick selects 18 combinations and runs one trial per combination.
It uses the moderate set size emphasized in Section 8.2 and retains every
metric and dimension, with the two threshold endpoints.

### Full reproduction

```bash
./shell_config_network.sh lan
./shell_run_bench_fpsi.sh --preset full
```

Full is the default when no preset is specified. It runs all 81 combinations
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
requires matching baseline runs; baseline integration remains pending.

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
