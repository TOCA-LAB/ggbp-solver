# ggbp-solver

C++ implementations of exact methods, constructive heuristics, and metaheuristics for the **Generalized Graph Burning Problem (GGBP)**, together with auxiliary optimization problems used to study lower bounds.

The project provides command-line executables for solving individual instances and exporting JSON results. It is research software under active development. Benchmark datasets, experiment recipes, and data-analysis tools are maintained in separate repositories.

## The problem

Given a simple undirected graph and an integer threshold `r >= 1`, a burning process starts with every vertex unburned. At each round:

1. One vertex is selected as an external fire source.
2. Every vertex with at least `r` neighbors burned by the end of the previous round also burns.
3. Burned vertices remain burned.

Propagation is synchronous: vertices that burn in the current round can trigger propagation only in the next round. The objective is to minimize the number of rounds needed to burn the entire graph. This minimum is denoted by `b_r(G)`; `r = 1` gives the classical Graph Burning Problem (GBP).

## Available methods

Executables are generated in `build/bin/`.

| Executable | Method | Role | Uses Gurobi at runtime |
|---|---|---|---|
| `naive` | Time-indexed integer formulation; decision/search and optimization variants; optional AC and RC constraints | Exact GGBP method or LP relaxation | Yes |
| `cover` | Covering formulation with roots consisting of `r` vertices | Exact GGBP method | Yes |
| `unburned_neighbors` | Constructive heuristic based on unburned neighbors | GGBP solution | No |
| `spontaneous_neighbors` | Constructive heuristic based on potential spontaneous propagation | GGBP solution | No |
| `dominating` | Constructive heuristic guided by a dominating set obtained through integer programming | GGBP solution | Yes |
| `grasp` | GRASP with randomized construction and local search | GGBP solution | No |
| `brkga` | Biased Random-Key Genetic Algorithm with a graph-burning decoder | GGBP solution | No |
| `prym` | Implementation of PRYM for the classical GBP | Auxiliary problem, fixed `r = 1` | Yes |
| `pap` | Integer formulation for the Perfect Awareness Problem with constant threshold `r` | Auxiliary problem | Yes |

Heuristic solution lengths are upper bounds on `b_r(G)`, provided that the sequences are feasible. They are not optimality certificates. The optimum of the classical GBP is a lower bound on `b_r(G)`; an optimal PAP seed-set size is another lower bound. A merely feasible PAP seed set does **not** provide that guarantee.

The additional `sandbox` executable is for development and is not part of the solver interface described here.

## Requirements

- A C++23 compiler and standard library.
- CMake 3.16 or newer, as declared by the project.
- Git and a C++ build tool, such as Make or Ninja.
- Boost with its CMake package configuration.
- OpenMP support.
- Gurobi Optimizer, including the C and C++ development libraries. A valid license is required to run the Gurobi-backed methods.

The current CMake configuration requires Gurobi even when building only the heuristic applications. Follow the [official Gurobi installation and license setup instructions](https://support.gurobi.com/hc/en-us/articles/14799677517585-Getting-Started-with-Gurobi-Optimizer).

CMake downloads these dependencies during the first configuration:

| Dependency | Version configured in CMake | Purpose |
|---|---|---|
| [nlohmann/json](https://github.com/nlohmann/json) | 3.12.0 | JSON output |
| [argparse](https://github.com/p-ranav/argparse) | v3.2 | Command-line parsing |
| [Catch2](https://github.com/catchorg/Catch2) | v3.11.0 | Optional tests |

The BRKGA implementation is bundled under `include/brkga-lib/`. Network access is needed for dependency retrieval unless the sources are already available to CMake.

## Build

The following commands assume a Linux installation of Gurobi. Replace the example path with the directory containing its `include/` and `lib/` subdirectories.

```sh
git clone https://github.com/TOCA-LAB/ggbp-solver.git
cd ggbp-solver

export GUROBI_HOME="/path/to/gurobi/linux64"
export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTS=OFF \
  -DGUROBI_DIR="${GUROBI_HOME}"

cmake --build build --parallel 2
```

`GUROBI_HOME` can also be used by the CMake finder directly. `GUROBI_DIR` makes the chosen installation explicit. Configure the license separately according to your Gurobi installation.

To build only selected applications:

```sh
cmake --build build --target unburned_neighbors grasp brkga --parallel 2
```

Runtime contract checks are enabled by default. They can be disabled at configuration time with `-DENABLE_CONTRACTS=OFF`.

## Quick start

Create a small triangle instance and an output directory:

```sh
cat > example.mtx <<'EOF'
%%MatrixMarket matrix coordinate pattern symmetric
% A triangle with three vertices and three edges
3 3 3
1 2
2 3
1 3
EOF

mkdir -p results
```

Run a constructive heuristic with threshold `r = 2`:

```sh
./build/bin/unburned_neighbors example.mtx 2 \
  --json results/unburned.json
```

Inspect the options of any application with `--help`:

```sh
./build/bin/naive --help
./build/bin/grasp --help
./build/bin/brkga --help
```

Most applications use the following invocation pattern:

```text
./build/bin/<executable> <instance-file> <r> [options]
```

`prym` takes only the instance file: it always solves the classical GBP.

### Exact methods and relaxations

The default `naive` mode performs a binary search using a decision formulation. `--optimization` selects the optimization formulation.

```sh
# Decision formulation with binary search
./build/bin/naive example.mtx 2 \
  --timeout 60 --json results/naive-dec.json

# Optimization formulation
./build/bin/naive example.mtx 2 --optimization \
  --timeout 60 --json results/naive-opt.json

# Optimization with both AC and RC constraints
./build/bin/naive example.mtx 2 --optimization \
  --opt-const-active-neighborhood --opt-const-radius-cut \
  --timeout 60 --json results/naive-all.json

# LP relaxation: --optimization is required
./build/bin/naive example.mtx 2 --optimization --linear-relax \
  --timeout 60 --json results/naive-lp.json

# Covering formulation
./build/bin/cover example.mtx 2 \
  --timeout 60 --json results/cover.json
```

AC and RC can also be enabled individually. The same constraint flags can be combined with the LP relaxation.

| Option | Applications | Meaning |
|---|---|---|
| `--timeout SECONDS` | `naive`, `cover`, `prym`, `pap`, `dominating` | Solver time-limit setting |
| `--mem-lim GB` | Same applications | Gurobi memory-limit setting |
| `--soft-mem-lim GB` | Same applications | Gurobi soft memory-limit setting |
| `--lower-bound L` | `naive`, `cover` | Supplied GGBP lower bound |
| `--upper-bound U` | `naive`, `cover` | Supplied GGBP upper bound/horizon |
| `--optimization` (`-opt`) | `naive` | Select the optimization formulation |
| `--linear-relax` (`-lr`) | `naive` | Solve the continuous relaxation; requires `--optimization` |
| `--opt-const-active-neighborhood` (`-ac`) | `naive` | Enable Active Neighborhood constraints |
| `--opt-const-radius-cut` (`-rc`) | `naive` | Enable Radius Cut constraints |

Supplied bounds must be valid for the instance and threshold. Without them, `naive` and `cover` compute initial bounds internally. Time-limit options do not replace an external limit on the entire process, which also includes graph loading and preprocessing.

Cover enumerates roots of size `r`, whose number is `binomial(n, r)`. This can make model construction expensive even before optimization starts.

### Heuristics and metaheuristics

```sh
./build/bin/dominating example.mtx 2 \
  --timeout 60 --json 

./build/bin/spontaneous_neighbors example.mtx 2 \
  --json results/unburned.json

./build/bin/grasp example.mtx 2 --iter 300 --seed 42 \
  --json results/grasp.json

./build/bin/brkga example.mtx 2 \
  --pop 200 --gens 150 --elite 0.15 --mutant 0.10 \
  --rhoe 0.5 --n_pops 2 --MAXT 1 --seed 42 \
  --json results/brkga.json
```

GRASP accepts `--iter` (`-T`) to set its iteration budget. Pass `--seed` explicitly for repeatable experiments; otherwise, GRASP derives its seed from the system clock.

| BRKGA option | Meaning | CLI default |
|---|---|---|
| `--pop` (`-p`) | Size of each population | 100 |
| `--gens` (`-G`) | Number of generations | 100 |
| `--elite` (`-e`) | Elite fraction | 0.20 |
| `--mutant` (`-m`) | Mutant fraction | 0.10 |
| `--rhoe` (`-r`) | Probability of inheriting an allele from the elite parent | 0.70 |
| `--n_pops` (`-K`) | Number of populations | 2 |
| `--MAXT` (`-T`) | Parallel decoding thread setting | 4 |
| `--seed` (`-S`) | Random seed | 0 |

BRKGA currently exchanges its two best individuals between populations every ten generations. Use positive iteration/generation counts and population settings compatible with that exchange. `--MAXT` controls BRKGA decoding; it does not configure Gurobi's thread count.

### Auxiliary problems

```sh
# Classical graph burning: no r argument
./build/bin/prym example.mtx \
  --timeout 60 --json results/prym.json

# Perfect Awareness Problem with constant threshold 2
./build/bin/pap example.mtx 2 \
  --timeout 60 --json results/pap.json
```

PRYM returns a classical burning sequence. PAP returns a seed set for an awareness process, not a GGBP burning sequence. Keep these objective values distinct when comparing results or constructing bounds.

## Input formats

Use `--format mtx` for the default coordinate input, or `--format g6` for graph6.

### Coordinate input (`mtx`)

The current reader accepts the following subset of Matrix Market coordinate files:

```text
% Optional leading comment lines, including a Matrix Market banner
n n m
u_1 v_1
u_2 v_2
...
u_m v_m
```

- Vertex identifiers are **one-based**, from `1` to `n`.
- Each edge record contains exactly two integer endpoints. Use a pattern matrix, without a third weight/value column.
- The graph is interpreted as undirected; duplicate edges and self-loops are ignored.
- Comments belong before the dimension line. Avoid blank lines and comments interspersed with edge records.
- Provide a square dimension line and the correct number of records; the reader is not a general Matrix Market validator.

The small example above follows this format. Benchmark `.mtx` files can also omit the Matrix Market banner and begin directly with the dimensions.

### graph6 (`g6`)

Provide one headerless graph6 record per file, without a `>>graph6<<` prefix or multiple graph records. For example, `Bw` encodes a triangle:

```sh
printf 'Bw\n' > example.g6
./build/bin/unburned_neighbors example.g6 2 --format g6 \
  --json results/unburned-g6.json
```

## JSON results

Use `--json PATH` to write a result file. Create its parent directory beforehand. The applications also print a human-readable summary to standard output.

The JSON schema varies by method; the following fields are useful when processing results:

| Field | Interpretation |
|---|---|
| `instance`, `order`, `size` | Instance filename stem, number of vertices, and number of edges |
| `algorithm` | Method identifier |
| `burning_number` | Length of the returned sequence; the field name alone does not establish optimality |
| `burning_sequence` | Selected vertices, using **zero-based** identifiers |
| `execution_time` | Elapsed program time in seconds, including graph loading and method-specific preprocessing |
| `params` | Metaheuristic settings exported by GRASP or BRKGA |
| `iters` | Best recorded solution length after each GRASP iteration or BRKGA generation |
| `first_best` | Zero-based index of the first recorded step attaining the final best value |
| `completed`, `error_message` | Completion/error information in the Gurobi-backed applications |
| `lower_bound`, `upper_bound` | Bounds exported by the applicable method; auxiliary-problem bounds refer to that auxiliary problem |
| `linear_relaxation` | Fractional bound from a successful `naive` LP relaxation |
| `seed_size`, `seed_set` | PAP objective value and seed vertices, with zero-based identifiers |

Console sequences are printed with **one-based** identifiers, unlike JSON sequences. Some console summaries use the label `b(G)` even for heuristics; interpret the value according to the selected method and its solver status.

Inspect completion, solver status, available solutions, and bounds together. A timeout may leave partial information, and an empty sequence or zero sentinel is not a valid optimum for a nonempty graph. A relaxation value is a lower bound, not an integer burning sequence. Heuristic outputs do not include the same completion/status fields as solver-backed outputs.

For experiment archives, preserve the full command, source revision, input files, seeds, environment, stdout/stderr, and JSON together: JSON alone does not record every setting. The companion experiment repository handles this organization through recipes executed by [xLab](https://github.com/TOCA-LAB/xlab).

## Current limitations
Input support is limited to the formats described above, and result fields are not identical across executables. Preserve the source revision with each experiment so that changes to the implementation or output format remain traceable.

## Tests

Tests are opt-in and use Catch2; the current regression configuration also uses a Python 3 interpreter.

```sh
cmake -S . -B build-tests \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_TESTS=ON \
  -DGUROBI_DIR="${GUROBI_HOME}"

cmake --build build-tests --parallel 2
ctest --test-dir build-tests --output-on-failure
```

Useful subsets are `ctest --test-dir build-tests -L regression-core --output-on-failure` and `ctest --test-dir build-tests -L regression-gurobi --output-on-failure`. Gurobi model tests require a usable installation and license; the suite also contains regression checks for outstanding issues, so a successful build does not imply every regression is resolved.

## Repository layout

```text
app/                    Command-line entry points
include/ggbp/           Graph, simulation, and algorithm interfaces
include/brkga-lib/      Bundled BRKGA implementation and RNG
include/utils/          Shared utilities
src/                    Algorithm and graph implementations
tests/                  Unit and regression tests
FindGUROBI.cmake         Gurobi discovery for CMake
```

## Related research artifacts

- [ggbp-benchmark](https://github.com/TOCA-LAB/ggbp-benchmark): benchmark instances and their source metadata.
- [ggbp-experiments](https://github.com/TOCA-LAB/ggbp-experiments): experiment and analysis recipes.
- [ggbp-analisis](https://github.com/TOCA-LAB/ggbp-analisis): processing of raw results, tables, and figures.

To reproduce a published experiment, use the solver revision and environment specified by its recipe rather than assuming the current branch reproduces historical results.

## Research background and attribution

This repository was developed in the undergraduate research project *Queimando Grafos* at the Federal University of ABC (UFABC), by Lucas Guido under the supervision of Maycon Sambinelli, with support from FAPESP, process **2025/04352-6**.

Relevant references include:

- Yinkui Li, Xiaoxiao Qin, and Wen Li. *The generalized burning number of graphs*. Applied Mathematics and Computation, 411, 126306, 2021. [DOI](https://doi.org/10.1016/j.amc.2021.126306).
- Pereira et al. *A Row Generation Algorithm for Finding Optimal Burning Sequences of Large Graphs*. ESA 2024. [Paper](https://doi.org/10.4230/LIPIcs.ESA.2024.94).
- Rodrigo F. Toso and Mauricio G. C. Resende. *brkgaAPI: A C++ application programming interface for biased random-key genetic algorithms*. [Software and technical report](https://mauricio.resende.info/src/brkgaAPI/).

When using this software in research, identify the repository and the exact commit used, and cite the relevant algorithmic references.

## License

This project is licensed under the **GNU Lesser General Public License v3.0 or later** (LGPL-3.0-or-later). See [`LICENSE`](LICENSE) for the full license texts.

Copyright (C) 2026 Lucas Guido, Maycon Sambinelli / Universidade Federal do ABC

### Third-party components

- **Gurobi** is proprietary software and is **not** included in this repository or covered by this license. To run the ILP components you need your own Gurobi license (free academic licenses are available from Gurobi).
- Any bundled third-party library code remains under its own license.
