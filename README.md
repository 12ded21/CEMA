# CEMA

CEMA is an experimental solver for the Edge Interdiction Clique Problem (EICP). Given an undirected graph and an edge deletion budget, the program searches for a set of edges to remove so that the maximum clique size in the remaining graph is minimized.

The implementation combines graph reductions, clique-based preprocessing, upper/lower bound initialization, maximum clique subroutines, and an IBM ILOG CPLEX integer programming model with lazy constraints.

## Repository Structure

```text
.
|-- main.cpp              # command-line entry point and experiment pipeline
|-- ilp.cpp / ilp.h       # CPLEX ILP model and lazy constraint callback
|-- cliques.cpp / .h      # clique heuristics, lower bounds, preprocessing
|-- get_mxclq.cpp / .h    # maximum clique subroutine
|-- EdgeListGraph.h       # graph structure, DIMACS parser, reductions
|-- utils.cpp / .h        # utility functions and bound helper functions
|-- run_k_experiments.py  # batch experiments with fixed edge budgets k
|-- run_d_experiments.py  # batch experiments with budget density d
|-- data/                 # benchmark graph instances
|-- results/              # the results of our experiments
`-- CMakeLists.txt        # CMake build configuration
```

Main data directories:

```text
data/ER_1 ... data/ER_5      # generated Erdos-Renyi graphs
data/DIMACS-10               # DIMACS-10 style benchmark instances
data/2nd_dimacs_clique       # DIMACS clique benchmark instances
```

## Requirements

- Linux x86-64
- CMake 3.12 or later
- C++17 compiler, such as `g++`
- IBM ILOG CPLEX Optimization Studio, including CPLEX and Concert
- Python 3, only for the batch experiment scripts

CPLEX is an external dependency and is not stored in this repository. Set its root directory before configuring CMake:

```bash
export CPLEX_ROOT=/path/to/CPLEX
```

Alternatively, pass it directly to CMake with `-DCPLEX_ROOT=/path/to/CPLEX`.

## Build

From the repository root:

```bash
mkdir -p build
cd build
cmake ..
make
```

The executable is generated as:

```text
build/main
```

## Input Format

The solver reads DIMACS edge-list text format:

```text
p edge <number_of_vertices> <number_of_edges>
e <u> <v>
e <u> <v>
...
```

Example:

```text
p edge 4 3
e 1 2
e 2 3
e 3 4
```

Notes:

- The parser skips lines before the first `p` line.
- After the `p` line, graph edges should be provided as `e u v`.
- Self-loops are ignored.
- Duplicate edges are removed.
- The current parser remaps vertices that appear in at least one edge, so isolated vertices in the header are not preserved.

## Single-Instance Usage

```bash
./build/main -l <time_limit> (-k <budget> | -d <density>) [-m <lb_method>] [-e <edge_method>] <input_graph> [output_file]
```

Required arguments:

- `-l <time_limit>`: time limit in seconds.
- `-k <budget>`: integer edge deletion budget.
- `-d <density>`: budget ratio in `[0, 1]`; the program sets `K = ceil(d * |E|)`.
- `<input_graph>`: input graph in DIMACS edge-list format.

`-k` and `-d` are mutually exclusive. Exactly one of them must be provided.

Optional arguments:

- `[output_file]`: if provided, the final tab-separated summary row is appended to this file. Otherwise, the summary row is printed to standard error.
- `-m <lb_method>`: ordering strategy for the initial lower-bound heuristic.
- `-e <edge_method>`: edge selection strategy for clique-based edge banning.

Lower-bound methods for `-m`:

| Value | Method |
| --- | --- |
| `1` | Degree ascending |
| `2` | Degree descending, default |
| `3` | Degeneracy order |
| `4` | Reverse degeneracy order |
| `5` | Random order |

Edge selection methods for `-e`:

| Value | Method |
| --- | --- |
| `1` | Common-neighbor count descending |
| `2` | Common-neighbor count ascending |
| `3` | Degree-sum descending |
| `4` | Degree-sum ascending |
| `5` | Random order, default |

Examples:

```bash
./build/main -l 605 -k 10 -m 2 -e 5 data/2nd_dimacs_clique/c-fat200-1.clq build/example-k.txt
./build/main -l 605 -d 0.001 -m 2 -e 5 data/DIMACS-10/soc-karate build/example-d.txt
```

The program prints detailed progress logs during preprocessing, ILP solving, callbacks, and solution checking.

## Output Columns

Each completed run appends one tab-separated summary row:

| Column | Meaning |
| --- | --- |
| `graph` | input graph file name |
| `initial_n` | number of vertices after reading/remapping |
| `initial_m` | number of undirected edges after reading/removing duplicates |
| `red_n` | vertices after degree reduction |
| `red_m` | edges after degree reduction |
| `red_time` | reduction time in seconds |
| `color_ban_m` | edges remaining after color-based banning |
| `color_ban_time` | color-based banning time in seconds |
| `clique_ban_m` | edges remaining after clique-based banning |
| `clique_ban_time` | clique-based banning time in seconds |
| `all_constraints` | total generated ILP constraints |
| `cut_num` | number of CPLEX user cuts reported by the solver |
| `ILP_time` | ILP solving time in seconds |
| `LB` | initial lower bound |
| `LB_time` | lower-bound computation time in seconds |
| `ans` | final maximum clique size after interdiction |
| `global_time` | total running time in seconds |

Example row:

```text
brock200_1.clq	200	14834	200	14834	0.000185	14834	0.169735	14834	6.99083	1207	468	1.04458	15	0.002741	19	13.5458
```

## Batch Experiments

The Python scripts launch multiple solver runs in parallel.

Fixed-budget experiments:

```bash
cd <project-root>
python3 run_k_experiments.py --time-limit 605 --threads 64 --m_methods 2 --edge_methods 5
```

Density-budget experiments:

```bash
cd <project-root>
python3 run_d_experiments.py --time-limit 605 --threads 16 --m_methods 2 --edge_methods 5
```

Default settings in the scripts:

- `run_k_experiments.py`: `k_values = [10, 15, 20, 25, 30]`
- `run_d_experiments.py`: `d_values = [0.0001, 0.0005, 0.001, 0.0015, 0.002]`
- Both scripts use `data/DIMACS-10` by default.
- Relative `--data-dir` values are resolved from the project root.
- The default executable is `build/main`.

For example, to run on another bundled dataset:

```bash
python3 run_k_experiments.py --data-dir data/2nd_dimacs_clique --time-limit 605 --threads 16
```

To reproduce a different experiment group, edit the budget list, thread count, or output file naming template in the corresponding script.

## Algorithm Pipeline

For each graph instance, `main.cpp` performs the following steps:

1. Read the graph in DIMACS edge-list format.
2. Compute an initial lower bound using a fast greedy clique-coloring heuristic.
3. Apply degree-based graph reduction.
4. Ban edges using color-based and clique-based preprocessing rules.
5. Initialize an upper bound by repeatedly deleting edges from current maximum cliques.
6. Build and solve the ILP model with CPLEX.
7. Add violated clique constraints through a lazy constraint callback.
8. Verify the final solution by recomputing the maximum clique after deleting the selected edges.

Compile-time switches controlling major components are defined near the top of `main.cpp` and `cliques.h`, including `BAN_2`, `DEGREE`, `COLOR`, `EXTRA`, and `CLIQUE`.

## Reproducibility Notes

- Random choices are used in the default edge selection method and upper-bound initialization.
- The current seed is initialized with `srand(time(0))` in `main.cpp`.
- For deterministic experiments, replace this with a fixed seed and record the seed in the result table.

## Citation

If this code is used for a paper, cite the corresponding paper or thesis that introduces the algorithm and experimental setup.
