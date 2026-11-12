# Width Calculator

This project is licensed under the [MIT License](LICENSE).

## Description

This repository contains a console tool for computing graph parameters related to database join queries and subgraph isomorphism such as:
- treewidth (tw),
- clique embedding power (clemb),
- submodular width (subw),
- #-submodular width (#subw),
- m-width.

This tool can compute these parameters for specified graphs or for all graphs with a given number of nodes or edges.
Besides computing the parameters themselves, one application of the tool is to find examples of graphs for which the computed values among clemb, subw, #subw, and m-width are not all equal; we call such graphs _interesting_.
Treewidth is not included in this comparison.

This tool is supplementary material for the [following paper](TODO).
See the paper for all the definitions of the terms used in this README.


## Installation

1. Clone the repository.
2. Install the Google OR-Tools Python package via `pip` or the Python package manager of your choice.
For example, run `pip install ortools`.
We use the mixed integer linear program solver from this OR-Tools package.
3. Compile a C++ helper file to enhance performance.
Run `g++ -std=c++20 -shared -o helper.so -fPIC helper.cpp -O3`.
You may use a C++ compiler of your choice instead of `g++`.
Minimum C++ version required is C++14.
Please do not change any arguments except for the compiler and the C++ version (i.e., `-shared -o helper.so -fPIC helper.cpp -O3` part should stay the same).
The `helper.cpp` file contains an algorithm for computing all minimal tree decompositions of a given graph.
This is one of the bottlenecks of our algorithm for large graphs, and we use C++ to speed up the computation.
If you do not compile `helper.cpp` into the `helper.so` file, the program will compute minimal tree decompositions in python which may significantly reduce performance for larger graphs.
4. If you decide to use the `--save_pictures_of_interesting_graphs_to` argument, you will also need the `networkx` and `matplotlib` python packages.
You can install them via the `pip install networkx matplotlib` command.
See the next section for the details on the `--save_pictures_of_interesting_graphs_to` argument.


## Usage

The main entry point of the tool is `width-calculator.py`.
Run `python width-calculator.py [options]`, replacing `[options]` with the command line arguments of your choice.
Most of the interaction with the tool is carried via the command line arguments you add to the command.

In the big picture, when calling the tool, you specify all the graph parameters you want to compute and the list of graphs for which you want to compute the parameters.
After that, the tool does all the computation.
It prints the results into the console and a more detailed log into the `calculator.log` file.
(See the `--log_file` argument below to change the log file and the `--nologs` argument to disable the logs.)
You can add multiple sources of input graphs, and all of them will be merged.
If you do not specify any graph input source, you will be asked for a single graph in form of an edge list to be typed into the console.
You have to specify at least one graph parameter to compute or request tree decompositions with `--print_minimal_tree_decompositions` or `--print_interesting_tree_decompositions`.

By default, all graphs with clique separators (i.e., subgraphs that form a clique whose removal disconnects the graph) are skipped because all the graph parameters we compute split over clique separators, i.e., the value for the whole graph is equal to the maximum of the value over all its induced clique-separator-free parts.
(You can disable this behavior by adding the `--allow_clique_separators` argument.)

If for some graph the computed values among clemb, subw, #subw, and m-width are not all equal, such a graph is considered to be _interesting_ and attention is brought to such graphs.
(See the `--ignore_...` family of arguments below to disable this behavior.)

In what follows we describe all the command line arguments in more detail.
You can include any subset of the arguments and specify them in any order.
Please carefully examine the list of all arguments and read about the ones that interest you as (almost) all interaction with the program is carried using the command line arguments.
See the last section of this file for some full examples of how to combine these arguments to compute interesting results.

#### General-Purpose Arguments

- `-h, --help`: Show the help message with the short descriptions of all arguments from this list.
- `-v, --verbosity`: Specify how detailed the console output should be by repeating this command.
Ranges from 1 (least amount of detail) to 5 (most amount of detail).
The default value is 3.
When verbosity is set to 1, only errors and interesting graphs are printed; verbosity 2 adds a progress bar; verbosity 3 adds some basic information and logs about all processed graphs; and verbosity 4 and 5 add more log output.
The easiest way to set this argument is to add one of the following arguments: `-v`, `-vv`, `-vvv`, `-vvvv`, `-vvvvv`.
Here the number of `v`s specifies the verbosity.
Regardless of the verbosity value, you can always find the fully detailed log in the log file (`calculator.log` by default), unless specified otherwise using the `--nologs` argument.
Hence, even if you set your verbosity too low, you can still restore all necessary information.


#### Graph Parameters to Compute

- `--clemb`: Calculate clique embedding power and the corresponding optimal clique embeddings of all input graphs.
(You can see in the logs that sometimes clique embedding power is computed even if this argument was not specified, this is to enhance performance.)
- `--subw`: Calculate submodular width of all input graphs.
- `--sharp_subw`: Calculate #-submodular width of all input graphs.
- `--m_width`: Calculate m-width and the corresponding worst degree exponents of all input graphs.
(You can see in the logs that sometimes m-width is computed even if this argument was not specified, this is to enhance performance.)
- `--treewidth`, `--tree_width`: Calculate treewidth of all input graphs.

#### Graphs to Consider

- `-n, --number_of_nodes`: Specify the number of nodes, and all connected graphs with this number of nodes will be considered.
(Recall that graphs with clique separators are skipped.)
Graphs with any number of nodes between 2 and 8 are available.
The graphs are stored in the `graphs_from_bdm_database` folder.
The graph files included in this repository were generated from [Brendan McKay's graph database](https://users.cecs.anu.edu.au/~bdm/data/graphs.html).
You can specify multiple different node numbers and even ranges.
For example, by specifying `-n 6-8 5-7 2 4`, all graphs with 2, 4, 5, 6, 7, and 8 nodes will be considered.
See also `--force_recalculate_clique_separator_free`.
- `-m, --number_of_edges`: Specify the number of edges, and all connected graphs with this number of edges will be considered.
(Recall that graphs with clique separators are skipped.)
Graphs with any number of edges between 1 and 12 are available.
See `--number_of_nodes` for details.
- `-f, --graph_file`: Specify one or more files from which you want to load graphs.
Every file should contain one or multiple descriptions of graphs.
Each graph description should follow the following format:
In the first line, the number of nodes $n$ and the number of edges $m$ are given; in the next $m$ lines the edges of the graph are given as pairs of numbers.
Node labels may be any distinct integers.
Empty lines anywhere in the file are allowed.
It is also allowed (but not required) to specify the number of graphs in the file in the first line of the file.
If this number is provided, the tool checks that the file contains exactly this many graphs.
See `graphs_from_bdm_database/graphs-with-3-edges.txt` for an example of how to format graphs in a file.
- `-a, --add_graph`: Specify a graph (or graphs) to consider as a string.
Multiple different formats are supported.
The simplest one is an edge list in the form of a python list.
Alternatively, some convenient graph-theoretical names are allowed.
I.e., cliques, multipartite cliques, cycles, parallel path graphs, perfect matching graphs, graph differences (edge removal).
One can combine multiple of the formats to specify multiple graphs.
Each graph description should be enclosed by quotation marks.
For example, by specifying `-a "[(1, 2), (2, 3), (1, 3)]" "K_4" "K_{7}" "K_{3, 3, 3}" "C_5" "P(3, 2, 2)" "K_6\PM_3"`, we add cliques on 3, 4, and 7 nodes, a tripartite clique with 3 nodes in each part, a cycle on 5 nodes, parallel path graph $P(3, 2, 2)$, and a complete graph on 6 nodes from which 3 disjoint edges are removed.
When specifying a graph via a graph-theoretical name, some natural node numbering is chosen.
- `--allow_clique_separators`: By default, all graphs with clique separators are skipped.
By adding this argument, you disable this filtering.

#### Graphs to Ignore

- `--remove_graph`: Specify a graph (or graphs) that should be removed from the list of graphs to consider.
See `--add_graph` for the format.
Graphs are compared up to isomorphism, so node labels do not affect removal.
- `--node_lower_bound`: Specify an integer $n$, and no graph with less than $n$ nodes will be considered.
- `--node_upper_bound`: Specify an integer $n$, and no graph with more than $n$ nodes will be considered.
- `--edge_lower_bound`: Specify an integer $m$, and no graph with less than $m$ edges will be considered.
- `--edge_upper_bound`: Specify an integer $m$, and no graph with more than $m$ edges will be considered.
- `--skip_series_parallel_graphs`: If you add this argument, only graphs with treewidth at least 3 will be considered.
- `--skip_non_series_parallel_graphs`: If you add this argument, only graphs with treewidth at most 2 will be considered.
- `--skip_parallel_path_graphs`: If you add this argument, no parallel path graph will be considered.
- `--skip_non_parallel_path_graphs`: If you add this argument, only parallel path graphs will be considered.
- `--start_from`: If an interruption occurred while running the tool for some list of graphs, this argument can be added to pick up the computation from where it stopped.
You can specify either the graph index to start from or the percentage point.
The index is 1-based; a percentage skips the first $\lfloor pN/100\rfloor$ graphs, where $N$ is the length of the filtered graph list and $p$ is the percentage.
See the last occurrence of the line "current graph index is ..." in the log file for the graph index or "reached ...%" for the last percentage point.
For example, `--start_from 137` or `--start_from 17%`.
- `--skip_graphs_with_clemb_at_least`: Specify some real value $x$.
Clique embedding power will be computed for all graphs, and for graphs with clique embedding power at least $x$, no further width notions will be computed.
- `--skip_graphs_with_clemb_at_most`: Specify some real value $x$.
Clique embedding power will be computed for all graphs, and for graphs with clique embedding power at most $x$, no further width notions will be computed.

#### Interesting Graphs

- `--save_interesting_graphs_to`: Specify a file, and all interesting graphs will be saved to it.
If the file does not exist, it will be created.
Otherwise, it will be overwritten.
- `--save_pictures_of_interesting_graphs_to`: Specify a folder, and the pictures of all interesting graphs will be saved to it.
If some interesting graph is planar, a planar layout will be chosen, otherwise an arbitrary one.
This argument requires the `networkx` and `matplotlib` python packages to be installed.
You can install them via the `pip install networkx matplotlib` command.
The name of each individual file is `x.pdf`, where `x` is the index of the corresponding interesting graph in the graph list.
If the specified folder does not exist, it will be created.
- `--ignore_subw_vs_clemb`: By default, graphs for which the computed values among clemb, subw, #subw, and m-width are not all equal are considered interesting and are highlighted.
By specifying this argument, graphs for which the only inequality is that submodular width is larger than clique embedding power will not be considered interesting.
- `--ignore_sharp_subw_vs_subw`: Analogous to `--ignore_subw_vs_clemb`.
If multiple `--ignore_...` arguments are specified, the definition of an interesting graph is narrowed even further.
- `--ignore_m_width_vs_sharp_subw`: Analogous to `--ignore_subw_vs_clemb`.

#### Tree Decompositions

- `--print_minimal_tree_decompositions`: If this argument is specified, all minimal tree decompositions of each considered graph will be printed into the log file.
- `--print_interesting_tree_decompositions`: If this argument is specified, all interesting tree decompositions of each considered graph will be printed into the log file.
These are the minimal collections of bags left after the easy-bag filtering described below; they need not themselves be tree decompositions.
- `--minimal_tree_decompositions_given`: Specify a file with the tree decompositions of the considered graph.
This argument is allowed only if exactly one graph is provided in the input.
The tree decompositions should be provided as a 3-dimensional python list (list of decompositions, which are lists of bags, which are lists of nodes).
This argument should be used only for larger graphs, for which the computation of all minimal tree decompositions is a costly task which should be done only once and then reused.
No guarantee is given if the provided list of tree decompositions is invalid.
You can use the output of the `--print_minimal_tree_decompositions` argument as the input here.
- `--disable_filtering_FEC_easy_bags`: To speed up the computation of width parameters, by default all bags $B$ with fractional edge-cover number $\mathrm{FEC}(H[B])$ at most the clique embedding power of the whole graph $H$ are removed from the minimal tree decompositions.
Here $H[B]$ is the induced subgraph on $B$; the cover assigns nonnegative weights to its edges and to single nodes so that every node is covered by total weight at least one, and minimizes the sum of the weights.
(See the paper for the proof of why this optimization is valid.)
Specify this argument to disable this behavior.

#### Logs

- `--log_file`: Specify a file to write detailed logs to.
By default, `calculator.log` is used.
Note that the log file is overwritten.
- `--nologs`: If this argument is specified, logs will not be written to the log file.

#### Miscellaneous

- `--print_optimal_submodular_function_for_subw`: If this argument is specified, when submodular and #-submodular widths are computed, all $2^n$ values of the optimal solution of the MILP are printed in lexicographic order (i.e., the i-th value corresponds to the subset of nodes with the indices indicated by the binary representation of i).
- `--fix_degrees`: Fix all or some of the values of the $d$ vector in the MILP for m-width.
This argument is allowed only if exactly one graph is provided in the input.
After this argument, specify exactly $n$ values (one per each node of the graph).
Each value should be either a real number between 0 and 1 or -1 if the current node's degree exponent should not be fixed.
Invalid values and inequality bounds outside this range are rejected, not rounded or clamped.
For the values that are specified, either a decimal or a normal fraction can be used.
Inequalities of the form ">=" or "<=" are allowed as constraints.
For example: `--fix_degrees 0.33333333 0.5 -1 "<=1/2" 1/3 0 3/4 ">=0.15" -1 1`.
(Use quotation marks if you use inequalities.)
By specifying this argument, you restrict the MILP, so only a lower bound on m-width will be found, not the exact value.
Furthermore, when later computing subw and #subw, we also get only lower bounds on these parameters.
Be careful!
- `--check_worst_m_width_degrees_for_ds2_width`: If this argument is specified, after computing the worst degree exponents of a graph from m-width, we check what time complexity exponent they would give as the degree-splitting exponents for ds2-width.
The resulting ds2-width upper bound is printed and recorded in the log.
If the resulting exponent is equal to m-width, this serves as a proof of ds2-width=m-width.
However, if the resulting exponent is not equal to m-width, it does not mean anything as the worst degree exponents for m-width are not necessarily the optimal degree-splitting exponents for ds2-width.
With `--fix_degrees`, the computed m-width value is only a lower bound, while the printed ds2-width value is still an upper bound. Compare this upper bound with m-width from a separate unrestricted run; a valid splitting certificate need not maximize m-width at the same degree vector.
The former misspelling `--check_worst_m_width_degees_for_ds2_width` remains available as a deprecated alias.
- `--force_recalculate_clique_separator_free`: When specifying the `-n` or `-m` arguments, the graphs containing clique separators are filtered out, and the clique-separator-free graphs are cached into the corresponding `...-no-clique-separators.txt` files to improve performance of the subsequent runs.
This process tries to automatically handle mistakes and interruptions but if something goes wrong, you can use this argument to force the tool to recalculate the clique-separator-free graphs from the `graphs_from_bdm_database` folder.
- `--shuffle_graphs`: If this argument is specified, the list of graphs to consider is randomly shuffled before processing.
This argument can be used if processing all specified graphs is infeasible but one wants to at least consider some examples from different parts of the provided list of graphs.
- `--prefer_decimal`: By default, real values are displayed as rational approximations.
Specify this argument to use decimal approximations instead.


## Full examples

`python width-calculator.py --clemb --subw --sharp_subw --m_width -n 3-6 -m 2 3 4 -vv` — find clique embedding power, submodular width, #-submodular width, and m-width for all graphs on 3, 4, 5, or 6 nodes or 2, 3, or 4 edges without clique separators.
Set verbosity to only display the errors, the progress bar, and the interesting graphs.
This finds a single interesting graph on 6 nodes.

`python width-calculator.py --clemb --m_width --subw --sharp_subw -n 7 --start_from 50% -v --save_pictures_of_interesting_graphs_to pictures_7_nodes` — consider the second half of the list of graphs on 7 nodes.
Save all interesting graphs to the `pictures_7_nodes/` folder. Set verbosity to not print anything unnecessary. (Though, without the progress bar, you will have to wait for a long time without any output.)

`python width-calculator.py --clemb --m_width --subw --sharp_subw -n 8 -vv --ignore_subw_vs_clemb` — consider the graphs on 8 nodes.
Only consider a graph interesting if subw < #subw or #subw < m-width.


`python width-calculator.py --m_width --sharp_subw -vvvvv -a "P(9, 2, 2)" --fix_degrees 1/9 1/9 1/9 1/9 2/9 2/9 2/9 2/9 1/9 1/9 5/9 5/9 --minimal_tree_decompositions_given p922-tree-decompositions.txt` — run the computation of m-width and #subw for $P(9, 2, 2)$.
As it cannot be done in a reasonable timespan, we help by providing worst degree exponents that we guessed and the precomputed minimal tree decompositions.
This gives a lower bound of 2-1/9 on m-width and #subw of $P(9, 2, 2)$.

`python width-calculator.py --clemb --m_width -vvvv` — compute clique embedding power and m-width of a graph specified in the console.

`python width-calculator.py --clemb --m_width -f graphs.txt --skip_series_parallel_graphs` — compute clique embedding power and m-width of all graphs specified in the `graphs.txt` file.
Only consider graphs with treewidth at least 3.
