import os
import sys
import argparse
import re
from math import isfinite
from fractions import Fraction
from time import time
from random import shuffle
from ctypes import c_int, c_char_p, CDLL, CFUNCTYPE, POINTER
from ast import literal_eval
from ortools.linear_solver import pywraplp


# define constants
INF = 1_000_000_000
EPSILON = 1e-6
OUTPUT_ROUNDING_DIGITS = 5
MAX_NODES_FOR_HELPER = 20
# type of the function through which helper.cpp passes its messages: (importance, message)
LOG_CALLBACK = CFUNCTYPE(None, c_int, c_char_p)


# set up console argument parser
def create_parser():
    parser = argparse.ArgumentParser(prog='width-calculator', description='Width Calculator for Graphs. See README for a more detailed description of the arguments.')

    parser.add_argument('--verbosity', '-v', action='count', default=0, help='How verbose the output should be. Value from 1 to 5. Default is 3. Example: -vvvv.')

    parser.add_argument('--clemb', action='store_true', help='Calculate clique embedding power.')
    parser.add_argument('--subw', action='store_true', help='Calculate submodular width.')
    parser.add_argument('--sharp_subw', action='store_true', help='Calculate #-submodular width.')
    parser.add_argument('--m_width', action='store_true', help='Calculate m-width.')
    parser.add_argument('--treewidth', '--tree_width', dest='treewidth', action='store_true', help='Calculate treewidth.')

    parser.add_argument('-n', '--number_of_nodes', nargs='+', help='Process all connected graphs with a given number of nodes.')
    parser.add_argument('-m', '--number_of_edges', nargs='+', help='Process all connected graphs with a given number of edges.')
    parser.add_argument('-f', '--graph_file', nargs='+', help='Input files for reading the graphs.')
    parser.add_argument('-a', '--add_graph', nargs='+', help='Specify input graph(s) as an edge list or as a graph-theoretical name. See README for more details. Example: --add_graph "[(1, 2), (2, 3), (3, 4), (4, 1)]" "K_{3, 3}" "C_5" "P(3, 2, 2)" (remember quotation marks).')
    parser.add_argument('--allow_clique_separators', action='store_true', help='Process graphs with clique separators.')

    parser.add_argument('--remove_graph', nargs='+', help='Specify graph(s) that should be removed from the input list as an edge list or as a graph-theoretical name. See README for more details. Example: --remove_graph "[(1, 2), (2, 3), (3, 4), (4, 1)]" "K_{3, 3}" "C_5" "P(3, 2, 2)" (remember quotation marks).')
    parser.add_argument('--node_lower_bound', type=int, default=0, help='Skip all graphs that have fewer nodes than the parameter.')
    parser.add_argument('--node_upper_bound', type=int, default=INF, help='Skip all graphs that have more nodes than the parameter.')
    parser.add_argument('--edge_lower_bound', type=int, default=0, help='Skip all graphs that have fewer edges than the parameter.')
    parser.add_argument('--edge_upper_bound', type=int, default=INF, help='Skip all graphs that have more edges than the parameter.')
    parser.add_argument('--skip_series_parallel_graphs', action='store_true', help='Skip all graphs with treewidth at most 2.')
    parser.add_argument('--skip_non_series_parallel_graphs', action='store_true', help='Skip all graphs with treewidth at least 3.')
    parser.add_argument('--skip_parallel_path_graphs', action='store_true', help='Skip all parallel path graphs.')
    parser.add_argument('--skip_non_parallel_path_graphs', action='store_true', help='Skip all graphs that are not parallel path graphs')
    parser.add_argument('--start_from', help='Skip all graphs in the filtered graph list that precede the parameter. Use either a 1-based index or a percentage point. Examples: 47 or 16%%.')
    parser.add_argument('--skip_graphs_with_clemb_at_least', help='Do not calculate other width parameters of graphs with clique embedding power at least the parameter. Use either a decimal or a normal fraction.')
    parser.add_argument('--skip_graphs_with_clemb_at_most', help='Do not calculate other width parameters of graphs with clique embedding power at most the parameter. Use either a decimal or a normal fraction.')

    parser.add_argument('--save_interesting_graphs_to', help='Edge lists of all interesting graphs will be stored in the provided file. If the file does not exist, it will be created.')
    parser.add_argument('--save_pictures_of_interesting_graphs_to', help='Pictures (pdf) of all interesting graphs will be stored in the specified folder. If the folder does not exist, it will be created. Requires networkx and matplotlib python packages.')

    parser.add_argument('--ignore_subw_vs_clemb', action='store_true', help='Do not consider interesting the graphs with subw > clemb.')
    parser.add_argument('--ignore_sharp_subw_vs_subw', action='store_true', help='Do not consider interesting the graphs with #subw > subw.')
    parser.add_argument('--ignore_m_width_vs_sharp_subw', action='store_true', help='Do not consider interesting the graphs with m-width > #subw.')

    parser.add_argument('--print_minimal_tree_decompositions', action='store_true', help='Print all minimal tree decompositions of all processed patterns.')
    parser.add_argument('--print_interesting_tree_decompositions', action='store_true', help='Print all interesting tree decompositions of all processed patterns H, i.e., minimal after removing bags B with FEC(H[B]) <= clemb(H).')
    parser.add_argument('--minimal_tree_decompositions_given', help='Provide minimal tree decompositions of a pattern as a 3-dimensional python list in a file.')
    parser.add_argument('--disable_filtering_FEC_easy_bags', action='store_true', help='Do not ignore bags B in tree decompositions that have FEC(H[B]) <= clemb(H).')

    parser.add_argument('--log_file', default='calculator.log', help='Specify the file to write logs to. By default, calculator.log is chosen.')
    parser.add_argument('--nologs', action='store_true', help='Do not write logs to the log file.')

    parser.add_argument('--print_optimal_submodular_function_for_subw', action='store_true', help='Print all values of the optimal solutions for the MILPs for subw and #subw in lexicographic order.')
    parser.add_argument('--fix_degrees', nargs='+', help='Fix all or some values of the d vector in the MILP for m-width. If this parameter is specified, the computed values for m-width, #subw, and subw are only lower bounds on the actual parameters. Only works for a single graph as input. Use either decimal or normal fractions. Inequalities of the form ">=" or "<=" are allowed as constraints. Use -1 for unspecified values. Specify exactly n values. Example: --fix_degrees 0.33333333 0.5 -1 "<=1/2" 1/3 0 3/4 ">=0.15" -1 1.')
    parser.add_argument('--check_worst_m_width_degrees_for_ds2_width', '--check_worst_m_width_degees_for_ds2_width', dest='check_worst_m_width_degrees_for_ds2_width', action='store_true', help='Check whether splitting along the worst degree exponents for m-width leads to a proof of ds2-width = m-width.')
    parser.add_argument('--force_recalculate_clique_separator_free', action='store_true', help='Filter out graphs with clique separators even if they were already calculated before. Use it if something goes wrong.')
    parser.add_argument('--shuffle_graphs', action='store_true', help='Process graphs in a random order.')
    parser.add_argument('--prefer_decimal', action='store_true', help='By default, when printing numbers, a rational approximation is being found. Use this argument to print a decimal approximation instead.')

    return parser


# in the code we store graphs as triplets: (n, edge_list, node_inds),
# where n is the number of nodes, node_inds is a sorted list of node labels,
# and edge_list is a compressed edge list: that is, an edge list where all node
# labels are replaced with values in [0..n - 1].
# this function converts a general edge list into this graph representation.
def generate_graph_from_edge_list(edge_list, n=None):
    node_inds = sorted(list(set(v for e in edge_list for v in e)))
    if n is None:
        n = len(node_inds)
    if len(node_inds) > n:
        raise ValueError(f'{len(node_inds)} different node indices in the edge list of a graph with {n} nodes')
    if len(node_inds) < n:
        raise ValueError(f'{len(node_inds)} different node indices in the edge list of a graph with {n} nodes; isolated nodes are not allowed')
    edge_list = sorted(list(set((min(a, b), max(a, b)) for a, b in edge_list if a != b))) # delete multi-edges and self-loops
    nodes_dict = dict()
    for i, v in enumerate(node_inds):
        nodes_dict[v] = i
    new_edge_list = [(nodes_dict[a], nodes_dict[b]) for a, b in edge_list]
    return n, new_edge_list, node_inds


# given a graph description,
# generates the corresponding multipartite clique
# or returns None if the description is invalid
def try_decode_multipartite_clique_from_description(graph_description):
    pattern = r'^K_\{(\d+(,\d+)*)\}$'
    match = re.match(pattern, graph_description.replace(' ', ''))
    if not match:
        match = re.match(r'^K_(\d+)$', graph_description.replace(' ', ''))
        if not match:
            return None
    part_sizes = list(map(int, match.group(1).split(',')))
    if min(part_sizes) < 1:
        return None
    n = sum(part_sizes)
    if len(part_sizes) == 1:
        cols = [i for i in range(n)]
    else:
        cols = [c for c, part_size in enumerate(part_sizes) for _ in range(part_size)]
    assert len(cols) == n
    edge_list = [(u + 1, v + 1) for u in range(n) for v in range(u) if cols[u] != cols[v]]
    return generate_graph_from_edge_list(edge_list, n)


# given a graph description,
# generates the corresponding cycle
# or returns None if the description is invalid
def try_decode_cycle_from_description(graph_description):
    match = re.match(r'^C_\{(\d+)\}$', graph_description.replace(' ', ''))
    if not match:
        match = re.match(r'^C_(\d+)$', graph_description.replace(' ', ''))
        if not match:
            return None
    l = int(match.group(1))
    if l < 3:
        return None
    edge_list = [(v + 1, (v + 1) % l + 1) for v in range(l)]
    return generate_graph_from_edge_list(edge_list, l)


# given a graph description,
# generates the corresponding parallel path graph
# or returns None if the description is invalid
def try_decode_parallel_path_graph_from_description(graph_description):
    pattern = r'^P\((\d+(,\d+)*)\)$'
    match = re.match(pattern, graph_description.replace(' ', ''))
    if not match:
        return None
    lens = list(map(int, match.group(1).split(',')))
    if min(lens) < 1 or lens.count(1) > 1:
        return None
    n = sum(lens) - len(lens) + 2
    edge_list = []
    curv = 3
    for l in lens:
        if l == 1:
            edge_list.append((1, 2))
        else:
            edge_list.append((1, curv))
            for i in range(l - 2):
                edge_list.append((curv, curv + 1))
                curv += 1
            edge_list.append((curv, 2))
            curv += 1
    assert curv == n + 1
    return generate_graph_from_edge_list(edge_list, n)


# given a graph description,
# generates the corresponding perfect matching
# or returns None if the description is invalid
def try_decode_perfect_matching_from_description(graph_description):
    match = re.match(r'^PM_\{(\d+)\}$', graph_description.replace(' ', ''))
    if not match:
        match = re.match(r'^PM_(\d+)$', graph_description.replace(' ', ''))
        if not match:
            return None
    k = int(match.group(1))
    n = 2 * k
    edge_list = [(v + 1, v + k + 1) for v in range(k)]
    return generate_graph_from_edge_list(edge_list, n)


# given graphs G and H, compute graph (V(G), E(G) \ E(H))
def graph_difference(g1, g2):
    if g1[0] != g2[0] or g1[2] != g2[2]:
        return None
    return g1[0], [e for e in g1[1] if not (e in g2[1])], g1[2]


# given a string representing a graph, compute the graph representation;
# returns None if the string is not in any supported format and
# raises ValueError if it is in a supported format but does not describe a valid graph
def try_decode_graph_from_description(graph_description):
    if graph_description.count('\\') == 1:
        g1, g2 = map(try_decode_graph_from_description, graph_description.split('\\'))
        if g1 and g2:
            return graph_difference(g1, g2)
    for decoding_method in [try_decode_multipartite_clique_from_description, try_decode_cycle_from_description, try_decode_parallel_path_graph_from_description, try_decode_perfect_matching_from_description]:
        graph = decoding_method(graph_description)
        if graph:
            return graph
    try:
        edge_list = literal_eval(graph_description)
    except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
        return None
    try:
        return generate_graph_from_edge_list(edge_list)
    except (ValueError, TypeError) as error:
        raise ValueError('expected an edge list such as "[(1, 2), (2, 3)]"') from error


# decode a graph given in the command line or report an error
def decode_graph_argument(parser, option, graph_description):
    try:
        graph = try_decode_graph_from_description(graph_description)
    except ValueError as error:
        parser.error(f'{option}: invalid graph description {graph_description!r}: {error}')
    if graph is None:
        parser.error(f'{option}: invalid graph description {graph_description!r}; see README for the supported formats')
    return graph


# given a graph representation, returns its edge list with correct node labels
def get_uncompressed_edge_list(graph):
    n, edge_list, node_inds = graph
    return [(node_inds[a], node_inds[b]) for a, b in edge_list]


# given a graph representation, returns its adjacency list
def convert_graph_to_adjacency_list(graph):
    n, edge_list, _ = graph
    adj_list = [[] for _ in range(n)]
    for a, b in edge_list:
        adj_list[a].append(b)
        adj_list[b].append(a)
    return adj_list


# compare graph structures without depending on the node labels
def are_graphs_isomorphic(graph, other_graph):
    n, edge_list, _ = graph
    if n != other_graph[0] or len(edge_list) != len(other_graph[1]):
        return False
    if edge_list == other_graph[1]:
        return True
    adj_list = convert_graph_to_adjacency_list(graph)
    other_adj_list = convert_graph_to_adjacency_list(other_graph)
    degrees = [len(neighbours) for neighbours in adj_list]
    other_degrees = [len(neighbours) for neighbours in other_adj_list]
    if sorted(degrees) != sorted(other_degrees):
        return False
    adj_masks = [sum(1 << u for u in neighbours) for neighbours in adj_list]
    other_adj_masks = [sum(1 << u for u in neighbours) for neighbours in other_adj_list]
    candidates = [[u for u in range(n) if degrees[v] == other_degrees[u]] for v in range(n)]
    mapping = [-1 for _ in range(n)]

    def search(used_mask):
        best_v = None
        best_matches = None
        for v in range(n):
            if mapping[v] != -1:
                continue
            matches = [u for u in candidates[v] if not (used_mask & (1 << u)) and all(
                bool(adj_masks[v] & (1 << w)) == bool(other_adj_masks[u] & (1 << mapping[w]))
                for w in range(n) if mapping[w] != -1)]
            if not matches:
                return False
            if best_matches is None or len(matches) < len(best_matches):
                best_v, best_matches = v, matches
        if best_v is None:
            return True
        for u in best_matches:
            mapping[best_v] = u
            if search(used_mask | (1 << u)):
                return True
        mapping[best_v] = -1
        return False

    return search(0)


# test whether a given mask has a specified bit turned on
def test_bit(mask, bit):
    return ((mask >> bit) & 1) == 1


# converts a subset of [0 ... n - 1] into a bitmask between 0 and 2^n-1
def list_to_mask(lst):
    return sum([1 << v for v in lst])


# converts a bitmask between 0 and 2^n-1 into a subset of [0 ... n - 1]
def mask_to_list(mask):
    lst = []
    k = 0
    while (1 << k) <= mask:
        if test_bit(mask, k):
            lst.append(k)
        k += 1
    return lst


# check whether a given graph does not contain a clique separator
def is_graph_clique_separator_free(graph):
    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    adj_list_mask = [(1 << v) + sum([(1 << u) for u in adj_list[v]]) for v in range(n)]
    for clique_mask in range((1 << n) - 1):
        is_clique = True
        for v in mask_to_list(clique_mask):
            if (adj_list_mask[v] & clique_mask) != clique_mask:
                is_clique = False
                break
        if not is_clique:
            continue

        # check whether the graph without this clique is connected or not
        used = [test_bit(clique_mask, v) for v in range(n)]
        bfs_start_node = min([v for v in range(n) if not used[v]])
        que = [bfs_start_node]
        used[bfs_start_node] = True
        while len(que) > 0:
            v = que[0]
            que = que[1:]
            for u in adj_list[v]:
                if not used[u]:
                    used[u] = True
                    que.append(u)
        if not all(used):
            return False
    return True


# reads one graph as an edge list, ignoring empty lines
def read_graph(args, file=sys.stdin, interactive=False, first_line=None):
    try:
        if interactive:
            print('Specify the number of nodes and the number of edges in the graph in one line')
            print('Example: 5 7')
        line = first_line if first_line is not None else file.readline()
        while line and not line.strip():
            line = file.readline()
        if not line:
            return None
        n, m = map(int, line.split())
        if n < 0 or m < 0:
            raise ValueError('the numbers of nodes and edges must be nonnegative')
        edge_list = []
        if interactive:
            print('List edges by giving two integer node indices that should be connected per line')
        for i in range(m):
            line = file.readline()
            while line and line.strip() == '':
                line = file.readline()
            if not line:
                raise ValueError(f'unexpected end of input while reading edge {i + 1} of {m}')
            a, b = map(int, line.split())
            edge_list.append((a, b))
        return generate_graph_from_edge_list(edge_list, n)
    except Exception as error:
        raise ValueError(f'Invalid graph description: {error}') from error


# read all graphs and check the optional graph-count header
def read_graphs(args, file, require_count=False):
    line = file.readline()
    while line and not line.strip():
        line = file.readline()
    if not line:
        if require_count:
            raise ValueError('missing graph-count header')
        return
    declared_count = None
    if len(line.split()) == 1:
        try:
            declared_count = int(line.strip())
        except ValueError as error:
            raise ValueError('invalid graph-count header') from error
        if declared_count < 0:
            raise ValueError('the graph count must be nonnegative')
        line = None
    elif require_count:
        raise ValueError('missing graph-count header')
    count = 0
    while True:
        graph = read_graph(args, file, first_line=line)
        line = None
        if graph is None:
            break
        count += 1
        if declared_count is not None and count > declared_count:
            raise ValueError(f'expected {declared_count} graphs, found more')
        yield graph
    if declared_count is not None and count != declared_count:
        raise ValueError(f'expected {declared_count} graphs, found {count}')


# given a graph representation, print its edge list either into a file or to the console
def print_graph(graph, file=sys.stdout):
    n, edge_list, _ = graph
    file.write(f'{n} {len(edge_list)}\n')
    for a, b in get_uncompressed_edge_list(graph):
        file.write(f'{a} {b}\n')
    file.write('\n')


# check whether a given file has a number at the top, and then
# the corresponding number of graphs as edge lists. in other words,
# verify that no interruption happened when this file was created.
def is_valid_filtered_graph_file(path):
    try:
        if not os.path.isfile(path):
            return False
        with open(path, 'r') as file:
            for _ in read_graphs(None, file, require_count=True):
                pass
        return True
    except Exception as _:
        return False


# read graphs from {src}.txt and write to {src}-no-clique-separators.txt
# the ones that do not contain a clique separator
def filter_out_graphs_with_clique_separators(args, src):
    filtered_path = src + '-no-clique-separators.txt'
    if args.force_recalculate_clique_separator_free or not is_valid_filtered_graph_file(filtered_path):
        unfiltered_path = src + '.txt'
        cnt_graphs = 0
        with open(unfiltered_path, 'r') as input_file, open(filtered_path, 'w') as output_file:
            for graph in read_graphs(args, input_file):
                if is_graph_clique_separator_free(graph):
                    print_graph(graph, output_file)
                    cnt_graphs += 1
        with open(filtered_path, 'r+') as file:
            content = file.read()
            file.seek(0)
            file.write(str(cnt_graphs) + '\n\n' + content)
    return filtered_path


# check whether a given graph is a parallel path graph with no clique separator
def is_parallel_path_graph(graph):
    if not is_graph_clique_separator_free(graph):
        return False
    n, edge_list, _ = graph
    degs = [0 for v in range(n)]
    for a, b in edge_list:
        degs[a] += 1
        degs[b] += 1
    return len([v for v in range(n) if degs[v] != 2]) <= 2


# check whether a given graph has treewidth at most two
# in other words, search for a K_4 minor
def is_treewidth_at_most_two(graph):
    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    # consider 4^n colorings but w.l.o.g. assume that node n-1 is in the 0th part
    for mask in range(1 << (2 * n - 2)):
        k4_adj_matrix = [[i == j for j in range(4)] for i in range(4)]
        cols = [((mask >> (2 * i)) & 3) for i in range(n)]
        for v in range(n):
            for u in adj_list[v]:
                k4_adj_matrix[cols[v]][cols[u]] = True
        is_4_clique = all([all(lst) for lst in k4_adj_matrix])
        if not is_4_clique:
            continue

        # now check that each color forms a connected subset of nodes
        used = [False for _ in range(n)]
        for col in range(4):
            bfs_start = min([v for v in range(n) if cols[v] == col])
            que = [bfs_start]
            used[bfs_start] = True
            while len(que) > 0:
                v = que[0]
                que = que[1:]
                for u in adj_list[v]:
                    if cols[u] == col and not used[u]:
                        used[u] = True
                        que.append(u)
        if all(used):
            return False
    return True


# check whether a given graph is connected
def is_connected(graph):
    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    used = [False for _ in range(n)]
    que = [0]
    used[0] = True
    while len(que) > 0:
        v = que[0]
        que = que[1:]
        for u in adj_list[v]:
            if not used[u]:
                used[u] = True
                que.append(u)
    return all(used)


# check whether a given graph satisfies the conditions
# required to be in the list of graphs to check
def graph_satisfies_simple_conditions(graph, args):
    if graph[0] <= 1 or not is_connected(graph):
        return False
    if not args.allow_clique_separators and not is_graph_clique_separator_free(graph):
        return False
    if args.skip_series_parallel_graphs and is_treewidth_at_most_two(graph):
        return False
    if args.skip_non_series_parallel_graphs and not is_treewidth_at_most_two(graph):
        return False
    if args.skip_parallel_path_graphs and is_parallel_path_graph(graph):
        return False
    if args.skip_non_parallel_path_graphs and not is_parallel_path_graph(graph):
        return False
    return args.edge_lower_bound <= len(graph[1]) and len(graph[1]) <= args.edge_upper_bound and args.node_lower_bound <= graph[0] and graph[0] <= args.node_upper_bound


# print a message into a log file if it exists;
# furthermore, if the message is important enough,
# print it into the console
def print_log(args, logs_file, importance, *vargs):
    if logs_file:
        print(*vargs, file=logs_file, flush=True)
    if importance <= args.verbosity:
        print(*vargs, flush=True)


# set up types of arguments and return types of functions from helper.cpp;
# return None (and fall back to python) if helper.so is missing or cannot be loaded
def setup_ctypes(args, log_file):
    if not os.path.isfile('helper.so'):
        return None
    try:
        libc = CDLL(os.path.abspath('helper.so'))
        libc.enumerate_minimal_tree_decompositions.argtypes = [c_int, c_int, POINTER(c_int), LOG_CALLBACK]
        libc.enumerate_minimal_tree_decompositions.restype = POINTER(c_int)
        libc.free_memory.argtypes = [POINTER(c_int)]
    except (OSError, AttributeError) as error:
        print_log(args, log_file, 1, f'WARNING: could not load helper.so ({error}); computing tree decompositions in python instead, which may be slow. Recompile helper.cpp as described in README to fix this.')
        return None
    return libc


# put a pdf with a picture of a given graph into a given directory
def draw_graph(graph, graph_index, short_graph_description, directory, nx, plt):
    n, _, node_inds = graph
    nxgraph = nx.Graph()
    for v in node_inds:
        nxgraph.add_node(v)
    for a, b in get_uncompressed_edge_list(graph):
        nxgraph.add_edge(a, b)

    pos = nx.planar_layout(nxgraph) if nx.check_planarity(nxgraph)[0] else nx.spring_layout(nxgraph)
    nx.draw(nxgraph, pos, with_labels=True, node_color='skyblue', node_size=500, font_size=11, font_color='white', font_weight='bold')
    plt.text(0.02, 1.04, short_graph_description, fontsize=11, ha='left', va='top', transform=plt.gca().transAxes)
    plt.savefig(os.path.join(directory, f'{graph_index}.pdf'), bbox_inches='tight')
    plt.clf()


# given a floating value, find the closest fraction
# with denominator at most max_denom
def float_to_frac_approximation(val, max_denom=1000):
    if val < 0:
        numerator, denominator = float_to_frac_approximation(-val, max_denom)
        return -numerator, denominator

    opt_numerator = 0
    opt_denominator = 1
    for denominator in range(1, max_denom + 1):
        approx_numerator = int(val * denominator)
        for numerator in range(max(0, approx_numerator - 5), approx_numerator + 5):
            if numerator * opt_denominator != denominator * opt_numerator and abs(numerator / denominator - val) < abs(opt_numerator / opt_denominator - val):
                opt_numerator = numerator
                opt_denominator = denominator
    return opt_numerator, opt_denominator


# convert a floating number into a string to print;
# depending on the settings, returns
# either a decimal number or a fraction
def float_to_string(val, args):
    if args.prefer_decimal:
        return str(round(val, OUTPUT_ROUNDING_DIGITS))
    else:
        a, b = float_to_frac_approximation(val)
        if b == 1:
            return str(a)
        else:
            return f'{a}/{b}'


# convert a string into a floating number;
# works both with decimal numbers and fractions as input
def string_to_float(s):
    s = s.replace(' ', '')
    if '/' in s:
        numerator, denominator = map(float, s.split('/'))
        return numerator / denominator
    else:
        return float(s)


# parse a fixed degree, an inequality, or the -1 sentinel
def parse_degree_constraint(raw_str):
    value_str = ''.join(raw_str.split())
    relation = '='
    if value_str.startswith(('>=', '<=')):
        relation, value_str = value_str[:2], value_str[2:]
    try:
        value = string_to_float(value_str)
    except (ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError(f'invalid degree constraint {raw_str!r}') from error
    if relation == '=' and value == -1:
        return None, None
    if not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'degree constraint {raw_str!r} must use a value between 0 and 1; use -1 for an unspecified degree')
    return relation, value


# a percentage skips the corresponding initial fraction of the graph list
def start_index(start_from, graph_count):
    if not start_from:
        return 0
    if start_from.endswith('%'):
        try:
            percent = Fraction(start_from[:-1])
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError('the starting percentage must be between 0% and 100%') from error
        if not 0 <= percent <= 100:
            raise ValueError('the starting percentage must be between 0% and 100%')
        return int(graph_count * percent / 100)
    try:
        index = int(start_from)
    except ValueError as error:
        raise ValueError('the starting graph index must be a positive integer') from error
    if index < 1:
        raise ValueError('the starting graph index must be a positive integer')
    return index - 1


# converts a list of strings that represent either positive integers or ranges of positive integers (e.g. "3-5") into a list of unique integers
def list_of_ranges_to_list_of_values(lst):
    ans = []
    for s in lst:
        s = s.replace(' ', '')
        try:
            if '-' in s:
                l, r = map(int, s.split('-'))
            else:
                l = r = int(s)
        except ValueError as error:
            raise ValueError(f'invalid number or range {s!r}; use, e.g., 5 or 3-7') from error
        if l < 0 or l > r:
            raise ValueError(f'invalid range {s!r}; the lower end must be nonnegative and at most the upper end')
        ans.extend(range(l, r + 1))
    return sorted(list(set(ans)))


# print multiple representations of a floating number into the logs
def print_graph_parameter(param_name, param_value, args, log_file, importance=3):
    numerator, denominator = float_to_frac_approximation(param_value)
    print_log(args, log_file, importance, f'{param_name} ≈ {round(param_value, OUTPUT_ROUNDING_DIGITS)} ≈', str(numerator) if denominator == 1 else f'{numerator}/{denominator} = {numerator // denominator} + {numerator % denominator}/{denominator}')


# solve a mixed integer linear program to proven optimality
def solve_to_optimality(solver):
    parameters = pywraplp.MPSolverParameters()
    parameters.SetDoubleParam(pywraplp.MPSolverParameters.RELATIVE_MIP_GAP, 0.0)
    return solver.Solve(parameters)


# given a graph H and a subset B of its nodes,
# compute fractional edge covering of H[B]
def fractional_edge_covering(bag, graph, args, log_file):
    if len(bag) == 0:
        return 0
    n, edge_list, node_inds = graph
    filtered_edge_list = [(a, b) for (a, b) in edge_list if a in bag and b in bag] + [(v, v) for v in bag]

    solver = pywraplp.Solver.CreateSolver('SCIP')
    edge_value = [solver.NumVar(0.0, 1.0, f'edge_value_{i}') for i in range(len(filtered_edge_list))]
    for v in bag:
        affecting_vars = []
        for var, edge in zip(edge_value, filtered_edge_list):
            if v in edge:
                affecting_vars.append(var)
        solver.Add(sum(affecting_vars) >= 1.0)
    solver.Minimize(sum(edge_value))

    status = solve_to_optimality(solver)
    if status != pywraplp.Solver.OPTIMAL:
        raise Exception(f'Failed while computing fractional edge covering of the graph {get_uncompressed_edge_list((n, filtered_edge_list, node_inds))}')
    return solver.Objective().Value()


# weight of a maximum-weight spanning tree of the bags,
# where two bags are joined by an edge whose weight is the size of their intersection
def max_spanning_tree_weight(bag_masks):
    k = len(bag_masks)
    in_tree = [False for _ in range(k)]
    best = [-1 for _ in range(k)]
    best[0] = 0
    total = 0
    for _ in range(k):
        i = max((j for j in range(k) if not in_tree[j]), key=lambda j: best[j])
        in_tree[i] = True
        total += best[i]
        for j in range(k):
            if not in_tree[j]:
                best[j] = max(best[j], bin(bag_masks[i] & bag_masks[j]).count('1'))
    return total


# read tree decompositions of a given graph from a file given as a python list
# and check that each of them is a valid tree decomposition of the graph;
# the line printed by --print_minimal_tree_decompositions is also accepted
def read_tree_decompositions(input_path, graph):
    with open(input_path, 'r') as input_file:
        content = input_file.read().strip()
    prefix = 'minimal tree decompositions ='
    if content.startswith(prefix):
        content = content[len(prefix):]
    try:
        tree_decompositions = literal_eval(content.strip())
    except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError) as error:
        raise ValueError('expected a 3-dimensional python list (decompositions, bags, nodes)') from error
    if not isinstance(tree_decompositions, (list, tuple)) or len(tree_decompositions) == 0:
        raise ValueError('expected a nonempty list of tree decompositions')
    n, edge_list, _ = graph
    result = []
    for i, tree_decomposition in enumerate(tree_decompositions, 1):
        if not isinstance(tree_decomposition, (list, tuple)) or len(tree_decomposition) == 0 or not all(isinstance(bag, (list, tuple)) for bag in tree_decomposition):
            raise ValueError(f'tree decomposition {i} is not a nonempty list of bags')
        for bag in tree_decomposition:
            for v in bag:
                if not isinstance(v, int) or isinstance(v, bool) or not 0 <= v < n:
                    raise ValueError(f'tree decomposition {i} contains node {v!r}; nodes must be 0-based indices 0, ..., {n - 1} of the sorted node labels, not the labels themselves')
        bags = [sorted(set(bag)) for bag in tree_decomposition]
        bag_masks = [list_to_mask(bag) for bag in bags]
        for a, b in edge_list:
            if not any(test_bit(mask, a) and test_bit(mask, b) for mask in bag_masks):
                raise ValueError(f'tree decomposition {i} has no bag containing the edge ({a}, {b})')
        node_counts = [sum(1 for mask in bag_masks if test_bit(mask, v)) for v in range(n)]
        if 0 in node_counts:
            raise ValueError(f'tree decomposition {i} has no bag containing node {node_counts.index(0)}')
        # the bags can be arranged into a tree in which the bags containing each node form a subtree
        # if and only if a maximum-weight spanning tree has weight sum_v (number of bags containing v - 1)
        if max_spanning_tree_weight(bag_masks) != sum(node_counts) - n:
            raise ValueError(f'the bags of tree decomposition {i} cannot be arranged into a tree decomposition')
        result.append(bags)
    return result


# sort a list and remove repetitions
def make_unique(lst):
    ans = []
    for elem in sorted(lst):
        if len(ans) == 0 or ans[-1] != elem:
            ans.append(elem)
    return ans


# compute all minimal tree decompositions of a graph
# call a C++ function from helper.cpp for better performance
def compute_minimal_tree_decompositions(graph, args, log_file, libc):
    print_log(args, log_file, 3, 'computation of minimal tree decompositions started')
    n, edge_list, _ = graph
    if libc is not None and n > MAX_NODES_FOR_HELPER:
        print_log(args, log_file, 1, f'WARNING: helper.so supports at most {MAX_NODES_FOR_HELPER} nodes; computing tree decompositions of this graph with {n} nodes in python instead')
    if libc is not None and n <= MAX_NODES_FOR_HELPER:
        c_int_edge_list = (c_int * (2 * len(edge_list)))(*[v for edge in edge_list for v in edge])
        log_callback = LOG_CALLBACK(lambda importance, message: print_log(args, log_file, importance, message.decode()))
        res_ptr = libc.enumerate_minimal_tree_decompositions(n, len(edge_list), c_int_edge_list, log_callback)
        res_len = res_ptr[0]
        res = [res_ptr[i] for i in range(res_len)]
        libc.free_memory(res_ptr)

        # decode
        tree_decompositions = [[] for _ in range(res[1])]
        ind = 2
        for i in range(len(tree_decompositions)):
            tree_decompositions[i] = [[] for _ in range(res[ind])]
            ind += 1
            for j in range(len(tree_decompositions[i])):
                tree_decompositions[i][j] = mask_to_list(res[ind])
                ind += 1
        assert ind == res_len

        print_log(args, log_file, 3, 'computation of minimal tree decompositions finished')
        return tree_decompositions

    # if helper.so does not exist, compute tree decompositions directly in python
    adj_list = convert_graph_to_adjacency_list(graph)
    tree_decompositions = []
    permutation = [i for i in range(n)]
    while True:
        adj_mtx = [list_to_mask(adj_list[v]) for v in range(n)]
        tree_decomposition = []
        # generate a tree decomposition
        alive_mask = (1 << n) - 1
        for v in permutation:
            tree_decomposition.append(((adj_mtx[v] | (1 << v)) & alive_mask))
            alive_mask ^= (1 << v)
            for u in range(n):
                if test_bit(adj_mtx[v], u):
                    adj_mtx[u] = ((adj_mtx[u] | adj_mtx[v]) ^ (1 << u))
        # remove bag repetitions and bags that are subsets of other bags
        tree_decomposition = make_unique(tree_decomposition)
        filtered_tree_decomposition = []
        for bag in tree_decomposition:
            superset_exists = False
            for bag2 in tree_decomposition:
                if bag != bag2 and (bag & bag2) == bag:
                    superset_exists = True
                    break
            if not superset_exists:
                filtered_tree_decomposition.append(bag)
        tree_decompositions.append(filtered_tree_decomposition)

        # generate next permutation
        i = n - 1
        while i > 0 and permutation[i] < permutation[i - 1]:
            i -= 1
        if i == 0:
            break
        j = i - 1
        while i + 1 < n and permutation[i + 1] > permutation[j]:
            i += 1
        permutation[i], permutation[j] = permutation[j], permutation[i]
        for k in range(j + 1, n):
            l = n + j - k
            if k < l:
                permutation[k], permutation[l] = permutation[l], permutation[k]

    tree_decompositions = make_unique(tree_decompositions)

    # leave only minimal tree decompositions
    minimal_tree_decompositions = []
    for tree_dec in tree_decompositions:
        is_superdec = False
        for other_tree_dec in tree_decompositions:
            if tree_dec == other_tree_dec:
                continue
            is_superdec_cur = True
            for bag in other_tree_dec:
                is_subbag = False
                for bag2 in tree_dec:
                    if (bag & bag2) == bag:
                        is_subbag = True
                        break
                if not is_subbag:
                    is_superdec_cur = False
                    break
            if is_superdec_cur:
                is_superdec = True
                break
        if not is_superdec:
            minimal_tree_decompositions.append([mask_to_list(bag) for bag in tree_dec])
    print_log(args, log_file, 3, 'computation of minimal tree decompositions finished')
    return minimal_tree_decompositions


# delete all easy bags from the tree decompositions in the list
# and compute the minimal tree decompositions in the resulting list
def compute_interesting_tree_decompositions_from_minimal_tree_decompositions(tree_decompositions, easy_bags):
    filtered_tree_decompositions = make_unique([sorted([list_to_mask(bag) for bag in tree_dec if not easy_bags[list_to_mask(bag)]]) for tree_dec in tree_decompositions])
    minimal_filtered_tree_decompositions = []
    for tree_dec in filtered_tree_decompositions:
        is_superdec = False
        for other_tree_dec in filtered_tree_decompositions:
            if other_tree_dec == tree_dec:
                continue
            is_superdec_cur = True
            for bag in other_tree_dec:
                is_subbag = False
                for bag2 in tree_dec:
                    if (bag & bag2) == bag:
                        is_subbag = True
                        break
                if not is_subbag:
                    is_superdec_cur = False
                    break
            if is_superdec_cur:
                is_superdec = True
                break
        if not is_superdec:
            minimal_filtered_tree_decompositions.append(tree_dec)
    return [[mask_to_list(bag) for bag in tree_dec] for tree_dec in sorted(minimal_filtered_tree_decompositions)]


# find an optimal clique embedding into a given graph
def find_clemb(graph, args, log_file):
    print_log(args, log_file, 3, 'computation of clique embedding power started')

    n, edge_list, node_inds = graph
    adj_list = convert_graph_to_adjacency_list(graph)

    # which subsets of nodes induce connected subgraphs?
    induces_connected_subgraph = [False for _ in range(1 << n)]
    for mask in range(1, 1 << n):
        bfs_start_node = min(mask_to_list(mask))
        que = [bfs_start_node]
        used = [not test_bit(mask, v) for v in range(n)]
        used[bfs_start_node] = True
        while len(que) > 0:
            v = que[0]
            que = que[1:]
            for u in adj_list[v]:
                if not used[u]:
                    used[u] = True
                    que.append(u)
        induces_connected_subgraph[mask] = all(used)

    solver = pywraplp.Solver.CreateSolver('SCIP')
    # weights embedded into all subgraphs
    weight = [solver.NumVar(0.0, 1.0, f'weight_{mask}') for mask in range(1 << n)]
    # is weight strictly positive?
    used = [solver.IntVar(0, 1, f'used_{mask}') for mask in range(1 << n)]
    OBJ = solver.NumVar(0.0, 1.0, 'OBJ') # 1/clemb (smallest weak edge depth if weights sum up to 1)

    solver.Add(sum(weight) == 1.0)
    for mask in range(1 << n):
        solver.Add(weight[mask] <= used[mask]) # if used=0, nothing can be embedded
        if not induces_connected_subgraph[mask]:
            solver.Add(weight[mask] == 0.0)
            continue
        for mask2 in range(mask + 1, 1 << n):
            if not induces_connected_subgraph[mask2]:
                continue
            subsets_touch = ((mask & mask2) != 0)
            for v in mask_to_list(mask):
                subsets_touch |= any([test_bit(mask2, u) for u in adj_list[v]])
            if not subsets_touch:
                # all pairs of subsets with something embedded should touch
                solver.Add(used[mask] + used[mask2] <= 1)
    for u, v in edge_list:
        # OBJ >= weak depth of (u, v)
        solver.Add(OBJ >= sum([weight[mask] for mask in range(1 << n) if test_bit(mask, v) or test_bit(mask, u)]))
    solver.Minimize(OBJ)

    status = solve_to_optimality(solver)
    print_log(args, log_file, 3, 'computation of clique embedding power finished')
    if status != pywraplp.Solver.OPTIMAL:
        raise Exception(f'Failed while computing clique embedding power of the graph {get_uncompressed_edge_list(graph)}')

    if args.clemb:
        def embedding_string(mask):
            nodes = [node_inds[v] for v in mask_to_list(mask)]
            return f'(nodes={nodes}, weight={float_to_string(weight[mask].solution_value(), args)})'
        selected = {embedding_string(mask) for mask in range(1 << n) if weight[mask].solution_value() > EPSILON}
        selected_string = str(selected).replace("'", "")
        print_log(args, log_file, 3, f'number of selected subgraphs for an optimal clique embedding: {len(selected)} \nselected subgraphs for an optimal clique embedding: {selected_string}')
    return 1 / solver.Objective().Value()


# compute phi_d_hat from the definition of m-width for a fixed vector d
def compute_phi_d_hat_for_fixed_degrees(graph, tree_decompositions, degrees):
    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    MAX_VAL = n + 1
    phi_d = [MAX_VAL for _ in range(1 << n)]
    phi_d[0] = 0
    for mask in range(1 << n):
        for v in mask_to_list(mask):
            # introduce node
            phi_d[mask] = min(phi_d[mask], phi_d[mask ^ (1 << v)] + (1.0 - degrees[v]))
            # search for v as a neighbour of u
            for u in adj_list[v]:
                if test_bit(mask, u):
                    phi_d[mask] = min(phi_d[mask], phi_d[mask ^ (1 << v)] + degrees[u])
    phi_d_hat = [MAX_VAL for _ in range(1 << n)]
    for mask in range(1 << n):
        for mask2 in range(1 << n):
            if (mask & mask2) == mask:
                # monotonicity
                phi_d_hat[mask] = min(phi_d_hat[mask], phi_d[mask2])
    return phi_d_hat


# compute max_i f_d^i for the fixed splitting degree exponents d
def compute_ds2_width_upper_bound_for_fixed_splitting_degrees(graph, tree_decompositions, degrees):
    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    MAX_VAL = n + 1
    OBJ = 0.0
    for hl_mask in range(1 << n):
        try_phi = [MAX_VAL for _ in range(1 << n)]
        try_phi[0] = 0.0
        for bag_mask in range(1, 1 << n):
            for v in mask_to_list(bag_mask):
                #introduce node
                try_phi[bag_mask] = min(try_phi[bag_mask], try_phi[bag_mask ^ (1 << v)] + (1.0 - (degrees[v] if test_bit(hl_mask, v) else 0.0)))
                for u in adj_list[v]:
                    if test_bit(bag_mask, u):
                        # introduce edge
                        try_phi[bag_mask] = min(try_phi[bag_mask], try_phi[bag_mask ^ (1 << v) ^ (1 << u)] + 1.0)
                        if not test_bit(hl_mask, u):
                            # search for v as a neighbour of u
                            try_phi[bag_mask] = min(try_phi[bag_mask], try_phi[bag_mask ^ (1 << v)] + degrees[u])
        for bag_mask in range(1 << n):
            for bag_mask2 in range(1 << n):
                if (bag_mask & bag_mask2) == bag_mask:
                    # enforce monotonicity
                    try_phi[bag_mask] = min(try_phi[bag_mask], try_phi[bag_mask2])
        best_tree_dec_value = min([max([try_phi[list_to_mask(bag)] for bag in tree_dec]) for tree_dec in tree_decompositions])
        OBJ = max(OBJ, best_tree_dec_value)
    return OBJ


myind = 0
# sets var to be the maximum of the variables from varlist
def add_max_equality(solver, MAX_VAL, var, varlist):
    global myind
    myind += 1
    assert len(varlist) > 0
    for var2 in varlist:
        solver.Add(var >= var2)
    intvars = [solver.IntVar(0, 1, f'add_max_equality_bool_var_{myind}_{i}') for i in range(len(varlist))]
    solver.Add(sum(intvars) == 1)
    for i in range(len(varlist)):
        solver.Add(var <= varlist[i] + MAX_VAL * (1 - intvars[i]))


# find m-width of a given graph and the corresponding worst degree exponents
def find_m_width(graph, tree_decompositions, args, log_file, lower_bound_on_m_width):
    print_log(args, log_file, 3, 'computation of m-width started')
    if tree_decompositions == [[]]:
        print_log(args, log_file, 3, 'there exists a fractional-edge-covering-easy tree decomposition')
        print_log(args, log_file, 3, 'computation of m-width finished')
        return lower_bound_on_m_width, None

    n = graph[0]
    adj_list = convert_graph_to_adjacency_list(graph)
    MAX_VAL = n + 1

    solver = pywraplp.Solver.CreateSolver('SCIP')
    degexp = [solver.NumVar(0, 1, f'degexp_{v}') for v in range(n)]
    phi = [solver.NumVar(0, MAX_VAL, f'phi_{mask}') for mask in range(1 << n)]
    # maximum of phi_hat over all bags of each tree decomposition
    tree_dec_value = [solver.NumVar(0, MAX_VAL, f'tree_dec_value_{i}') for i in range(len(tree_decompositions))]
    # m-width
    OBJ = solver.NumVar(0 if args.fix_degrees else lower_bound_on_m_width, MAX_VAL, 'OBJ')

    simple_success = False
    worst_degrees = None
    m_width = None
    # if --fix_degrees argument is set, fix degree exponents of specified nodes
    if args.fix_degrees:
        cnt_fixed_degrees = 0
        degree_constraints = [parse_degree_constraint(raw_str) for raw_str in args.fix_degrees]
        for v, (relation, value) in enumerate(degree_constraints):
            if relation == '>=':
                solver.Add(degexp[v] >= value)
            elif relation == '<=':
                solver.Add(degexp[v] <= value)
            elif relation == '=':
                solver.Add(degexp[v] == value)
                cnt_fixed_degrees += 1
        if cnt_fixed_degrees == n: # if all degrees are fixed, phi_d_hat is defined uniquely
            worst_degrees = [value for _, value in degree_constraints]
            phi_d_hat_vals = compute_phi_d_hat_for_fixed_degrees(graph, tree_decompositions, worst_degrees)
            m_width = max(lower_bound_on_m_width, min([max([phi_d_hat_vals[list_to_mask(bag)] for bag in decomposition]) for decomposition in tree_decompositions]))
            simple_success = True

    if not simple_success:
        for bag_mask in range(1 << n):
            if bag_mask == 0:
                solver.Add(phi[bag_mask] == 0)
            else:
                for v in mask_to_list(bag_mask):
                    # introduce node
                    solver.Add(phi[bag_mask] <= phi[bag_mask ^ (1 << v)] + (1.0 - degexp[v]))
                    # search for v as a neighbour of u
                    for u in adj_list[v]:
                        if test_bit(bag_mask, u):
                            solver.Add(phi[bag_mask] <= phi[bag_mask ^ (1 << v)] + degexp[u])
            for v in range(n):
                if not test_bit(bag_mask, v):
                    # enforce monotonicity
                    solver.Add(phi[bag_mask] <= phi[bag_mask | (1 << v)])
        for i, tree_dec in enumerate(tree_decompositions):
            # set tree_dec_value[i] to be equal to the maximum of phi_hat over all bags of tree_dec
            add_max_equality(solver, MAX_VAL, tree_dec_value[i], [phi[list_to_mask(bag)] for bag in tree_dec])
            # OBJ is the minimum over all tree decompositions
            solver.Add(OBJ <= tree_dec_value[i])
        solver.Maximize(OBJ)
        status = solve_to_optimality(solver)
        if status == pywraplp.Solver.INFEASIBLE: # after deleting easy bags, some tree decomposition became too easy
            worst_degrees = [0 for _ in range(n)]
            m_width = lower_bound_on_m_width
        elif status != pywraplp.Solver.OPTIMAL:
            raise Exception(f'Failed while computing m-width of the graph {get_uncompressed_edge_list(graph)}')
        else:
            worst_degrees = [di.solution_value() for di in degexp]
            m_width = max(solver.Objective().Value(), lower_bound_on_m_width)

    print_log(args, log_file, 3, 'computation of m-width finished')
    if args.m_width:
        print_log(args, log_file, 4, 'worst degree exponents:', ', '.join([float_to_string(d, args) for d in worst_degrees]))
    if args.check_worst_m_width_degrees_for_ds2_width:
        ds2_bound = compute_ds2_width_upper_bound_for_fixed_splitting_degrees(graph, tree_decompositions, worst_degrees)
        print_graph_parameter('ds2-width upper bound', ds2_bound, args, log_file)
        if abs(ds2_bound - m_width) < EPSILON:
            print_log(args, log_file, 3, 'worst degree exponents give a proof of m-width=ds2-width')
        elif args.fix_degrees:
            print_log(args, log_file, 3, 'fixed-degree lower bound and splitting upper bound do not certify equality')
        else:
            print_log(args, log_file, 2, f'worst degree exponents DO NOT give a proof of m-width=ds2-width for graph {get_uncompressed_edge_list(graph)}')
    return m_width, worst_degrees


# find submodular or #submodular width of a given graph
def find_subw(graph, tree_decompositions, is_sharp_subw, worst_m_degrees, args, log_file, lower_bound_on_subw, upper_bound_on_subw):
    width_name = '#submodular' if is_sharp_subw else 'submodular'
    print_log(args, log_file, 3, f'computation of {width_name} width started')
    if tree_decompositions == [[]] or (not args.print_optimal_submodular_function_for_subw and upper_bound_on_subw - lower_bound_on_subw < EPSILON):
        print_log(args, log_file, 3, f'computation of {width_name} width finished')
        return lower_bound_on_subw
    n, edge_list, _ = graph
    adj_list = convert_graph_to_adjacency_list(graph)
    MAX_VAL = n + 1
    simple_success = False
    opt_submodular_function = None
    ans = None
    target_width = max(lower_bound_on_subw, upper_bound_on_subw)
    # which intersections are allowed for submodular constraints
    good_intersection = [True for _ in range(1 << n)]
    if is_sharp_subw:
        good_intersection = [False for _ in range(1 << n)]
        good_intersection[0] = True
        for v in range(n):
            good_intersection[1 << v] = True
        for u, v in edge_list:
            good_intersection[(1 << v) | (1 << u)] = True

    if worst_m_degrees:
        # check whether submodular function from m-width computation also works for this case
        phi_d_hat = compute_phi_d_hat_for_fixed_degrees(graph, tree_decompositions, worst_m_degrees)
        is_submodular = True
        for mask in range(1 << n):
            for mask2 in range(mask + 1, 1 << n):
                if good_intersection[mask & mask2]:
                    is_submodular &= (phi_d_hat[mask & mask2] + phi_d_hat[mask | mask2] <= phi_d_hat[mask] + phi_d_hat[mask2] + EPSILON)
        # A previously computed bound can exceed this candidate's width,
        # particularly when the degree exponents were fixed.
        candidate_width = min(max(phi_d_hat[list_to_mask(bag)] for bag in decomposition) for decomposition in tree_decompositions)
        if is_submodular and abs(candidate_width - target_width) < EPSILON:
            print_log(args, log_file, 3, f'computation of {width_name} width finished')
            simple_success = True
            opt_submodular_function = [val for val in phi_d_hat]
            ans = candidate_width

    if not simple_success:
        # first attempt fixes phi({v}) = 1 - worst_m_degrees[v]
        # and tries to prove that the current width is equal to m-width;
        # if it doesn't succeed, the second attempt runs a general solver
        for attempt in range(2):
            if attempt == 0 and not worst_m_degrees:
                continue
            solver = pywraplp.Solver.CreateSolver('SCIP')
            phi = [solver.NumVar(0, MAX_VAL, f'phi_{mask}') for mask in range(1 << n)]
            # maximum of phi over all bags of each tree decomposition
            tree_dec_value = [solver.NumVar(0, MAX_VAL, f'tree_dec_value_{i}') for i in range(len(tree_decompositions))]
            # subw
            OBJ = solver.NumVar(0 if args.fix_degrees else lower_bound_on_subw, MAX_VAL, 'OBJ')

            # add degree constraints
            if attempt == 0:
                for v in range(n):
                    solver.Add(phi[1 << v] == 1.0 - worst_m_degrees[v])
            solver.Add(phi[0] == 0)
            # enforce edge-domination
            for u, v in edge_list:
                solver.Add(phi[(1 << u) | (1 << v)] <= 1.0)
            # enforce submodular constraints
            for mask in range(1 << n):
                if is_sharp_subw:
                    for mask2 in range(mask + 1, 1 << n):
                        if good_intersection[mask & mask2]:
                            solver.Add(phi[mask & mask2] + phi[mask | mask2] <= phi[mask] + phi[mask2])
                for v in mask_to_list(mask):
                    solver.Add(phi[mask] >= phi[mask ^ (1 << v)]) # enforce monotonicity constraint
                    if not is_sharp_subw:
                        for u in mask_to_list(mask):
                            if u < v:
                                solver.Add(phi[mask ^ (1 << v) ^ (1 << u)] + phi[mask] <= phi[mask ^ (1 << v)] + phi[mask ^ (1 << u)])
            for i, tree_dec in enumerate(tree_decompositions):
                # set tree_dec_value[i] to be equal to the maximum of phi over all bags of tree_dec
                add_max_equality(solver, MAX_VAL, tree_dec_value[i], [phi[list_to_mask(bag)] for bag in tree_dec])
                # OBJ is the minimum over all tree decompositions
                solver.Add(OBJ <= tree_dec_value[i])
            solver.Maximize(OBJ)

            status = solve_to_optimality(solver)
            if attempt == 0 and (status != pywraplp.Solver.OPTIMAL or solver.Objective().Value() < target_width - EPSILON):
                continue
            print_log(args, log_file, 3, f'computation of {width_name} width finished')
            if status == pywraplp.Solver.INFEASIBLE:
                return lower_bound_on_subw # after deleting easy bags, some tree decomposition became too easy
            elif status != pywraplp.Solver.OPTIMAL:
                raise Exception(f'Failed while computing {width_name} width of the graph {get_uncompressed_edge_list(graph)}')
            opt_submodular_function = [var.solution_value() for var in phi]
            ans = solver.Objective().Value()
            break

    if args.print_optimal_submodular_function_for_subw:
        print_log(args, log_file, 5, f'optimal {width_name} function values:', str([float_to_string(val, args) for val in opt_submodular_function]).replace("'", ""))
    return max(ans, lower_bound_on_subw)


# main entry point
def main():
    start_time = time()

    parser = create_parser()
    args = parser.parse_args()
    if args.verbosity == 0:
        args.verbosity = 3
    if not 1 <= args.verbosity <= 5:
        parser.error('verbosity must be between 1 and 5')
    if args.fix_degrees:
        try:
            for raw_str in args.fix_degrees:
                parse_degree_constraint(raw_str)
        except ValueError as error:
            parser.error(f'--fix_degrees: {error}')
    # parse the clemb thresholds
    for option in ['skip_graphs_with_clemb_at_least', 'skip_graphs_with_clemb_at_most']:
        raw_str = getattr(args, option)
        if raw_str is None:
            continue
        try:
            value = string_to_float(raw_str)
        except (ValueError, ZeroDivisionError, OverflowError):
            value = None
        if value is None or not isfinite(value):
            parser.error(f'--{option}: expected a decimal or a fraction such as 1.5 or 17/9, got {raw_str!r}')
        setattr(args, option, value)

    if not args.clemb and not args.subw and not args.sharp_subw and not args.m_width and not args.treewidth and not args.print_interesting_tree_decompositions and not args.print_minimal_tree_decompositions:
        print('\n\nSpecify at least one graph parameter to compute or request tree decompositions\n\n')
        parser.print_help()
        print('\n\nSpecify at least one graph parameter to compute or request tree decompositions\n\n')
        return

    # collect all input files
    input_paths = []
    if args.graph_file:
        for file in args.graph_file:
            input_paths.append(file)
    # graphs from bdm database are once filtered from graphs with clique separators,
    # and then the result is saved and reused in later calls
    try:
        if args.number_of_nodes:
            for n in list_of_ranges_to_list_of_values(args.number_of_nodes):
                src = os.path.join('graphs_from_bdm_database', f'graphs-with-{n}-nodes')
                input_paths.append(src + '.txt' if args.allow_clique_separators else filter_out_graphs_with_clique_separators(args, src))
        if args.number_of_edges:
            for m in list_of_ranges_to_list_of_values(args.number_of_edges):
                src = os.path.join('graphs_from_bdm_database', f'graphs-with-{m}-edges')
                input_paths.append(src + '.txt' if args.allow_clique_separators else filter_out_graphs_with_clique_separators(args, src))
    except (OSError, ValueError) as error:
        parser.error(str(error))

    log_file = None if args.nologs else open(args.log_file, 'w')
    # collect all graphs to consider
    print_log(args, log_file, 3, f'Collecting the graphs')
    graphs = []
    for input_path in input_paths:
        try:
            with open(input_path, 'r') as input_file:
                graphs.extend(read_graphs(args, input_file))
        except (OSError, ValueError) as error:
            parser.error(f'{input_path}: {error}')
    if args.add_graph:
        for graph_description in args.add_graph:
            graphs.append(decode_graph_argument(parser, '--add_graph', graph_description))
    if len(input_paths) == 0 and not args.add_graph:
        cur_time = time() - start_time
        try:
            graph = read_graph(args, interactive=True)
        except ValueError as error:
            parser.error(str(error))
        start_time = time() - cur_time # do not include console input time into the running time
        if graph is None:
            print_log(args, log_file, 2, 'No graph provided. Nothing to process.')
            return
        graphs.append(graph)
    graphs_to_remove = []
    if args.remove_graph:
        for graph_description in args.remove_graph:
            graphs_to_remove.append(decode_graph_argument(parser, '--remove_graph', graph_description))
    unsupported_graph_count = sum(1 for graph in graphs if graph[0] <= 1 or not is_connected(graph))
    if unsupported_graph_count > 0:
        print_log(args, log_file, 1, f'WARNING: skipping {unsupported_graph_count} input graph(s) that are disconnected or have fewer than two nodes; only connected graphs with at least two nodes are processed')
    graphs = [graph for graph in graphs if graph_satisfies_simple_conditions(graph, args) and not any(are_graphs_isomorphic(graph, other_graph) for other_graph in graphs_to_remove)]
    if args.fix_degrees and (len(graphs) != 1 or len(args.fix_degrees) != graphs[0][0]):
        parser.error('--fix_degrees requires exactly one input graph and exactly one value per node')
    given_tree_decompositions = None
    if args.minimal_tree_decompositions_given:
        if len(graphs) != 1:
            parser.error('--minimal_tree_decompositions_given requires exactly one input graph')
        try:
            given_tree_decompositions = read_tree_decompositions(args.minimal_tree_decompositions_given, graphs[0])
        except (OSError, ValueError) as error:
            parser.error(f'--minimal_tree_decompositions_given: {error}')
    try:
        index_to_start = start_index(args.start_from, len(graphs))
    except ValueError as error:
        parser.error(f'--start_from: {error}')
    if args.shuffle_graphs:
        shuffle(graphs)

    if args.save_pictures_of_interesting_graphs_to:
        try:
            import networkx as nx
        except ImportError:
            raise Exception('networkx package is required for the --save_pictures_of_interesting_graphs_to argument but is not installed')
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise Exception('matplotlib package is required for the --save_pictures_of_interesting_graphs_to argument but is not installed')
        os.makedirs(args.save_pictures_of_interesting_graphs_to, exist_ok=True)
    interesting_graphs_file = open(args.save_interesting_graphs_to, 'w') if args.save_interesting_graphs_to else None
    interesting_graph_descriptions = []
    libc = setup_ctypes(args, log_file)

    print_log(args, log_file, 2, f'Collected {len(graphs)} graphs. Starting to process them.')
    processed_graph_count = 0
    for cur_index, graph in enumerate(graphs):
        if cur_index < index_to_start:
            continue
        # Report completed positions, rounding down so a resume cannot skip
        # the graph that is about to be processed.
        cur_percent = 100 * cur_index // len(graphs)
        if cur_index == index_to_start or 100 * (cur_index - 1) // len(graphs) < cur_percent:
            if args.verbosity > 2:
                print_log(args, log_file, 1, '\n')
            print_log(args, log_file, 2, f'reached {cur_percent}%')
        processed_graph_count += 1
        n, m = graph[0], len(graph[1])
        print_log(args, log_file, 3, f'\n\ncurrent graph index is {cur_index + 1}')
        print_log(args, log_file, 4, f'nodes = {n}, edges = {m}, edge list = {get_uncompressed_edge_list(graph)}')

        forbid_FEC_filtering = args.disable_filtering_FEC_easy_bags or args.print_optimal_submodular_function_for_subw or args.minimal_tree_decompositions_given or args.check_worst_m_width_degrees_for_ds2_width or args.fix_degrees
        use_FEC_filtering = (args.subw or args.sharp_subw or args.m_width or args.print_interesting_tree_decompositions) and not forbid_FEC_filtering
        lower_bound_on_all_parameters = -1.0
        clemb = None
        if args.clemb or use_FEC_filtering or args.skip_graphs_with_clemb_at_least is not None or args.skip_graphs_with_clemb_at_most is not None:
            tmp_clemb = find_clemb(graph, args, log_file)
            lower_bound_on_all_parameters = max(lower_bound_on_all_parameters, tmp_clemb)
            if args.skip_graphs_with_clemb_at_least is not None and tmp_clemb >= args.skip_graphs_with_clemb_at_least - EPSILON:
                print_log(args, log_file, 3, 'skipping this graph as its clemb is at least the specified threshold')
                continue
            if args.skip_graphs_with_clemb_at_most is not None and tmp_clemb <= args.skip_graphs_with_clemb_at_most + EPSILON:
                print_log(args, log_file, 3, 'skipping this graph as its clemb is at most the specified threshold')
                continue
            if args.clemb:
                clemb = tmp_clemb
                print_graph_parameter('clemb', clemb, args, log_file)

        # bags with FEC-value at most clemb
        easy_bags = [fractional_edge_covering(mask_to_list(mask), graph, args, log_file) <= lower_bound_on_all_parameters + EPSILON for mask in range(1 << n)] if use_FEC_filtering else [False for _ in range(1 << n)]
        interesting_tree_decompositions = None
        treewidth = None
        if args.subw or args.sharp_subw or args.m_width or args.treewidth or args.print_interesting_tree_decompositions or args.print_minimal_tree_decompositions:
            minimal_tree_decompositions = given_tree_decompositions if given_tree_decompositions is not None else compute_minimal_tree_decompositions(graph, args, log_file, libc)
            if args.print_minimal_tree_decompositions:
                print_log(args, log_file, 5, f'minimal tree decompositions = {minimal_tree_decompositions}')
            if args.treewidth:
                treewidth = min([max(map(len, decomposition)) for decomposition in minimal_tree_decompositions]) - 1
                print_graph_parameter('treewidth', treewidth, args, log_file)
            interesting_tree_decompositions = compute_interesting_tree_decompositions_from_minimal_tree_decompositions(minimal_tree_decompositions, easy_bags)
            if args.print_interesting_tree_decompositions:
                print_log(args, log_file, 5, f'interesting tree decompositions = {interesting_tree_decompositions}')

        m_width = None
        worst_m_degrees = None
        upper_bound_on_sharp_subw = None
        # overall, we compute m-width as its value helps to compute other parameters
        if args.m_width or args.subw or args.sharp_subw:
            tmp_m_width, worst_m_degrees = find_m_width(graph, interesting_tree_decompositions, args, log_file, lower_bound_on_all_parameters)
            upper_bound_on_sharp_subw = tmp_m_width
            if args.m_width:
                m_width = tmp_m_width
                print_graph_parameter('m-width', m_width, args, log_file)

        subw = None
        lower_bound_on_sharp_subw = lower_bound_on_all_parameters
        if args.subw:
            subw = find_subw(graph, interesting_tree_decompositions, False, worst_m_degrees, args, log_file, lower_bound_on_all_parameters, upper_bound_on_sharp_subw)
            lower_bound_on_sharp_subw = subw
            print_graph_parameter('subw', subw, args, log_file)

        sharp_subw = None
        if args.sharp_subw:
            sharp_subw = find_subw(graph, interesting_tree_decompositions, True, worst_m_degrees, args, log_file, lower_bound_on_sharp_subw, upper_bound_on_sharp_subw)
            print_graph_parameter('#subw', sharp_subw, args, log_file)


        params_string = ', '.join([f'{name} = {float_to_string(value, args)}' for name, value in [('clemb', clemb), ('subw', subw), ('#subw', sharp_subw), ('m-width', m_width), ('treewidth', treewidth)] if value])
        full_graph_description = f'nodes = {n}, edges = {m}, edge list = {get_uncompressed_edge_list(graph)}, graph index = {cur_index + 1}, ' + params_string
        short_graph_description = f'n = {n}, m = {m}, idx = {cur_index + 1}, ' + params_string
        print_log(args, log_file, 5, full_graph_description)

        params = [clemb, subw, sharp_subw, m_width]
        skips = [args.ignore_subw_vs_clemb, args.ignore_sharp_subw_vs_subw, args.ignore_m_width_vs_sharp_subw]
        interesting_graph = False
        for i in range(len(params)):
            if params[i]:
                for j in range(i + 1, len(params)):
                    if params[j]:
                        if abs(params[i] - params[j]) > EPSILON and not all(skips[i:j]):
                            interesting_graph = True
                        break
        if interesting_graph:
            interesting_graph_descriptions.append(full_graph_description)
            print_log(args, log_file, 1, f'\n\nINTERESTING GRAPH:\n\n{full_graph_description}\n\n')
            if args.save_pictures_of_interesting_graphs_to:
                draw_graph(graph, cur_index + 1, short_graph_description, args.save_pictures_of_interesting_graphs_to, nx, plt)
            if args.save_interesting_graphs_to:
                print_graph(graph, interesting_graphs_file)
                interesting_graphs_file.flush()

    if graphs:
        print_log(args, log_file, 2, 'reached 100%')
    print_log(args, log_file, 3, f'\n\nall {len(graphs)} graphs processed' if processed_graph_count == len(graphs) else f'\n\n{processed_graph_count} of {len(graphs)} graphs processed')
    print_log(args, log_file, 3, f'{len(interesting_graph_descriptions)} interesting graphs found:' if len(interesting_graph_descriptions) > 0 else 'no interesting graphs found')
    for graph_description in interesting_graph_descriptions:
        print_log(args, log_file, 3, graph_description)

    end_time = time()
    print_log(args, log_file, 2, f'total time: {round(end_time - start_time, 2)} sec.')
    if interesting_graphs_file:
        interesting_graphs_file.close()
    if log_file:
        log_file.close()


if __name__ == '__main__':
	main()
