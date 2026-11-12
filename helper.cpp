#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <array>
#include <cstdint>
#include <functional>
#include <unordered_set>

using namespace std;

extern "C" {
    // messages are passed to the caller (width-calculator.py),
    // which prints them to the console and to the log file depending on their importance
    typedef void (*log_callback)(int importance, const char* message);

    int* enumerate_minimal_tree_decompositions(int n, int m, const int* edge_list, log_callback log) {
        const int BIG_GRAPH = 12;
        if (n >= BIG_GRAPH) {
            log(3, ("WARNING: generating tree decompositions of a graph with " + to_string(n) + " vertices could take a long time").c_str());
        }
        vector<int> adj_mtx(n, 0);
        for (int i = 0; i < m; i++) {
            int a = edge_list[2 * i], b = edge_list[2 * i + 1];
            adj_mtx[a] |= (1 << b);
            adj_mtx[b] |= (1 << a);
        }
        vector<int> permutation(n);
        iota(permutation.begin(), permutation.end(), 0);
        int64_t iters = 1;
        for (int i = 1; i <= n; i++) {
            iters *= i;
        }
        int64_t checkpoint = iters / 100 + 1, curiter = 0, curpercent = 0;
        const uint64_t MOD = 1'000'000'007, MULT = (uint64_t(1) << n) * 3 + 17;
        auto vector_hash = [&](const vector<int>& vec) {
            uint64_t hash_value = 0;
            for (int val : vec) {
                hash_value = (hash_value * MULT + val + 1) % MOD;
            }
            return size_t(hash_value);
        };
        auto vector_equal = equal_to<vector<int>>();
        unordered_set<vector<int>, decltype(vector_hash), decltype(vector_equal)> tree_decompositions_set(0, vector_hash, vector_equal);
        do {
            // report progress after every completed percent of the permutations
            if (++curiter % checkpoint == 0) {
                curpercent++;
                if (n >= BIG_GRAPH && (curpercent % 5 == 0 || n > BIG_GRAPH)) {
                    log(4, ("\tprogress: " + to_string(curpercent) + "%").c_str());
                }
            }
            vector<int> cur_adj_mtx = adj_mtx;
            int alive_mask = (1 << n) - 1;
            vector<int> tree_decomposition;
            for (int v : permutation) {
                int bag = ((cur_adj_mtx[v] | (1 << v)) & alive_mask);
                tree_decomposition.push_back(bag);
                alive_mask ^= (1 << v);
                for (int u = 0; u < n; u++) {
                    if ((cur_adj_mtx[v] >> u) & 1) {
                        cur_adj_mtx[u] = ((cur_adj_mtx[u] | cur_adj_mtx[v]) ^ (1 << u));
                    }
                }
            }
            // remove bag repetitions
            sort(tree_decomposition.begin(), tree_decomposition.end());
            tree_decomposition.resize(unique(tree_decomposition.begin(), tree_decomposition.end()) - tree_decomposition.begin());
            // remove bags that are subsets of other bags
            vector<int> filtered_tree_decomposition;
            for (int bag : tree_decomposition) {
                bool superset_exists = false;
                for (int bag2 : tree_decomposition) {
                    if (bag != bag2 && (bag & bag2) == bag) {
                        superset_exists = true;
                        break;
                    }
                }
                if (!superset_exists) {
                    filtered_tree_decomposition.push_back(bag);
                }
            }
            tree_decompositions_set.insert(filtered_tree_decomposition);
        } while (next_permutation(permutation.begin(), permutation.end()));
        if (n >= BIG_GRAPH) {
            log(3, "Tree decompositions generated. Filtering of tree decompositions started.");
        }
        vector<vector<int>> tree_decompositions(tree_decompositions_set.begin(), tree_decompositions_set.end());

        // remove tree decompositions that are supersets of other tree decompositions
        vector<vector<int>> minimal_decompositions;
        for (const auto& tree_dec : tree_decompositions) {
            bool is_superdec = false;
            for (const auto& other_tree_dec : tree_decompositions) {
                if (tree_dec == other_tree_dec) {
                    continue;
                }
                bool is_superdec_cur = true;
                for (int bag : other_tree_dec) {
                    bool is_subbag = false;
                    for (int bag2 : tree_dec) {
                        if ((bag & bag2) == bag) {
                            is_subbag = true;
                            break;
                        }
                    }
                    if (!is_subbag) {
                        is_superdec_cur = false;
                        break;
                    }
                }
                if (is_superdec_cur) {
                    is_superdec = true;
                    break;
                }
            }
            if (!is_superdec) {
                minimal_decompositions.push_back(tree_dec);
            }
        }
        sort(minimal_decompositions.begin(), minimal_decompositions.end());
        if (n >= BIG_GRAPH) {
            log(3, "Filtering of tree decompositions completed");
        }

        vector<int> ans;
        ans.push_back(minimal_decompositions.size());
        for (const auto& decomposition : minimal_decompositions) {
            ans.push_back(decomposition.size());
            for (int bag : decomposition) {
                ans.push_back(bag);
            }
        }
        int* ans_array = new int[ans.size() + 1];
        ans_array[0] = ans.size() + 1;
        copy(ans.begin(), ans.end(), ans_array + 1);
        return ans_array;
    }

    void free_memory(int* ptr) {
        delete[] ptr;
    }
}
