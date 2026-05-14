#pragma once
#include <vector>
#include <algorithm>
#include <set>
#include <cassert>
#include "EdgeListGraph.h"
#include "get_mxclq.h"
#include <random>

#define EXTRA
#define CLIQUE

extern EdgeListGraph Graph;
extern int LB;
extern int K;
extern clock_t start_time;
extern clock_t end_time;


// 传入 method 参数 (1:度数正序 2:度数逆序 3:degen正序 4:degen逆序 5:随机)
// 默认值为 2，保持原有代码逻辑一致
int LB_with_fast_color(int method = 2);
// 找到每个点的最大独立集，用于ban边
void Max_clqs_of_nodes(std::vector<int>&);
// 找到包含u,v两点邻居的最大团，相当于找到包含边（u, v）的团。
std::vector<int> neg_clq_edges(int u,int v);
// 找到包含点u邻居的最大团。
std::vector<int> neg_clq_nodes(int u);
// 找到最后的最大团
std::vector<int> final_max_clique_nodes();

// 找到每条边的最大独立集，用于ban边
void Max_clqs_of_edges(std::vector<int>&);
// 找到 (u, v) 的公共邻居
int common_neighbor_num(int u, int v);
// 找到（公共邻居数目，边号）的映射
std::vector<std::pair<int, int>> common_neighbor_vector();
// 找公共邻居染色
std::vector<std::pair<int, int>> color_vector();
