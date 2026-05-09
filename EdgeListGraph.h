#pragma once
#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <map>
#include <cstring>
#include <ilcplex/ilocplex.h>

#include "MBitSet.h"
#include "utils.h"
#include "LinearHeap.h"
//using namespace std;


// cishu是全局变量
extern int cishu;
extern int LB;

static int fileSuffixPos(const char* filepath) {
	int j = strlen(filepath) - 1;
	while (filepath[j] != '.')
		j--;
	return j + 1;
}

static std::string integer_to_string(long long number) {
	std::vector<int> sequence;
	if (number == 0) sequence.push_back(0);
	while (number > 0) {
		sequence.push_back(number % 1000);
		number /= 1000;
	}

	char buf[5];
	std::string res;
	for (unsigned int i = sequence.size(); i > 0; i--) {
		if (i == sequence.size()) sprintf(buf, "%u", sequence[i - 1]);
		else sprintf(buf, ",%03u", sequence[i - 1]);
		res += std::string(buf);
	}
	return res;
}

struct MatrixGraph {
	int n;
	ept m;
	int* matrix;
	bool isEdge(int i, int j) {
		assert(i < n && j < n && i != j); // no selfloop
		if (matrix == nullptr) {
			fprintf(stderr, "Graph not initialized!\n");
			exit(1);
		}
		return matrix[i * n + j] != 0;
	}

};

// 定义带边编号的无向边（需要在读图的时候就定义好变量的）
typedef struct UEWT{
	int u;
	int v;
	int tag;
	//自定义去重, 只要 u 和 v 分别相同，那么他们相同。
	bool operator == (const UEWT& other) const {
        return u == other.u && v == other.v;
    }
}Undir_Edge_With_Tag;

// 根据边的顶点大小排序
inline bool cmp_with_node(const Undir_Edge_With_Tag & a, const Undir_Edge_With_Tag & b){
	if(a.u == b.u)
		if(a.v == b.v)
			return a.tag < b.tag;
		else
			return a.v < b.v;
	else
		return a.u < b.u;
}

// 根据边的编号排序
inline bool cmp_with_tag(const Undir_Edge_With_Tag & a, const Undir_Edge_With_Tag & b){
	if(a.tag == b.tag)
		if(a.u == b.u)
			return a.v < b.v;
		else
			return a.u < b.u;
	else
		return a.tag < b.tag;
}

// 定义带有编号的边，（结点号，边编号）
typedef struct EWT{
	int v;
	int tag;
}Edge_With_Tag;

// 根据顶点进行排序
inline bool cmp_with_v(const Edge_With_Tag & a, const Edge_With_Tag & b){
	return a.v < b.v;
}

struct EdgeListGraph {
	int n; // number of vertices
	ept m; // number of edges (undirected, so m = 2 * number of edges)
	ept* pstart = nullptr; //offset of neighbors of nodes
	ept* pend = nullptr;	//used in search
	// 存储当前结点 i 的度数
	int* degree = nullptr;
	// 存储根据结点度数的排序
	int* degree_rank = nullptr;
	// int* edges;
	// 这个 Edge_With_Tag 里面存储到达结点和边的编号
	Edge_With_Tag * edges = nullptr;
	// 存储虚禁用数组；
	int* ban_edge = nullptr;
	// 禁用顶点；
	int* ban_node = nullptr;
	// 从边的编号到两个顶点的映射
	std::pair<int, int>* edge_to_nodes;
	//
	int max_core;
	//oriented graph
	ept *pstart_o = nullptr;
	//oriented graph 
	int *edges_o = nullptr; 

	int *current_clique = nullptr;
	// size of max_core, used for sorting vertices w.r.t. #colors
	int *head = nullptr; 
	// size of max_core, used for sorting vertices w.r.t. #colors
	int *next = nullptr; 
	int *id = nullptr;
	int *mapping = nullptr;
	//adjacency matrix of an induced subgraph
	unsigned char *matrix = nullptr;

	std::vector<int> vs_buf;
	std::vector<int> color_buf;

	int *rid = nullptr;
	char *vis = nullptr;

	int mapping_n;
	int matrix_len;

	std::vector<int> contractions;

	long long branches;
	int max_depth;

	std::vector<std::pair<int, int> > changes;

	std::vector<int> del, degree_one, degree_two, degree_three;

	// 存储这个图中当前的clqs
	std::vector<std::vector<int>> clqs_edge;
	std::vector<std::vector<int>> clqs_node;
	// 用来记录z变量的编号
	std::vector<std::vector<int>> clqs_z;
	
	std::vector<int> max_clique;

	// 进行图染色进行ban边
	int* col = nullptr;
	// degen_order
	int* peel_sequence = nullptr;
	//
	int* core = nullptr;
	// 图的密度
	double R;
	//
	int* out_mapping;
	// 记录删除掉的边
	std::vector<int> global_delete_edges;

	EdgeListGraph() : n(0), m(0), max_core(0), pstart(nullptr), pend(nullptr), edges(nullptr), degree(nullptr), degree_rank(nullptr),
		ban_edge(nullptr), ban_node(nullptr), edge_to_nodes(nullptr), col(nullptr), peel_sequence(nullptr), core(nullptr),
		out_mapping(nullptr), rid(nullptr){}
	EdgeListGraph(const EdgeListGraph& g) : max_core(g.max_core), n(g.n), m(g.m) {
        pstart = g.pstart ? new ept[n + 1] : nullptr;
        pend   = g.pend   ? new ept[n + 1] : nullptr;
        edges  = g.edges  ? new Edge_With_Tag[m]     : nullptr;
		degree = g.degree ? new int[n] : nullptr;
		degree_rank = g.degree_rank ? new int[n] : nullptr;
		ban_edge = g.ban_edge ? new int[m / 2] : nullptr;
		edge_to_nodes = g.edge_to_nodes ? new std::pair<int, int> [m / 2] : nullptr;
		col = g.col ? new int [n] : nullptr;
		peel_sequence = g.peel_sequence ? new int [n] : nullptr;
		core = g.core ? new int [n] : nullptr;
		out_mapping = g.out_mapping ? new int [n] : nullptr;
		rid = g.rid ? new int [n] : nullptr;
        if (pstart) std::copy(g.pstart, g.pstart + n + 1, pstart);
        if (pend)   std::copy(g.pend,   g.pend + n + 1, pend);
        if (edges)  std::copy(g.edges,  g.edges + m,    edges);
		if (degree) std::copy(g.degree,  g.degree + n, degree);
		if(degree_rank) std::copy(g.degree_rank, g.degree_rank + n, degree_rank);
		if(ban_edge) std::copy(g.ban_edge, g.ban_edge + m / 2, ban_edge);
		if(ban_node) std::copy(g.ban_node, g.ban_node + n, ban_node);
		if(edge_to_nodes) std::copy(g.edge_to_nodes, g.edge_to_nodes + m / 2, edge_to_nodes);
		if(col) std::copy(g.col, g.col + n, col);
		if(peel_sequence) std::copy(g.peel_sequence, g.peel_sequence + n, peel_sequence);
		if(core) std::copy(g.core, g.core + n, core);
		if(out_mapping) std::copy(g.out_mapping, g.out_mapping + n, out_mapping);
		if(rid) std::copy(g.rid, g.rid + n, rid);
    }

    EdgeListGraph(EdgeListGraph&& g) noexcept
        : max_core(g.max_core), n(g.n), m(g.m), pstart(g.pstart), pend(g.pend), edges(g.edges), degree(g.degree), degree_rank(g.degree_rank), ban_edge(g.ban_edge), edge_to_nodes(g.edge_to_nodes),
			col(g.col), core(g.core){
        g.pstart = nullptr;
        g.pend = nullptr;
        g.edges = nullptr;
		g.degree = nullptr;
		g.degree_rank = nullptr;
		g.ban_edge = nullptr;
		g.ban_node = nullptr;
		g.edge_to_nodes = nullptr;
		g.col = nullptr;
		g.peel_sequence = nullptr;
		g.core = nullptr;
		g.out_mapping = nullptr;
		g.rid = nullptr;
    }
    EdgeListGraph& operator=(EdgeListGraph&& g) noexcept {
        if (this != &g) {
            delete[] pstart;
            delete[] pend;
            delete[] edges;
			delete[] degree;
			delete[] degree_rank;
			delete[] ban_edge;
			delete[] ban_node;
			delete[] edge_to_nodes;
			delete[] col;
			delete[] peel_sequence;
			delete[] core;
			delete[] out_mapping;
			delete[] rid;
            n = g.n;
            m = g.m;
			max_core = g.max_core;
            pstart = g.pstart;
            pend = g.pend;
            edges = g.edges;
			degree = g.degree;
			degree_rank = g.degree_rank;
			ban_edge = g.ban_edge;
			ban_node = g.ban_node;
			edge_to_nodes = g.edge_to_nodes;
			col = g.col;
			peel_sequence = g.peel_sequence;
			core = g.core;
			out_mapping = g.out_mapping;
			rid = g.rid;
            g.pstart = nullptr;
            g.pend = nullptr;
            g.edges = nullptr;
			g.degree = nullptr;
			g.degree_rank = nullptr;
			g.ban_edge = nullptr;
			g.edge_to_nodes = nullptr;
			g.ban_node = nullptr;
			g.col = nullptr;
			g.peel_sequence = nullptr;
			g.core = nullptr;
			g.out_mapping = nullptr;
			g.rid = nullptr;
        }
        return *this;
    }
	//release
	~EdgeListGraph() {
		if (pstart != nullptr) delete[] pstart;
		if (pend != nullptr) delete[] pend;
		if (edges != nullptr) delete[] edges;
		if (degree != nullptr) delete[] degree;
		if(degree_rank != nullptr) delete[] degree_rank;
		if(ban_edge != nullptr) delete[] ban_edge;
		if(ban_node != nullptr) delete[] ban_node;
		if(edge_to_nodes != nullptr) delete[] edge_to_nodes;
		if(col != nullptr) delete[] col;
		if(peel_sequence != nullptr) delete[] peel_sequence;
		if(core != nullptr) delete[] core;
		if(out_mapping != nullptr) delete[] out_mapping;
		if(rid != nullptr) delete[] rid;
 		pstart = nullptr;
		pend = nullptr;
		edges = nullptr;
		degree = nullptr;
		degree_rank = nullptr;
		ban_edge = nullptr;
		ban_node = nullptr;
		edge_to_nodes = nullptr;
		col = nullptr;
		peel_sequence = nullptr;
		core = nullptr;
		out_mapping = nullptr;
		rid = nullptr;
	}


	// 用边进行子图初始化
	void sub_init(std::vector<std::pair<int, int>> the_subedges){
		std::vector<Undir_Edge_With_Tag> temp_edge;
		int cnt = 0;
		for(auto [a, b] : the_subedges){
			temp_edge.push_back((Undir_Edge_With_Tag){a, b, cnt});
			temp_edge.push_back((Undir_Edge_With_Tag){b, a, cnt});
			cnt++;
		}
		sort(temp_edge.begin(), temp_edge.end(), cmp_with_node);
		// std::cerr << "The num of temp_edge(in sub structure):" << temp_edge.size() << std::endl;
		// 统计边数和顶点数
		n = -1;
		m = temp_edge.size();
		for(int i = 0; i < temp_edge.size(); i++){
			if(temp_edge[i].u > n)
				n = temp_edge[i].u;
			if(temp_edge[i].v > n)
				n = temp_edge[i].v;
			//std::cerr << "debug:The temp_edge(u, v) is " << "(" << temp_edge[i].u << ", " << temp_edge[i].v << ")" << " n is:" << n <<std::endl;
		}
		n++;

		// std::cerr << "The num of vertice(in sub structure):" << n << std::endl;
		if(pstart != nullptr)
			delete[] pstart;
		pstart = new ept[n+1];
		if(edges != nullptr)
			delete[] edges;
		edges = new Edge_With_Tag[m];
		if(edge_to_nodes != nullptr)
			delete[] edge_to_nodes;
		edge_to_nodes = new std::pair<int, int>[m / 2];
		// 处理好边
		int idx = 0;
		for(int i = 0; i < n; i++){
			pstart[i] = idx;
			while(idx < m && temp_edge[idx].u == i){
				edges[idx] = (Edge_With_Tag){temp_edge[idx].v, temp_edge[idx].tag};
				idx++;
			}
		}
		pstart[n] = m;

		// 完成边的映射
		for(int i = 0; i < n; i++){
			for(int j = pstart[i]; j < pstart[i+1]; j++){
				edge_to_nodes[edges[j].tag].first = std::min(i, edges[j].v);
				edge_to_nodes[edges[j].tag].second = std::max(i, edges[j].v);
			}
		}
	}

	void heuristic_max_clique_max_degree(int processed_threshold) {
		assert(max_clique.empty());
		int *head = new int[n];
		int *next = new int[n];
		int *degree = new int[n];

		int *vis = new int[n];
		memset(vis, 0, sizeof(int)*n);
		++cishu;
		if (cishu==36)
		{
			cishu = cishu * 2;
		}
		int max_degree = 0;
		for(int i = 0;i < n;i ++) head[i] = n;
		for(int i = 0;i < n;i ++) {
			degree[i] = pstart[i+1]-pstart[i];
			if(degree[i] > max_degree) max_degree = degree[i];
			next[i] = head[degree[i]];
			head[degree[i]] = i;
		}
		// std::cerr << "n is:" << n << std::endl;
		// std::cerr << "max_degree is:" << max_degree << std::endl;
		for(int processed_vertices = 0;max_degree >= max_clique.size()&&processed_vertices < processed_threshold;processed_vertices ++) {
			int u = n;
			while(max_degree >= max_clique.size()&&u == n) {
				for(int v = head[max_degree];v != n;) {
					// std::cerr << "v is:" << v << std::endl;
					int tmp = next[v];
					if(degree[v] == max_degree) {
						u = v;
						head[max_degree] = tmp;
						break;
					}
					else if(degree[v] >= max_clique.size()) {
						next[v] = head[degree[v]];
						head[degree[v]] = v;
					}
					v = tmp;
				}
				if(u == n) {
					head[max_degree] = n;
					-- max_degree;
				}
			}
			if(u == n) break;

			vis[u] = 1;
			for(int k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k].v]) -- degree[edges[k].v];

			std::vector<int> vs;
			for(int j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j].v]) vs.push_back(edges[j].v);

			std::vector<int> vs_deg(vs.size());
			for(int j = 0;j < vs.size();j ++) vis[vs[j]] = 2;
			for(int j = 0;j < vs.size();j ++) {
				int v = vs[j], d = 0;
				for(int k = pstart[v];k < pstart[v+1];k ++) {
					// std::cerr << "pstart:" << pstart[v] << " " << pstart[v+1] << std::endl;
					// std::cerr << "m is:" << m << std::endl;
					// std::cerr << "debug:edges[k].v" << ":" << edges[k].v << " k:" << k << std::endl;
					if(vis[edges[k].v] == 2) ++ d;
				}
				vs_deg[j] = d;
			}
			for(int j = 0;j < vs.size();j ++) vis[vs[j]] = 0;

			std::vector<int> res; res.push_back(u);
			int vs_size = vs.size();
			while(vs_size > 0&&res.size() + vs_size > max_clique.size()) {
				int idx = 0;
				for(int j = 1;j < vs_size;j ++) {
					if(vs_deg[j] > vs_deg[idx]) idx = j;
					else if(vs_deg[j] == vs_deg[idx]&&degree[vs[j]] > degree[vs[idx]]) idx = j;
				}
				u = vs[idx];

				int new_size = 0;
				for(int j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j].v]) vis[edges[j].v] = 2;
				for(int j = 0;j < vs_size;j ++) if(vis[vs[j]]) {
					if(j != new_size) std::swap(vs[new_size], vs[j]);
					vs_deg[new_size] = vs_deg[j];
					++ new_size;
				}
				for(int j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j].v] == 2) vis[edges[j].v] = 0;

				res.pb(u);
				for(int k = 0;k < new_size;k ++) vis[vs[k]] = k+2;
				for(int j = new_size;j < vs_size;j ++) {
					int v = vs[j];
					for(int k = pstart[v];k < pstart[v+1];k ++) {
						if(vis[edges[k].v] >= 2) -- vs_deg[vis[edges[k].v]-2];
					}
				}
				for(int k = 0;k < new_size;k ++) vis[vs[k]] = 0;

				vs_size = new_size;
			}

			if(res.size() > max_clique.size()) max_clique = res;
		}

		delete[] vis;
		delete[] head;
		delete[] next;
		delete[] degree;
	#ifndef NDEBUG
		//printf("*** Heuristic clique size: %lu, time: %s (microseconds)\n", max_clique.size(), integer_to_string(t.elapsed()).c_str());
		printf("There is a debug!\n");
	#endif
	}

	int coloring_adj_list(const int *vs, const int vs_size, const int original_size, int *color, char *vis, const int start_idx, const int start_color) {
		assert(start_color <= vs_size&&original_size >= vs_size);
		for(int i = 0;i < original_size;i ++) color[vs[i]] = n;
		for(int i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;

		int max_color = 0;
		for(int i = vs_size - start_color;i > start_idx;i --) {
			int u = vs[i-1];
			for(int j = pstart[u];j < pend[u];j ++) {
				int c = color[edges[j].v];
				if(c != n) vis[c] = 1;
			}
			for(int j = 0;;j ++) if(!vis[j]) {
				color[u] = j;
				if(j > max_color) max_color = j;
				break;
			}
			for(int j = pstart[u];j < pend[u];j ++) {
				int c = color[edges[j].v];
				if(c != n) vis[c] = 0;
			}
		}

		return max_color + 1;
	}

	void get_higher_neighbors(const int u, int &vs_size, std::vector<int> &vs_buf, std::vector<int> &color_buf, const ept *pstart_o, const int *edges_o) {
		vs_size = 0;
		for(int j = pstart_o[u];j < pstart_o[u+1];j ++) {
			if(vs_buf.size() == vs_size) {
				vs_buf.pb(edges_o[j]);
				color_buf.pb(0);
			}
			else vs_buf[vs_size] = edges_o[j];
			++ vs_size;
		}
	}

	//Bintld the oriented graph
	void shrink_graph(int *&peel_sequence, int *&core, int *&color, int *&out_mapping, ept *&pstart_o, int *&edges_o) {
		int *rid = new int[n];
		for(int i = 0;i < n;i ++) rid[peel_sequence[i]] = i;

		int clique_size = max_clique.size();
		assert(pstart_o == nullptr);
		pstart_o = new ept[n+1];
		pstart_o[0] = 0;
		int cnt = 0;
		for(int i = 0;i < n;i ++) {
			pstart_o[i+1] = pstart_o[i];
			if(core[i] < clique_size) continue;
			++ cnt;
			ept &pos = pstart_o[i+1];
			for(ept j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j].v] > rid[i]) edges[pos ++] = edges[j];
		}

		assert(out_mapping == nullptr);
		out_mapping = new int[cnt];

		cnt = 0;
		for(int i = 0;i < n;i ++) if(core[i] >= clique_size) {
			out_mapping[cnt] = i;
			core[cnt] = core[i];
			color[cnt] = color[i];
			rid[i] = cnt ++;
		}
		int x = n - cnt;
		for(int i = x;i < n;i ++) peel_sequence[i-x] = rid[peel_sequence[i]];
		assert(edges_o == nullptr);
		edges_o = new int[pstart_o[n]];
		for(int i = 0;i < pstart_o[n];i ++) edges_o[i] = rid[edges[i].v];
		for(int i = 0;i < cnt;i ++) pstart_o[i] = pstart_o[out_mapping[i]];
		pstart_o[cnt] = pstart_o[n];
		n = cnt;

		delete[] rid;

		int *t_peel_sequence = new int[n];
		memcpy(t_peel_sequence, peel_sequence, sizeof(int)*n);
		delete[] peel_sequence; peel_sequence = t_peel_sequence;

		int *t_core = new int[n];
		memcpy(t_core, core, sizeof(int)*n);
		delete[] core; core = t_core;

		int *t_color = new int[n];
		memcpy(t_color, color, sizeof(int)*n);
		delete[] color; color = t_color;

		ept *t_pstart = new ept[n+1];
		memcpy(t_pstart, pstart_o, sizeof(ept)*(n+1));
		delete[] pstart_o; pstart_o = t_pstart;
	#ifndef NDEBUG
		printf("\tReduced graph size: |V|=%s, |E|=%s (undirected)\n", integer_to_string(cnt).c_str(), integer_to_string(pstart_o[n]).c_str());
	#endif
	}

	int color_bound(const int *vs, const int vs_size, const int *color, char *vis) {
	#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) assert(!vis[color[vs[i]]]);
	#endif
		int color_bound = 0;
		for(int i = 0;i < vs_size;i ++) if(!vis[color[vs[i]]]) {
			vis[color[vs[i]]] = 1;
			++ color_bound;
		}
		for(int i = 0;i < vs_size;i ++) vis[color[vs[i]]] = 0;
		return color_bound;
	}

	//greedily enlarge max_clique by including u if feasible
	char greedy_extend(const int u, const char *vis, char print) {
		for(auto v: max_clique) if(!vis[v]) return 0;
		max_clique.pb(u);
		if(print) printf("greedy_extend finds clique of size: %lu\n", max_clique.size());
		return 1;
	}

	// construct induced subgraph from pstart_o and edges_o
	void construct_induced_subgraph(const int *vs, const int vs_size, char *vis, int *degree, const ept *pstart_o, const int *edges_o) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif
		for(int j = 0;j < vs_size;j ++) degree[vs[j]] = 0;
		for(int j = 0;j < vs_size;j ++) {
			int v = vs[j];
			for(int k = pstart_o[v];k < pstart_o[v+1];k ++) {
				if(vis[edges_o[k]]) {
					++ degree[v];
					++ degree[edges_o[k]];
				}
			}
		}

		pstart[vs[0]] = 0;
		for(int j = 1;j < vs_size;j ++) pstart[vs[j]] = pstart[vs[j-1]] + degree[vs[j-1]];
		for(int j = 0;j < vs_size;j ++) pend[vs[j]] = pstart[vs[j]];
		for(int j = 0;j < vs_size;j ++) {
			int v = vs[j];
			for(int k = pstart_o[v];k < pstart_o[v+1];k ++) {
				int w = edges_o[k];
				if(vis[w]) {
					edges[pend[v] ++].v = w;
					edges[pend[w] ++].v = v;
				}
			}
		}
	}

	void kcore_reduction(int *vs, int &vs_size, char *vis, int *degree, const int K, int *queue) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif

		int queue_n = 0;
		for(int i = 0;i < vs_size;i ++) if(degree[vs[i]] < K) {
			vis[vs[i]] = 0;
			queue[queue_n ++] = vs[i];
		}
		for(int i = 0;i < queue_n;i ++) {
			int u = queue[i];
			for(int j = pstart[u];j < pend[u];j ++) if(vis[edges[j].v]) {
				if((-- degree[edges[j].v]) == K-1) {
					vis[edges[j].v] = 0;
					queue[queue_n ++] = edges[j].v;
				}
			}
		}
		int new_vs_size = 0;
		for(int i = 0;i < vs_size;i ++) if(vis[vs[i]]) {
			if(new_vs_size != i) std::swap(vs[new_vs_size], vs[i]);
			++ new_vs_size;
		}
		vs_size = new_vs_size;
	}

	// construct matrix from pstart_o and edges_o for vertices in vs
	void construct_matrix(int *vs, const int vs_size, int *mapping, char *vis, int *rdegree, const ept *pstart_o, const int *edges_o) {
		assert(rid != nullptr&&mapping != nullptr&&matrix != nullptr);
		assert(vs_size <= max_core);

		mapping_n = vs_size;
#ifdef _BITSET_
		matrix_len = (vs_size+7)/8;
#else
		matrix_len = vs_size;
#endif
		for(int j = 0;j < vs_size;j ++) {
			vis[vs[j]] = 1;
			rid[vs[j]] = j;
			mapping[j] = vs[j];
			vs[j] = j;
		}

		memset(matrix, 0, sizeof(unsigned char)*mapping_n*matrix_len);
		for(int j = 0;j < mapping_n;j ++) rdegree[j] = mapping_n - 1;
		for(int j = 0;j < mapping_n;j ++) for(int k = pstart_o[mapping[j]];k < pstart_o[mapping[j]+1];k ++) if(vis[edges_o[k]]) {
			int v = rid[edges_o[k]];
			assert(v >= 0&&v < mapping_n);
			set_bit(matrix + j*matrix_len, v);
			set_bit(matrix + v*matrix_len, j);
			-- rdegree[j]; -- rdegree[v];
		}
		for(int j = 0;j < mapping_n;j ++) set_bit(matrix + j*matrix_len, j);

		for(int j = 0;j < mapping_n;j ++) vis[mapping[j]] = 0;
	}

	void remove_vertex_idx(int *vs, int idx, int &vs_size, int *rdegree) {
		assert(idx >= 0&&idx < vs_size);

		unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
		vs[idx] = vs[-- vs_size];
		for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) -- rdegree[vs[j]];
	}

	void store_a_larger_clique(const int clique_size, const char *info, char print) {
		print = 0;
		assert(max_clique.size() < clique_size);
		while(max_clique.size() < clique_size) max_clique.pb(0);

		int contract_size = contractions.size();
		for(int k = clique_size-1;k > 0;k --) {
			if(current_clique[k] == n) {
				assert(contract_size > 0);
				if(contractions[contract_size-1] == 1) {
					contract_size -= 4;
					if(vis[contractions[contract_size]]) max_clique[k] = contractions[contract_size+1];
					else max_clique[k] = contractions[contract_size+2];
				}
				else if(contractions[contract_size-1] == 2) {
					contract_size -= 5;
					assert(!vis[contractions[contract_size]]||!vis[contractions[contract_size+1]]);
					if(!vis[contractions[contract_size]]&&!vis[contractions[contract_size+1]]) max_clique[k] = contractions[contract_size+3];
					else max_clique[k] = contractions[contract_size+2];
				}
				else if(contractions[contract_size-1] == 3) {
					contract_size -= 5;
					char in1 = (vis[contractions[contract_size]] != 0);
					char in2 = (vis[contractions[contract_size+1]] != 0);
					char in3 = (vis[contractions[contract_size+2]] != 0);
					assert(in1+in2+in3 < 3);
					if(in1+in2+in3 == 0) max_clique[k] = contractions[contract_size+3];
					else if(in1+in2+in3 == 1) {
						if(in1) max_clique[k] = contractions[contract_size+1];
						else if(in2) max_clique[k] = contractions[contract_size+2];
						else max_clique[k] = contractions[contract_size];
					}
					else {
						if(!in1) max_clique[k] = contractions[contract_size];
						else if(!in2) max_clique[k] = contractions[contract_size+1];
						else max_clique[k] = contractions[contract_size+2];
					}
				}
				else {
					printf("WA!\n");
				}
			}
			else max_clique[k] = current_clique[k];
			vis[max_clique[k]] = 1;
		}
		assert(contract_size == 0);
		for(int k = 1;k < max_clique.size();k ++) vis[max_clique[k]] = 0;
		for(int k = 1;k < max_clique.size();k ++) max_clique[k] = mapping[max_clique[k]];
		max_clique[0] = current_clique[0];
		if(print) printf("%s finds clique of size: %lu\n", info, max_clique.size());

#ifndef NDEBUG
		int total_edges = max_clique.size();
		for(int i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 1;
		for(int i = 0;i < max_clique.size();i ++) for(ept j = pstart_o[max_clique[i]];j < pstart_o[max_clique[i]+1];j ++) {
			if(vis[edges_o[j]]) total_edges += 2;
		}
		if(total_edges != max_clique.size()*max_clique.size()) printf("WA! Not a clique! %s\n", info);
		for(int i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 0;
#endif
	}

	void degree_one_two_reduction_with_folding_matrix(int &current_clique_size, int *vs, int &vs_size, int *rdegree, std::vector<int> &contractions) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) {
			int rd = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
			assert(test_bit(t_matrix, vs[i]));
			assert(rdegree[vs[i]] == rd);
		}
#endif

		char changed = 1;
		while(changed) {
			changed = 0;
			for(int i = 0;i < vs_size;) {
				int rd = rdegree[vs[i]];
				if(rd == 0) {
					changed = 1;
					current_clique[current_clique_size ++] = vs[i];
					vs[i] = vs[-- vs_size];
				}
				else if(vs_size - rd + current_clique_size <= max_clique.size()) {
					changed = 1;
					remove_vertex_idx(vs, i, vs_size, rdegree);
				}
				else if(rd == 1) {
					changed = 1;
					current_clique[current_clique_size ++] = vs[i];
					unsigned char *t_matrix = matrix + vs[i]*matrix_len;
					vs[i] = vs[-- vs_size];
					for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
						remove_vertex_idx(vs, j, vs_size, rdegree);
						break;
					}
				}
				else if(rd == 2) {
					changed = 1;
					int u = vs[i];
					unsigned char *t_matrix = matrix + u*matrix_len;
					vs[i] = vs[-- vs_size];
					int idx1 = n, idx2 = n;
					for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
						if(idx1 == n) idx1 = j;
						else {
							idx2 = j;
							break;
						}
					}
					assert(idx1 != n&&idx2 != n&&idx1 < idx2);

					if(!test_bit(matrix+vs[idx1]*matrix_len, vs[idx2])) {
						current_clique[current_clique_size ++] = u;
						remove_vertex_idx(vs, idx2, vs_size, rdegree);
						remove_vertex_idx(vs, idx1, vs_size, rdegree);
					}
					else {
						current_clique[current_clique_size ++] = n;
						-- rdegree[vs[idx1]];
						contractions.pb(vs[idx1]); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1);
						unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
						unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
						vs[idx2] = vs[-- vs_size];

						for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix2, vs[j])) {
							if(test_bit(t_matrix1, vs[j])) {
								reverse_bit(t_matrix1, vs[j]);
								reverse_bit(matrix + vs[j]*matrix_len, vs[idx1]);
								++ rdegree[vs[idx1]];
							}
							else -- rdegree[vs[j]];
						}
					}
				}
				else ++ i;
			}
		}

#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			int rd = 0;
			for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
			assert(rd == rdegree[vs[i]]);
			assert(vs_size - 1 - rd + current_clique_size >= max_clique.size());
			assert(rd >= 2);
		}
#endif

		if(current_clique_size > max_clique.size()) store_a_larger_clique(current_clique_size, "degree_one_two_with_folding", 0);
	}

	int degeneracy_maximal_clique_matrix(int current_clique_size, int *vs, const int vs_size, int *degree, char heuristic_gen, char print) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) {
			int d = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(int j = 0;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
			assert(test_bit(t_matrix, vs[i]));
			assert(degree[vs[i]]+1 == d);
		}
#endif

		int start_color = 0;
		for(int j = 0;j < vs_size;j ++) {
			int min_idx = j;
			for(int k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
			if(min_idx != j) std::swap(vs[min_idx], vs[j]);
			int v = vs[j];
			if(degree[v] + 1 + j == vs_size) {
				for(int k = j+1;k < vs_size;k ++) degree[vs[k]] = vs_size-1-k;
				start_color = degree[v];
				assert(start_color > 0&&start_color < vs_size);
				if(heuristic_gen&&degree[v] + 1 + current_clique_size > max_clique.size()) {
					for(int k = j;k < vs_size;k ++) current_clique[current_clique_size ++] = vs[k];
					store_a_larger_clique(current_clique_size, "degen", print);
				}
				break;
			}
			unsigned char *t_matrix = matrix + v*matrix_len;
			for(int k = j+1;k < vs_size;k ++) if(test_bit(t_matrix, vs[k])) -- degree[vs[k]];
		}
		return start_color;
	}

	//coloring a graph that is represented by matrix
	//return the number of colors used
	int coloring_matrix(const int *vs, const int vs_size, int *color, char *vis, const int start_idx, const int start_color) {
		assert(start_color > 0&&start_color <= vs_size);
		for(int i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;
		int max_color = start_color-1;

		for(int i = vs_size - start_color;i > start_idx;i --) {
			int u = vs[i-1];
			unsigned char *t_matrix = matrix + u*matrix_len;
			for(int j = i;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) vis[color[vs[j]]] = 1;
			for(int j = 0;;j ++) if(!vis[j]) {
				color[u] = j;
				if(j > max_color) max_color = j;
				break;
			}
			for(int j = i;j < vs_size;j ++) vis[color[vs[j]]] = 0;
		}

		return max_color + 1;
	}

		int degeneracy_maximal_clique_adj_list(const int current_clique_size, int *vs, const int vs_size, char *vis, int *degree, ListLinearHeap *heap) {
		assert(vs_size > 0);

		for(int j = 0;j < vs_size;j ++) vis[vs[j]] = 1;

		char sparse = 0;
		int total_edges = 0;
		for(int i = 0;i < vs_size;i ++) total_edges += degree[vs[i]];
		if(total_edges*10 < vs_size*(vs_size-1)) sparse = 1;

		if(sparse) heap->init(vs_size, vs_size-1, vs, degree);

		int start_color = 0;
		for(int j = 0;j < vs_size;j ++) {
			int v, key;
			if(sparse) {
				heap->pop_min(v, key);
				vs[j] = v;
			}
			else {
				int min_idx = j;
				for(int k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
				if(min_idx != j) std::swap(vs[min_idx], vs[j]);
				v = vs[j]; key = degree[v];
			}
			if(key + j + 1 == vs_size) {
				start_color = key + 1;
				if(sparse) {
					int new_size = j+1;
					heap->get_ids(vs, new_size);
					assert(new_size == vs_size);
				}
				if(key + 1 + current_clique_size > max_clique.size()) {
					//printf("Find clique of size %u after search %u egos\n", degree[v] + 2, n-i);
					max_clique.clear();
					max_clique.reserve(key + 1 + current_clique_size);
					for(int k = j;k < vs_size;k ++) max_clique.pb(vs[k]);
					for(int k = current_clique_size;k > 0;k --) max_clique.pb(current_clique[k-1]);

	#ifndef NDEBUG
					int total_edges = max_clique.size();
					for(int i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 1;
					for(int i = 0;i < max_clique.size();i ++) for(ept j = pstart_o[max_clique[i]];j < pstart_o[max_clique[i]+1];j ++) {
						if(vis[edges_o[j]]) total_edges += 2;
					}
					if(total_edges != max_clique.size()*max_clique.size()) printf("WA! Not a clique\n");
					for(int i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 0;
	#endif
				}
				break;
			}
			vis[v] = 0;
			if(sparse) {
				for(int k = pstart[v];k < pend[v];k ++) if(vis[edges[k].v] == 1) heap->decrement(edges[k].v, 1);
			}
			else {
				for(int k = pstart[v];k < pend[v];k ++) if(vis[edges[k].v] == 1) -- degree[edges[k].v];
			}
		}
		for(int j = 0;j < vs_size;j ++) vis[vs[j]] = 0;

		return start_color;
	}

	int color_bound(const int *vs, const int vs_size, const int *mapping, const int *color, char *vis) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) assert(!vis[color[mapping[vs[i]]]]);
#endif
		int color_bound = 0;
		for(int i = 0;i < vs_size;i ++) if(!vis[color[mapping[vs[i]]]]) {
			vis[color[mapping[vs[i]]]] = 1;
			++ color_bound;
		}
		for(int i = 0;i < vs_size;i ++) vis[color[mapping[vs[i]]]] = 0;
		return color_bound;
	}

	//construct a better maximal clique and a tighter upper bound by considering all ego-networks
	int ego_degen(const int *peel_sequence, const int *core, const int *color_, int *local_UBs, const int UB) {
	#ifndef NDEBUG
		for(int i = 0;i < n;i ++) assert(!vis[i]);
		assert(max_clique.size() >= 2&&max_clique.size() <= n);
	#endif

		//Timer t;

		int max_local_UB = 0, initial_size = max_clique.size();
		int *queue = new int[max_core];
		ListLinearHeap *heap = new ListLinearHeap(n, max_core);
		heap->init(0, 0, nullptr, nullptr);
		for(int i = n;i > 0;i --) {
			int u = peel_sequence[i-1];
			local_UBs[i-1] = 0;

			if(n-i+1 <= max_clique.size()) continue;
			if(core[u] < max_clique.size()) break;

			//get N^+(u)
			int vs_size = 0;
			get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);
			assert(vs_size <= max_core);

			//color-based prune
			if(vs_size < max_clique.size()||color_bound(vs_buf.data(), vs_size, color_, vis) < max_clique.size()) continue;

			for(int j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

			//test max_clique \cup \{u\}
			if(max_clique.size() > initial_size) greedy_extend(u, vis, 0);

			//construct G[N^+(u)] and reduce by k-core
			int old_size = vs_size, original_size = old_size;
			construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
			kcore_reduction(vs_buf.data(), vs_size, vis, degree, max_clique.size()-1, queue);
			for(int j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
			if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color_, vis) < max_clique.size()) continue;

			char kernel = 0; // by default set as 0
			int current_clique_size = 1; current_clique[0] = u;
			if(kernel) {
				//construct matrix
				int *rdegree = degree;
				construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

				contractions.clear();
				old_size = vs_size;
				degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
				if(max_clique.size() >= UB) break;

				//if(vs_size + current_clique_size <= max_clique.size()) continue;
				if(vs_size < old_size&&current_clique_size + color_bound(vs_buf.data(), vs_size, mapping, color_, vis) <= max_clique.size()) continue;

				//degeneracy-based maximal clique
				//sort(vs_buf.begin(), vs_buf.begin()+vs_size);
				for(int j = 0;j < vs_size;j ++) degree[vs_buf[j]] = vs_size - 1 - rdegree[vs_buf[j]];
				int start_color = degeneracy_maximal_clique_matrix(current_clique_size, vs_buf.data(), vs_size, degree, 1, 0);
				if(max_clique.size() >= UB) break;

				//color-based upper bound
				//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
				local_UBs[i-1] = current_clique_size + coloring_matrix(vs_buf.data(), vs_size, rid, vis, 0, start_color);
			}
			else {
				//degeneracy-based maximal clique
				int start_color = degeneracy_maximal_clique_adj_list(current_clique_size, vs_buf.data(), vs_size, vis, degree, heap);

				if(max_clique.size() >= UB) break;

				//color-based upper bound
				local_UBs[i-1] = current_clique_size + coloring_adj_list(vs_buf.data(), vs_size, original_size, rid, vis, 0, start_color);
			}

			if(local_UBs[i-1] > max_local_UB) max_local_UB = local_UBs[i-1];
		}

		delete[] queue;
		delete heap;

		if(max_local_UB > UB) max_local_UB = UB;

		int new_UB = max_clique.size();
		if(max_local_UB > new_UB) new_UB = max_local_UB;
	#ifndef NDEBUG
		//printf("*** ego_degen clique size: %lu, UB: %u, Time: %s (microseconds)\n", max_clique.size(), new_UB, integer_to_string(t.elapsed()).c_str());
		printf("There is a debug!\n");
	#endif
		return new_UB;
	}

	// degeneracy based clique, if opt also greedy search by maximum degree
	// return an upper bound of the maximum clique
	int degeneracy_maximal_clique_adjacency_list(int *peel_sequence, int *core, int *color, char opt, char greedy_extend = 0) {
		if(opt) heuristic_max_clique_max_degree(10);
		int threshold = max_clique.size();

		//Timer t;

		int *id_s = peel_sequence;
		int *degree = new int[n];
		char *vis = new char[n];
		memset(vis, 0, sizeof(char)*n);

		if(pend == nullptr) pend = new ept[n];
		for(int i = 0;i < n;i ++) {
			pend[i] = pstart[i+1];
			degree[i] = pend[i] - pstart[i];
		}

		int queue_n = 0, new_size = 0;
		for(int i = 0;i < n;i ++) if(degree[i] < threshold) id_s[queue_n ++] = i;
		for(int i = 0;i < queue_n;i ++) {
			int u = id_s[i]; degree[u] = 0;
			for(int j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j].v] > 0) {
				if((degree[edges[j].v] --) == threshold) id_s[queue_n ++] = edges[j].v;
			}
		}
		for(int i = 0;i < n;i ++) {
			if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
			else {
				vis[i] = 1;
				core[i] = 0;
			}
		}
		assert(queue_n + new_size == n);

		int UB = n;
		if(new_size == 0) UB = max_clique.size();
		else {
			ListLinearHeap *heap = new ListLinearHeap(n, new_size-1);
			heap->init(new_size, new_size-1, id_s+queue_n, degree);
			max_core = 0;
			std::vector<int> res;
			for(int i = 0;i < new_size;i ++) {
				int u, key;
				heap->pop_min(u, key);
				if(key > max_core) max_core = key;
				core[u] = max_core;
				id_s[queue_n + i] = u;
				if(key + i + 1 == new_size) {
					int x_size = i+1;
					heap->get_ids(id_s+queue_n, x_size);
					assert(x_size == new_size);
					for(int j = i;j < new_size;j ++) {
						core[id_s[queue_n+j]] = max_core;
						res.pb(id_s[queue_n+j]);
					}
					break;
				}
				vis[u] = 1;

				for(int j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j].v] == 0) {
					heap->decrement(edges[j].v, 1);
				}
			}
			delete heap;
	#ifndef NDEBUG
			//printf("*** Degeneracy clique size: %lu, max_core: %u, Time: %s (microseconds)\n", res.size(), max_core, integer_to_string(t.elapsed()).c_str());
			printf("There is a debug!\n");
	#endif
			if(res.size() > max_clique.size()) max_clique = res;

			if(max_clique.size() == max_core+1) UB = max_core + 1;
			else{
				memset(vis, 0, sizeof(char)*n);
				int start_idx = 0;
				while(start_idx < n&&core[id_s[start_idx]] < max_clique.size()) ++ start_idx;
				int num_color = coloring_adj_list(id_s, n, n, color, vis, start_idx, 0);
				assert(num_color <= UB);
				UB = num_color;
				if(max_clique.size() < UB&&greedy_extend) {
					for(int i = 0;i < res.size();i ++) vis[res[i]] = 1;
					for(int i = n-res.size();i > 0;i --) {
						int u = id_s[i-1], cnt = 0;
						if(core[u] < res.size()) break;
						for(int j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j].v]) ++ cnt;
						if(cnt == res.size()) {
							res.pb(u);
							vis[u] = 1;
						}
					}
	#ifndef NDEBUG
					//printf("*** Degen_greedy_extend clique size: %lu, num_colors: %u, Time: %s (microseconds)\n", res.size(), num_color, integer_to_string(t.elapsed()).c_str());
					printf("There is a debug!\n");
	#endif
				}
				if(res.size() > max_clique.size()) max_clique = res;
			}
		}
		memset(vis, 0, sizeof(char)*n);

		delete[] degree;
		delete[] vis;

		return UB;
	}

	void put_into_one_vector(int &current_clique_size, int &i, int *vs, int &vs_size, bool check, const int rd, int *rid) {
		if(vs_size - rd + current_clique_size <= max_clique.size()) {
			del.pb(vs[i]);
			//printf("insert %u\n", vs[i]);
		}
		else if(check&&rd <= 3) {
			switch(rd) {
			case 0: current_clique[current_clique_size ++] = vs[i];
					vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
					-- i;
					break;
			case 1: degree_one.pb(vs[i]);
					break;
			case 2: degree_two.pb(vs[i]);
					break;
			case 3: degree_three.pb(vs[i]);
			}
		}
	}

	void put_into_one_vector_eq(int &current_clique_size, int &i, int *vs, int &vs_size, bool check, const int rd, int *rid) {
		if(check) {
			if(rd <= 3&&vs_size - rd + current_clique_size > max_clique.size()) {
				switch(rd) {
				case 0: current_clique[current_clique_size ++] = vs[i];
						//printf("added %u to maximum clique\n", vs[i]);
						vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
						-- i;
						break;
				case 1: degree_one.pb(vs[i]);
					break;
				case 2: degree_two.pb(vs[i]);
					break;
				case 3: degree_three.pb(vs[i]);
				}
			}
		}
		else {
			if(vs_size - rd + current_clique_size == max_clique.size()) {
				del.pb(vs[i]);
				//printf("insert %u\n", vs[i]);
			}
		}
	}

	void degree_one_two_three_reduction_with_folding_matrix(int &current_clique_size, int *vs, int &vs_size, int *rdegree, int *rid) {
		assert(del.empty()&&degree_one.empty()&&degree_two.empty()&&degree_three.empty());

		for(int i = 0;i < vs_size;i ++) {
			rid[vs[i]] = i;
			put_into_one_vector(current_clique_size, i, vs, vs_size, true, rdegree[vs[i]], rid);
		}

		while(!del.empty()||!degree_one.empty()||!degree_two.empty()||!degree_three.empty()) {
			while(!del.empty()) {
				int u = del.back(); del.pop_back();
				rdegree[u] = 0;
				//printf("u: %u, rid[u]: %u, vs[rid[u]]: %u, vs_size: %u\n", u, rid[u], vs[rid[u]], vs_size);
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				int idx = rid[u];
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
				unsigned char *t_matrix = matrix + u*matrix_len;
				for(int i = 0;i < vs_size;i ++) {
					int &rd = rdegree[vs[i]];
					int old_rd = rd;
					if(!test_bit(t_matrix, vs[i])) -- rd;
					put_into_one_vector_eq(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}

#ifndef NDEBUG
				for(int i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
			while(del.empty()&&!degree_one.empty()) {
				int u = degree_one.back(); degree_one.pop_back();
				if(rdegree[u] != 1) continue;

				current_clique[current_clique_size ++] = u;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				int idx = rid[u]; rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				unsigned char *t_matrix = matrix + u*matrix_len; idx = n;
				for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					idx = j;
					break;
				}
				assert(idx != n);

				t_matrix = matrix + vs[idx]*matrix_len;
				rdegree[vs[idx]] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
				for(int i = 0;i < vs_size;i ++) {
					int &rd = rdegree[vs[i]];
					int old_rd = rd;
					if(!test_bit(t_matrix, vs[i])) -- rd;
					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}

#ifndef NDEBUG
				for(int i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
			while(del.empty()&&degree_one.empty()&&!degree_two.empty()) {
				int u = degree_two.back(); degree_two.pop_back();
				if(rdegree[u] != 2) continue;

				unsigned char *t_matrix = matrix + u*matrix_len;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				int idx = rid[u]; rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				int idx1 = n, idx2 = n;
				for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					if(idx1 == n) idx1 = j;
					else {
						idx2 = j;
						break;
					}
				}
				assert(idx1 != n&&idx2 != n&&idx1 < idx2);

				unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
				unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
				if(!test_bit(t_matrix1, vs[idx2])) {
					current_clique[current_clique_size ++] = u;
					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					rdegree[vs[idx1]] = 0;
					vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

					for(int i = 0;i < vs_size;i ++) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;
						if(!test_bit(t_matrix1, vs[i])) -- rd;
						if(!test_bit(t_matrix2, vs[i])) -- rd;

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
				}
				else {
					current_clique[current_clique_size ++] = n;
					int v = vs[idx1], old_degree_v = rdegree[v];
					-- rdegree[v];
					contractions.pb(v); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1); //type-1 contraction

					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;

					for(int i = 0;i < vs_size;i ++) if(vs[i] != v) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;

						if(!test_bit(t_matrix2, vs[i])) {
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								++ rdegree[v];
								changes.pb(std::make_pair(v, vs[i]));
							}
							else -- rd;
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
				}

#ifndef NDEBUG
				for(int i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}

			while(del.empty()&&degree_one.empty()&&degree_two.empty()&&!degree_three.empty()) {
				int u = degree_three.back(); degree_three.pop_back();
				if(rdegree[u] != 3) continue;

				unsigned char *t_matrix = matrix + u*matrix_len;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				int idx = rid[u];
				//std::swap(vs[idx], vs[vs_size-1]); rid[vs[idx]] = idx; rid[vs[vs_size-1]] = vs_size - 1;
				rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				int idx1 = n, idx2 = n, idx3 = n;
				for(int j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					if(idx1 == n) idx1 = j;
					else if(idx2 == n) idx2 = j;
					else {
						idx3 = j;
						break;
					}
				}
				assert(idx1 != n&&idx2 != n&&idx3 != n&&idx1 < idx2&&idx2 < idx3);

				unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
				unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
				unsigned char *t_matrix3 = matrix + vs[idx3]*matrix_len;
				char connected_1 = 0, connected_2 = 0, connected_3 = 0;
				if(test_bit(t_matrix1, vs[idx2])) connected_1 = 1;
				if(test_bit(t_matrix1, vs[idx3])) connected_2 = 1;
				if(test_bit(t_matrix2, vs[idx3])) connected_3 = 1;
				char total_connected = connected_1 + connected_2 + connected_3;
				if(total_connected == 0) { //isolation reduction
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = u;
					rdegree[vs[idx3]] = 0;
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					rdegree[vs[idx1]] = 0;
					vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

					for(int i = 0;i < vs_size;i ++) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;
						if(!test_bit(t_matrix1, vs[i])) -- rd;
						if(!test_bit(t_matrix2, vs[i])) -- rd;
						if(!test_bit(t_matrix3, vs[i])) -- rd;

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
				}
				else if(total_connected == 1) {
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;
					//the following ensures that vs[idx2] is the vertex to be deleted
					if(connected_3) {
						std::swap(vs[idx1], vs[idx2]);
						std::swap(t_matrix1, t_matrix2);
						rid[vs[idx1]] = idx1;
					}
					else if(connected_1) {
						std::swap(vs[idx2], vs[idx3]);
						std::swap(t_matrix2, t_matrix3);
					}

					int v = vs[idx1], old_degree_v = rdegree[v];
					rdegree[v] -= 2;
					contractions.pb(v); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(1); //type-1 contraction

					rdegree[vs[idx3]] = 0; rdegree[vs[idx2]] = 0;
					assert(idx2 != idx3);
					if(idx2 < idx3) {
						vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
						vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					}
					else {
						vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
						vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
					}

					for(int i = 0;i < vs_size;i ++) if(vs[i] != v) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;

						if(!test_bit(t_matrix2, vs[i])) -- rd;

						if(!test_bit(t_matrix3, vs[i])) {
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								changes.pb(std::make_pair(v, vs[i]));
								++ rdegree[v];
							}
							else -- rd;
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
				}
				else if(total_connected == 2) {
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;
					//the following ensures that vs[idx2] is the vertex to be deleted
					if(!connected_3) {
						std::swap(vs[idx1], vs[idx3]);
						std::swap(t_matrix1, t_matrix3);
						rid[vs[idx1]] = idx1;
					}
					else if(!connected_2) {
						std::swap(vs[idx2], vs[idx3]);
						std::swap(t_matrix2, t_matrix3);
						rid[vs[idx2]] = idx2;
					}

					int v = vs[idx1], old_degree_v = rdegree[v];
					int w = vs[idx2], old_degree_w = rdegree[w];
					-- rdegree[v]; -- rdegree[w];
					contractions.pb(v); contractions.pb(w); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(2); //type-2 contraction

					rdegree[vs[idx3]] = 0;
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;

					for(int i = 0;i < vs_size;i ++) if(vs[i] != v&&vs[i] != w) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;

						if(!test_bit(t_matrix3, vs[i])) {
							-- rd;
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								changes.pb(std::make_pair(v, vs[i]));
								++ rdegree[v];
								++ rd;
							}

							if(test_bit(t_matrix2, vs[i])) {
								reverse_bit(t_matrix2, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, w);
								changes.pb(std::make_pair(w, vs[i]));
								++ rdegree[w];
								++ rd;
							}
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
					put_into_one_vector(current_clique_size, rid[w], vs, vs_size, rdegree[w] != old_degree_w, rdegree[w], rid);
				}
				else {
					assert(total_connected == 3);
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;

					int v1 = vs[idx1], old_degree_v1 = rdegree[v1];
					int v2 = vs[idx2], old_degree_v2 = rdegree[v2];
					int v3 = vs[idx3], old_degree_v3 = rdegree[v3];
					-- rdegree[v3];
					reverse_bit(t_matrix1, v2); reverse_bit(t_matrix2, v1);
					changes.pb(std::make_pair(v1, v2));
					contractions.pb(v1); contractions.pb(v2); contractions.pb(v3); contractions.pb(u); contractions.pb(3); //type-3 contraction

					for(int i = 0;i < vs_size;i ++) if(vs[i] != v1&&vs[i] != v2&&vs[i] != v3) {
						int &rd = rdegree[vs[i]];
						int old_rd = rd;

						char conn1 = test_bit(t_matrix1, vs[i]);
						char conn2 = test_bit(t_matrix2, vs[i]);
						char conn3 = test_bit(t_matrix3, vs[i]);

						if(conn1&&!conn2) {
							reverse_bit(t_matrix1, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v1);
							changes.pb(std::make_pair(v1, vs[i]));
							++ rd;
							++ rdegree[v1];
						}

						if(conn2&&!conn3) {
							reverse_bit(t_matrix2, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v2);
							changes.pb(std::make_pair(v2, vs[i]));
							++ rd;
							++ rdegree[v2];
						}

						if(conn3&&!conn1) {
							reverse_bit(t_matrix3, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v3);
							changes.pb(std::make_pair(v3, vs[i]));
							++ rd;
							++ rdegree[v3];
						}

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v1], vs, vs_size, rdegree[v1] != old_degree_v1, rdegree[v1], rid);
					put_into_one_vector(current_clique_size, rid[v2], vs, vs_size, rdegree[v2] != old_degree_v2, rdegree[v2], rid);
					put_into_one_vector(current_clique_size, rid[v3], vs, vs_size, rdegree[v3] != old_degree_v3, rdegree[v3], rid);
				}

#ifndef NDEBUG
				for(int i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
		}

#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) assert(rdegree[vs[i]] > 3&&vs_size - rdegree[vs[i]] + current_clique_size > max_clique.size());
#endif
	}

	void search_triangle_matrix_color(const int *vs, const int vs_size, const int *color, int &clique_size) {
		if(clique_size+2 == max_clique.size()) {
			int idx1 = 1;
			while(color[idx1] == color[0]) ++ idx1;
			assert(idx1+1 < vs_size);
			int idx2 = idx1+1;
			while(color[idx2] == color[idx1]) ++ idx2;
			assert(idx2 < vs_size);
			for(int i = 0;i < idx1&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
				assert(rid != nullptr);
				int *neighbors = rid;
				int neighbors_n = 0;
				for(int j = idx2;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
					neighbors[neighbors_n ++] = vs[j];
				}
				for(int j = idx1;j < idx2&&clique_size < max_clique.size();j ++) if(test_bit(t_matrix1, vs[j])) {
					unsigned char *t_matrix2 = matrix + vs[j]*matrix_len;
					for(int k = 0;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
						current_clique[clique_size ++] = vs[i];
						current_clique[clique_size ++] = vs[j];
						current_clique[clique_size ++] = neighbors[k];
						break;
					}
				}
			}
		}
		else if(clique_size+1 == max_clique.size()) {
			int idx1 = 1;
			while(color[idx1] == color[0]) ++ idx1;
			assert(idx1 < vs_size);
			for(int i = 0;i < idx1&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix = matrix + vs[i]*matrix_len;
				for(int j = idx1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = vs[j];
					break;
				}
			}
		}
		else if(clique_size == max_clique.size()) current_clique[clique_size ++] = vs[0];
	}

	void search_triangle_matrix(const int *vs, const int vs_size, int &clique_size) {
		if(clique_size+2 == max_clique.size()) {
			for(int i = 0;i < vs_size&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
				assert(rid != nullptr);
				int *neighbors = rid;
				int neighbors_n = 0;
				for(int j = i+1;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
					neighbors[neighbors_n ++] = vs[j];
				}
				for(int j = 0;j < neighbors_n&&clique_size < max_clique.size();j ++) {
					unsigned char *t_matrix2 = matrix + neighbors[j]*matrix_len;
					for(int k = j+1;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
						current_clique[clique_size ++] = vs[i];
						current_clique[clique_size ++] = neighbors[j];
						current_clique[clique_size ++] = neighbors[k];
						break;
					}
				}
			}
		}
		else if(clique_size+1 == max_clique.size()) {
			for(int i = 0;i < vs_size&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix = matrix + vs[i]*matrix_len;
				for(int j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = vs[j];
					break;
				}
			}
		}
		else if(clique_size == max_clique.size()) current_clique[clique_size ++] = vs[0];
	}

	void obtain_degrees(const int *vs, const int vs_size, int *degree) {
		for(int i = 0;i < vs_size;i ++) degree[vs[i]] = 0;
		for(int i = 0;i < vs_size;i ++) {
			int &d = degree[vs[i]];
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(int j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
				++ degree[vs[j]];
				++ d;
			}
		}
	}

	int split_vs(int *vs, const int vs_size, int *color, const int threshold) {
		for(int i = 0;i < threshold;i ++) head[i] = n;
		for(int i = vs_size;i > 0;i --) if(color[i-1] < threshold) {
			int j = i-1;
			next[vs[j]] = head[color[j]];
			head[color[j]] = vs[j];
		}
		int end_idx = 0;
		for(int i = 0;i < vs_size;i ++) if(color[i] >= threshold) {
			color[end_idx] = color[i];
			vs[end_idx] = vs[i];
			++ end_idx;
		}
		int new_size = end_idx;
		for(int i = threshold;i > 0;i --) for(int u = head[i-1];u != n;u = next[u]) {
			vs[new_size] = u;
			color[new_size ++] = i-1;
		}
		assert(new_size == vs_size);

		return end_idx;
	}

	void move_min_cardinality_color_to_front(int *vs, const int vs_size, int *color) {
		int c, min_cnt = n, cnt = 0;
		for(int i = 0;i < vs_size;i ++) {
			++ cnt;
			if(i+1 == vs_size||color[i] != color[i+1]) {
				if(cnt < min_cnt) {
					min_cnt = cnt;
					c = color[i];
				}
				cnt = 0;
			}
		}
		int idx = vs_size;
		while(idx > 0&&color[idx-1] != c) -- idx;
		for(int i = idx;i > 0;i --) if(color[i-1] != c) {
			std::swap(vs[i-1], vs[-- idx]);
			std::swap(color[i-1], color[idx]);
		}
	}

	//coloring a graph that is represented by matrix, aiming to minimize the number of vertices with color >= threshold
	//return the number of colors used
	int coloring_matrix_advanced(const int *vs, const int vs_size, int *color, const int start_color, const int threshold) {
		assert(rid != nullptr&&head != nullptr); // rid is used to temporarily store the color
		for(int i = 0;i < vs_size;i ++) head[i] = n;

		assert(start_color > 0&&start_color <= vs_size);
		for(int i = vs_size - start_color;i < vs_size;i ++) {
			int c = vs_size - i - 1;
			rid[vs[i]] = c;
			id[i] = vs[i];
			next[i] = head[c];
			head[c] = i;
		}

		for(int i = vs_size - start_color;i > 0;) {
			int u = vs[-- i];
			unsigned char *t_matrix = matrix + u*matrix_len;
			rid[u] = n;
			for(int j = 0;j < threshold;j ++) {
				char ok = 1;
				for(int k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					ok = 0;
					break;
				}
				if(ok) {
					rid[u] = j;
					id[i] = vs[i];
					next[i] = head[j];
					head[j] = i;
					break;
				}
			}
#ifndef _RECOLOR_
			continue;
#endif
			if(rid[u] < threshold) continue;

			for(int j = 0;j < threshold;j ++) {
				int cnt = 0, idx;
				for(int k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					++ cnt;
					if(cnt == 1) idx = k;
					else break;
				}
				assert(cnt > 0);
				if(cnt != 1) continue;
				unsigned char *tt_matrix = matrix + id[idx]*matrix_len;
				for(int ii = threshold;ii > 0;) {
					-- ii;
					if(ii == j) continue;

					char ok = 1;
					for(int k = head[ii];k != n;k = next[k]) if(test_bit(tt_matrix, id[k])) {
						ok = 0;
						break;
					}
					if(ok) {
						rid[id[idx]] = ii;
						id[i] = id[idx];
						next[i] = head[ii];
						head[ii] = i;

						rid[u] = j;
						id[idx] = vs[i];
						break;
					}
				}
				if(rid[u] < threshold) break;
			}
			if(rid[u] < threshold) continue;

			for(int j = threshold;;j ++) {
				char ok = 1;
				for(int k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					ok = 0;
					break;
				}
				if(ok) {
					rid[u] = j;
					id[i] = vs[i];
					next[i] = head[j];
					head[j] = i;
					break;
				}
			}
		}

		int max_color = 0;
		for(int i = 0;i < vs_size;i ++) {
			color[i] = rid[vs[i]];
			if(color[i] > max_color) max_color = color[i];
		}

#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(int j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) assert(color[i] != color[j]);
		}
#endif
		return max_color+1;
	}

	char kernelization_color(int &clique_size, int *vs, int &vs_size, int *color) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size-1;i ++) if(color[i] != color[i+1]) {
			for(int j = i+1;j < vs_size;j ++) assert(color[i] != color[j]);
		}
#endif
		while(true) {
			int idx = n;
			for(int i = 0;i < vs_size;i ++) {
				if((i==0||color[i]!=color[i-1])&&(i+1==vs_size||color[i]!=color[i+1])) {
					idx = i;
					break;
				}
			}
			if(idx == n) break;

			current_clique[clique_size ++] = vs[idx];

			int new_size = 0;
			unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
			for(int i = 0;i < vs_size;i ++) {
				if(i == idx) continue;
				if(test_bit(t_matrix, vs[i])) {
					vs[new_size] = vs[i];
					color[new_size ++] = color[i];
				}
				else if(i+1 == vs_size||color[i+1] != color[i]) {
					if(new_size == 0||color[i] != color[new_size-1]) return 1;
				}
			}
			vs_size = new_size;
		}

		if(clique_size > max_clique.size()) {
			store_a_larger_clique(clique_size, "kernel_color", 1);
			return 1;
		}
		return 0;
	}

	char reduce(int *vs, int &vs_size, int *color, const int threshold, int clique_size) {
#ifndef NDEBUG
		for(int i = 0;i < vs_size;i ++) {
			int d = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(int j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
			assert(degree[vs[i]] == d);
		}
#endif

		assert(rid != nullptr);
		int start = 0;
		for(;start < vs_size;start ++) if(degree[vs[start]] >= threshold) {
			int *neighbors = rid;
			int neighbors_n = 0;
			unsigned char *t_matrix = matrix + vs[start]*matrix_len;
			int color_cnt = 0;
			for(int i = start+1;i < vs_size;i ++) if(test_bit(t_matrix, vs[i])) {
				neighbors[neighbors_n ++] = i;
				if(!vis[color[i]]) {
					vis[color[i]] = 1;
					++ color_cnt;
				}
			}
			assert(neighbors_n == degree[vs[start]]);
			for(int i = 0;i < neighbors_n;i ++) {
				vis[color[neighbors[i]]] = 0;
				neighbors[i] = vs[neighbors[i]];
			}
			if(color_cnt < threshold) continue;

			if(degree[vs[start]] >= threshold + 1) break;

			char ok = 1;
			for(int i = 0;i < neighbors_n&&ok;i ++) {
				t_matrix = matrix + neighbors[i]*matrix_len;
				for(int j = i+1;j < neighbors_n;j ++) if(!test_bit(t_matrix, neighbors[j])) {
					ok = 0;
					break;
				}
			}
			if(ok) {
				current_clique[clique_size ++] = vs[start];
				for(int i = 0;i < neighbors_n;i ++) current_clique[clique_size ++] = neighbors[i];
				assert(clique_size == max_clique.size()+1);
				store_a_larger_clique(clique_size, "reduce", 1);
				return 1;
			}
		}

		if(start) {
			for(int i = start;i < vs_size;i ++) {
				vs[i-start] = vs[i];
				color[i-start] = color[i];
			}
			vs_size -= start;
		}

		if(vs_size == 0) return 1;
		return 0;
	}

	void recursive_search_clique_color_with_kernelization(const int level, int clique_size, int vs_begin, int vs_size, char kernel) {
		assert(clique_size <= max_clique.size()&&vs_buf.size() == color_buf.size()&&vs_size);
		int *vs = vs_buf.data() + vs_begin;
		int *color = color_buf.data() + vs_begin;

	#ifndef _KERNEL_
		printf("Wrong invocation of recursive_search_clique_color_with_kernelization!\n");
	#endif

	#ifdef _STATISTIC_
		++ branches;
		if(level > max_depth) max_depth = level;
	#endif

	#ifndef NDEBUG
		if(!kernel) {
			int color_cnt = 1;
			for(int i = 1;i < vs_size;i ++) if(color[i] != color[i-1]) {
				++ color_cnt;
				for(int j = i+1;j < vs_size;j ++) assert(color[j] != color[i-1]);
			}
			assert(color_cnt + clique_size == max_clique.size()+1);
		}
		for(int i = 0;i < vs_size;i ++) for(int j = 1;j < clique_size;j ++) assert(vs[i] != current_clique[j]);
		for(int i = 0;i < vs_size;i ++) for(int j = i+1;j < vs_size;j ++) assert(vs[i] != vs[j]);
		for(int i = 1;i < clique_size;i ++) for(int j = i+1;j < clique_size;j ++) assert(current_clique[i] == n || current_clique[i] != current_clique[j]);
	#endif

		if(clique_size + 3 > max_clique.size()) {
			if(kernel) search_triangle_matrix(vs, vs_size, clique_size);
			else search_triangle_matrix_color(vs, vs_size, color, clique_size);
			if(clique_size > max_clique.size()) store_a_larger_clique(clique_size, "search_triangle", 1);
			return ;
		}

		kernel = 1;

		int end_idx = 0, old_max_clique_size = max_clique.size();
		if(kernel) {
			obtain_degrees(vs, vs_size, degree);
			if(level) {
				int *rdegree = degree;
				for(int i = 0;i < vs_size;i ++) rdegree[vs[i]] = vs_size-1-degree[vs[i]];
				degree_one_two_three_reduction_with_folding_matrix(clique_size, vs, vs_size, rdegree, rid);
				for(int i = 0;i < vs_size;i ++) degree[vs[i]] = vs_size-1-rdegree[vs[i]];
				if(clique_size > max_clique.size()) store_a_larger_clique(clique_size, "degree_one_two_three_with_folding", 1);

				if(!vs_size||max_clique.size() > old_max_clique_size) return ;
			}

			//int cc_cnt = compute_connected_components(vs, vs_size);
			//if(cc_cnt > 1) ++ cc_larger_than_one[level];

			//degeneracy-based maximal clique
			int start_color = degeneracy_maximal_clique_matrix(clique_size, vs, vs_size, degree, 1, 1);

			if(max_clique.size() > old_max_clique_size) return ;

			int threshold = max_clique.size() - clique_size;
			//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
			int color_cnt = coloring_matrix_advanced(vs, vs_size, color, start_color, threshold);
			if(color_cnt <= threshold) return ;

			//reduce to (threshold+1)-core
			if(reduce(vs, vs_size, color, threshold, clique_size)) return ;

			//reorganize vs
	#ifdef _BRACH_ON_COLOR_
			end_idx = split_vs_on_color(vs, vs_size, color, threshold);
	#else
			end_idx = split_vs(vs, vs_size, color, threshold);
	#endif

			if(color_cnt <= threshold+1) kernel = 0;
		}
		else {
			move_min_cardinality_color_to_front(vs, vs_size, color);

			end_idx = 1;
			while(end_idx < vs_size&&color[end_idx] == color[0]) ++ end_idx;
		}

		while(vs_buf.size() < vs_begin+vs_size+vs_size) {
			vs_buf.pb(0);
			color_buf.pb(0);
		}
		vs = vs_buf.data() + vs_begin;
		color = color_buf.data() + vs_begin;

		for(int i = end_idx;i > 0&&max_clique.size() == old_max_clique_size;i --) {
			int *tvs = vs + vs_size;
			int *tcolor = color + vs_size;
			int tvs_end = 0;
			unsigned char *t_matrix = matrix + vs[i-1]*matrix_len;

			int upper_bound = 0;
			for(int j = (kernel?i:end_idx);j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
				assert(vs_buf.size() > vs_begin + vs_size + tvs_end);

				if(tvs_end == 0||color[j] != tcolor[tvs_end-1]) ++ upper_bound;

				tvs[tvs_end] = vs[j];
				tcolor[tvs_end ++] = color[j];
			}
			// Note that a vertex u may not connected to a vertex with color i for i < color(u)
			if(upper_bound + 1 + clique_size <= max_clique.size()) continue;

			current_clique[clique_size] = vs[i-1];
			int new_clique_size = clique_size + 1;
			int old_size = tvs_end;

			if(!kernel) {
				assert(upper_bound + clique_size == max_clique.size());
				if(kernelization_color(new_clique_size, tvs, tvs_end, tcolor)) continue;
			}

			int old_contraction_size = contractions.size(), old_changes_size = changes.size();
			recursive_search_clique_color_with_kernelization(level+1, new_clique_size, vs_begin + vs_size, tvs_end, kernel);
			while(contractions.size() > old_contraction_size) contractions.pop_back();
			while(changes.size() > old_changes_size) {
				std::pair<int, int> p = changes.back(); changes.pop_back();
				reverse_bit(matrix+p.first*matrix_len, p.second);
				reverse_bit(matrix+p.second*matrix_len, p.first);
			}
			if(!kernel&&old_size + end_idx == vs_size) break;

			vs = vs_buf.data() + vs_begin;
			color = color_buf.data() + vs_begin;
		}
	}

	void search_oriented(const int *peel_sequence, const int *core, const int *color, const int *local_UBs) {
	#ifndef NDEBUG
		for(int i = 0;i < n;i ++) assert(!vis[i]);
		assert(max_clique.size() >= 2&&max_clique.size() <= n);
	#endif

	#ifdef _STATISTIC_
		int matrix_cnt = 0;
		double total_density = 0;
		double min_density = 1;
		long total_kernel_effect = 0;

		int search_ego_cnt = 0;
		double total_ego_density = 0;
		double min_ego_density = 1;

		branches = 0;
		max_depth = 0;
	#endif

		int initial_size = max_clique.size();
		int *queue = new int[max_core];

		for(int i = n;i > 0;i --) {
			int u = peel_sequence[i-1];

			if(local_UBs[i-1] <= max_clique.size()) continue;
			if(core[u] < max_clique.size()) break;

			int old_max_clique_size = max_clique.size();

			//get N^+(u)
			int vs_size = 0;
			get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);
			assert(vs_size <= max_core);

			//color-based prune
			if(vs_size < max_clique.size()||color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;

			for(int j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

			//test max_clique \cup \{u\}
			if(max_clique.size() > initial_size&&greedy_extend(u, vis, 1)) {
				for(int j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
				continue;
			}

	#ifdef _KERNEL_
			//construct G[N^+(u)] and reduce by k-core
			int old_size = vs_size;
			construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
			kcore_reduction(vs_buf.data(), vs_size, vis, degree, max_clique.size()-1, queue);
			for(int j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
			if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;
	#endif

			//construct matrix
			int *rdegree = degree;
			construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

	#ifdef _STATISTIC_
			++ matrix_cnt;
			int total_edges = 0;
			for(int j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
			double density = double(total_edges)/(vs_size*(vs_size-1));
			if(density < min_density) min_density = density;
			total_density += density;
	#endif

			int current_clique_size = 1; current_clique[0] = u;
			contractions.clear(); changes.clear();

	#ifdef _KERNEL_
			old_size = vs_size;
			//degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
			degree_one_two_three_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, rid);
			if(current_clique_size > max_clique.size()) store_a_larger_clique(current_clique_size, "outside kernelization", 1);
			total_kernel_effect += old_size - vs_size;
			if(max_clique.size() > old_max_clique_size||!vs_size) continue;

	#ifdef _STATISTIC_
			++ search_ego_cnt;
			total_edges = 0;
			for(int j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
			density = double(total_edges)/(vs_size*(vs_size-1));
			if(density < min_ego_density) min_ego_density = density;
			total_ego_density += density;
	#endif

			changes.clear();
			recursive_search_clique_color_with_kernelization(0, current_clique_size, 0, vs_size, 1);
			//recursive_search_clique_color_without_kernelization(current_clique_size, 0, vs_size);
	#else
			recursive_search_clique_color_without_kernelization(0, current_clique_size, 0, vs_size);
	#endif
		}

		delete[] queue;

	#ifdef _STATISTIC_
		if(matrix_cnt == 0) {
			// printf("No matrix is constructed!\n");
		}
		else {
	#ifndef NDEBUG
			printf("Number of matrix constructed: %s\n", integer_to_string(matrix_cnt).c_str());
			printf("Average density: %.4lf, min density: %.4lf, average kernel_effect: %.4lf\n", total_density/matrix_cnt, min_density, double(total_kernel_effect)/matrix_cnt);

			printf("Number of egos searched: %s, branches: %s\n", integer_to_string(search_ego_cnt).c_str(), integer_to_string(branches).c_str());
	#endif
			if(search_ego_cnt == 0) search_ego_cnt = 1;
	#ifndef NDEBUG
			printf("Average ego_density: %.4lf, min ego_density: %.4lf\n", total_ego_density/search_ego_cnt, min_ego_density);
	#endif
		}
	#endif
	}

		// maximum_clique_color 求解最大团
	void maximum_clique_color(char exact, char initial_by_MC_EGO) {
		int *peel_sequence = new int[n];
		int *core = new int[n];
		int *color = new int[n];
		// std::vector<int> max_clique;

		max_clique.clear();
		char opt = 1;
		if(!initial_by_MC_EGO) opt = 0;
		int UB = degeneracy_maximal_clique_adjacency_list(peel_sequence, core, color, opt);
		assert(max_clique.size() >= 2);
		//printf("max_core: %u, UB: %u\n", max_core, UB);
		if(max_clique.size() < UB&&(exact||initial_by_MC_EGO)) {
			int old_size = max_clique.size();

			int *out_mapping = nullptr;
			shrink_graph(peel_sequence, core, color, out_mapping, pstart_o, edges_o);
			assert(n > 0);

			if(pend != nullptr) delete[] pend;
			pend = new ept[n];

			assert(current_clique == nullptr); current_clique = new int[UB];
			assert(mapping == nullptr); mapping = new int[max_core];
			assert(head == nullptr); head = new int[max_core];
			assert(next == nullptr); next = new int[max_core];
			assert(id == nullptr); id = new int[max_core];
			assert(matrix == nullptr);
	#ifdef _BITSET_
			matrix = new unsigned char[max_core*((max_core+7)/8)];
	#else
			matrix = new unsigned char[max_core*max_core];
	#endif

			assert(degree == nullptr); degree = new int[n];
			assert(rid == nullptr); rid = new int[n];
			assert(vis == nullptr); vis = new char[n];
			memset(vis, 0, sizeof(char)*n);

			int *local_UBs = new int[n];
			if(initial_by_MC_EGO) {
				UB = ego_degen(peel_sequence, core, color, local_UBs, UB);
	#ifndef NDEBUG
				//printf("\tMC-EGO Time: %s\n", integer_to_string(t.elapsed()).c_str());
				printf("There is a debug!\n");
	#endif
			}
			else {
				for(int i = 0;i < n;i ++) local_UBs[i] = n;
			}

			if(exact&&UB > max_clique.size()) search_oriented(peel_sequence, core, color, local_UBs);

			if(max_clique.size() > old_size) {
				for(int i = 0;i < max_clique.size();i ++) {
					assert(max_clique[i] < n);
					max_clique[i] = out_mapping[max_clique[i]];
				}
			}

			delete[] out_mapping;
			delete[] local_UBs;
		}

		delete[] color;
		delete[] core;
		delete[] peel_sequence;
	#ifndef NDEBUG
		if(exact||UB <= max_clique.size()) 
			//printf("\tMaximum Clique Size: %lu, Max Depth: %u, Total Time: %s\n", max_clique.size(), max_depth, integer_to_string(t.elapsed()).c_str());
			printf("There is an if print!\n");
		else 
			//printf("\tHeuristic Clique Size: %lu, UB: %u, Total Time: %s\n", max_clique.size(), UB, integer_to_string(t.elapsed()).c_str());
			printf("There is an else print!\n");
	#endif
	}
	// 找到一个图的最大团
	std::vector<int> max_clq_nodes()
	{
		maximum_clique_color(1,1);
		return max_clique;
	}

	void degen_order()
	{
		int total_edges = 0;
		for (int i = 0; i < n; i++)
		{
			peel_sequence[i] = i;
			total_edges += degree[i];
		}
		ListLinearHeap* heap = new ListLinearHeap(n, n - 1);
		heap->init(n, n - 1, peel_sequence, degree);
		int max_core = 0;
		int idx = n;
		for (int i = 0; i < n; i++)
		{
			int u, key;
			heap->pop_min(u, key);
			if (key > max_core)
				max_core = key;
			core[u] = max_core;
			peel_sequence[i] = u;
			for (int j = pstart[u]; j < pstart[u + 1]; j++){
				if (vis[edges[j].v] == 0){
						total_edges -= 2;
						heap->decrement(edges[j].v, 1);
				}
			}
			
		}
	}

	// 染色求一个上界，用来ban边。
	void color_ban()
	{
		// 先进行染色
		int mxcol = 0;
		if (col == nullptr)
			col = new int[n];
		if(vis == nullptr)
			vis = new char[n];
		degen_order();
		for(int i = 0; i < n; i++){
			std::cerr << peel_sequence[i] << " ";
		}
		std::cerr << "\n";
		for (int i = 0; i < n; ++i)
			col[i] = n;
		for (int i = 0; i < n; ++i)
		{
			int u = peel_sequence[n - i - 1];
			for (int j = pstart[u]; j < pstart[u + 1]; j++)
			{
				int c = col[edges[j].v];
				if (c != n)
					vis[c] = 1;
			}
			int tmp = mxcol;
			for (int j = 0; j < mxcol; j++){
				if (!vis[j])
				{
					tmp = j;
					break;
				}	
			}
			col[u] = tmp;
			if (tmp == mxcol)
				vis[mxcol++] = 0;
			for (int j = pstart[u]; j < pstart[u + 1]; j++)
			{
				int c = col[edges[j].v];
				if (c != n)
					vis[c] = 0;
			}
		}

		// 求每条边的公共邻居，查看公共邻居中的数目
		for(int u = 0; u < n; u++){
			//结点 u 的邻居的集合
			std::set<int> u_neighbor;
			// 找到邻居，并且判断是否ban
			for(int i = pstart[u]; i < pstart[u + 1]; i++){
				// 如果通向这个点的边没有被 ban,则标记为u的邻居
				if(!ban_edge[edges[i].tag])
					u_neighbor.insert(edges[i].v);
			}
			// 遍历边
			for(int j = pstart[u]; j < pstart[u + 1]; j++){
				// 记录所有的颜色
				std::set<int> cols;
				// 找到另外一个点
				std::vector<int> neighbor;
				int v = edges[j].v;
				int v_tag = edges[j].tag;
				for(int i = pstart[v]; i < pstart[v + 1]; i++){
					if(!ban_edge[edges[i].tag] && u_neighbor.find(edges[i].v) != u_neighbor.end()){
						cols.insert(col[edges[i].v]);
						neighbor.push_back(edges[i].v);
					}
				}
				if(2 + cols.size() <= LB){
					ban_edge[v_tag] = 1;
				}
				std::cerr << "edge:(" << u << ", " << v << "): " << 2 + cols.size() << "(" << neighbor.size() << ")" << "\n";
			}
		}
	}

// degeneracy-based k-defective clique
// return an upper bound of the maximum k-defective clique size
	void degen(){
		ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
		if (vis == nullptr)
        	vis = new char[n];
		int threshold = LB - 1;
		int queue_n = 0, new_size = 0;
		for (int i = 0; i < n; i++)                                                                                                    
			if (degree[i] < threshold)
				peel_sequence[queue_n++] = i;
		for (int i = 0; i < queue_n; i++)
		{
			int u = peel_sequence[i];
			degree[u] = 0;
			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (degree[edges[j].v] > 0)
				{
					if ((degree[edges[j].v]--) == threshold)
						peel_sequence[queue_n++] = edges[j].v;
				}
		}
		ept total_edges = 0;
		memset(vis, 0, sizeof(char) * n);
		for (int i = 0; i < n; i++){
			if (degree[i] >= threshold)
			{
				peel_sequence[queue_n + (new_size++)] = i;
				total_edges += degree[i];
			}
			else
			{
				vis[i] = 1;
				core[i] = 0;
			}
		}
		assert(queue_n + new_size == n);

		if (new_size != 0)
		{
			heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
			int max_core = 0;
			int idx = n;
			for (int i = 0; i < new_size; i++)
			{
				int u, key;
				heap->pop_min(u, key);
				if (key > max_core)
					max_core = key;
				core[u] = max_core;
				peel_sequence[queue_n + i] = u;
				for (ept j = pstart[u]; j < pstart[u + 1]; j++)
					if (vis[edges[j].v] == 0)
					{

						total_edges -= 2;
						heap->decrement(edges[j].v, 1);
					}
			}
		}
	}
	// 度数规约
	void core_shrink_graph()
	{
		int cnt = 0;
		out_mapping = new int[n];
		for (int i = 0; i < n; i++)
			out_mapping[i] = i;
		rid = new int[n];
		for (int i = 0; i < n; i++){
			//如果是不规约掉，那么就记录下来
			// std::cerr << "degree[" << i << "] = " << degree[i] << "\n";
			if (degree[i] > LB - 1){
				rid[i] = cnt;
				out_mapping[cnt] = out_mapping[i];
				++cnt;
			}
		}
		std::vector<int> tag_record;
		if (cnt != n)
		{
			cnt = 0;
			ept pos = 0;
			for (int i = 0; i < n; i++)
				if (degree[i] > LB - 1)
				{
					ept t_start = pstart[i];
					pstart[cnt] = pos;
					for (ept j = t_start; j < pstart[i + 1]; j++)
						if (degree[edges[j].v] > LB - 1)
						{
							edges[pos].tag = edges[j].tag;
							tag_record.push_back(edges[j].tag);
							edges[pos].v = rid[edges[j].v];
							pos++;
						}
					++cnt;
				}
			pstart[cnt] = pos;
			sort(tag_record.begin(),tag_record.end());
			int* tag_map = new int [m + 1];
			for(int i = 0; i < m + 1; i++)
				tag_map[i] = -1;
			for(int i = 0; i < tag_record.size(); i += 2){
				tag_map[tag_record[i]] = i / 2;
			}
			for(int i = 0; i < pos; i++){
				edges[i].tag = tag_map[edges[i].tag];
			}
			delete[] tag_map;
			for (int i = 0; i < cnt; i++)
				degree[i] = pstart[i + 1] - pstart[i];
			for (int i = 0; i < cnt; i++)
			{
				peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			}
			n = cnt;
			m = pos;
			// 更改边到点的映射
			int cnt = 0;
			for(int i = 0; i < n; i++){
				for(int j = pstart[i]; j < pstart[i + 1]; j++){
					edge_to_nodes[edges[j].tag].first = std::min(i, edges[j].v);
					edge_to_nodes[edges[j].tag].second = std::max(i, edges[j].v);
				}
			}
			// std::cerr << "num of edge tag: " << cnt << std::endl;
		}
		// 重新计算密度
		if(m != 0)
			R = (double)m /(double)((double)n * ((double)n - 1));
		else	
			R = 0;
	}
	// 用来ban边, 没有ban点
	bool ban()
	{
		if(m == 0){
			R = 0;
			return 0;
		}
		bool b = 0;
		ept pos = 0;
		std::vector<int> tag_record;
		for (int i = 0; i < n; i++)
		{
			ept t_start = pstart[i];
			pstart[i] = pos;
			for (ept j = t_start; j < pstart[i + 1]; j++){
				if (!ban_edge[edges[j].tag])
				{
					edges[pos].tag = edges[j].tag;
					tag_record.push_back(edges[j].tag);
					edges[pos].v = edges[j].v;
					pos++;
				}
				else{
					b = 1;
				}	
			}
		}
		pstart[n] = pos;
		sort(tag_record.begin(),tag_record.end());
		int* tag_map = new int [m + 1];
		for(int i = 0; i < m + 1; i++)
			tag_map[i] = -1;
		for(int i = 0; i < tag_record.size(); i += 2){
			tag_map[tag_record[i]] = i / 2;
		}
		// int max_trans = 0;
		for(int i = 0; i < pos; i++){
			edges[i].tag = tag_map[edges[i].tag];
		}
		for (int i = 0; i < n; i++)
			degree[i] = pstart[i + 1] - pstart[i];
		m = pos;
		// 更改边到点的映射
		int cnt = 0;

		for(int i = 0; i < n; i++){
			// std::cerr << pstart[i] << " " << pstart[i + 1] << std::endl;
			for(int j = pstart[i]; j < pstart[i + 1]; j++){
				edge_to_nodes[edges[j].tag].first = std::min(i, edges[j].v);
				edge_to_nodes[edges[j].tag].second = std::max(i, edges[j].v);
			}
		}
		
		for(int i = 0; i < m / 2; i++){
			ban_edge[i] = 0;
		}
		
		// 如果发生了ban边，把已加入团变量进行编号更改
		if(b){
			for(auto &c : clqs_edge){
				for(int i = 0; i < c.size(); i++){
					c[i] = tag_map[c[i]];
				}
			}
		}
		delete[] tag_map;
		if(m != 0)
			R = (double)m /(double)((double)n * ((double)n - 1));
		else	
			R = 0;
		return b;
	}

	// 读图
	void readDIMACS2Text(const char* filepath) {
		std::ifstream infile;
		char buf[1024];
		char tmps1[1024], tmps2[1024];
		//std::vector<std::pair<int, int> > epairs;
		std::vector<Undir_Edge_With_Tag> etuples; 
		std::vector<int> nodes;
		infile.open(filepath, std::ios::in);
		if (!infile.is_open()) {
			fprintf(stderr, "can not find file %s\n", filepath);
			exit(1);
		}
		int from, to;
		infile.getline(buf, 1024);
		//BUG:buff can be potentially empty string
		while (buf[0] != 'p') infile.getline(buf, 1024);
		sscanf(buf, "%s %s %u %lu", tmps1, tmps2, &n, &m);
		int lno = 1;
		while (infile.getline(buf, 1024)) {
			if (strlen(buf) == 0) continue;
			if (buf[0] != 'e') {
				fprintf(stderr, "ERROR in line %d\n", lno);
				continue;
			}
			sscanf(buf, "%s %d %d", tmps1, &from, &to);
			//std::cerr << from << " " << to << std::endl;
			if (from != to) {
				etuples.push_back((Undir_Edge_With_Tag){from, to, lno - 1});
				etuples.push_back((Undir_Edge_With_Tag){to, from, lno - 1});
				nodes.push_back(from);
				nodes.push_back(to);
			}
			lno++;
		}
		sort(nodes.begin(), nodes.end());
		nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

		sort(etuples.begin(), etuples.end(), cmp_with_node);
		int temp_record = etuples.size();
		etuples.erase(unique(etuples.begin(), etuples.end()), etuples.end());
		// 这里可能会有删除，根据边编号进行排序，排完序i/2改完编号再把顺序排回去，
		if(temp_record > etuples.size()){
			sort(etuples.begin(), etuples.end(), cmp_with_tag);
			for(int i = 0; i < etuples.size(); i++){
				// 修改两边
				etuples[i].tag = i / 2;
			}
			// 按照边的两顶点重新排序
			sort(etuples.begin(), etuples.end(), cmp_with_node);
		}

		int contn = 1;
		std::map<int, int> idmp;
		for (int i = 0; i < nodes.size(); i++) {
			idmp[nodes[i]] = i;
			if (nodes[i] != i) {
				contn = 0;
			}
		}
		if (contn == 0) printf("Node ids are not preserved! \n");

		n = nodes.size();
		m = etuples.size();
		printf("n = %s, (undirected) m = %s\n",
			integer_to_string(n).c_str(),
			integer_to_string(m / 2).c_str());

		pstart = new ept [n + 1];
		edges = new Edge_With_Tag [m];
		degree = new int [n];
		ban_edge = new int [m / 2];
		ban_node = new int [n];
		edge_to_nodes = new std::pair<int, int> [m / 2];
		core = new int [n];
		peel_sequence = new int [n];
		for(int i = 0; i < m / 2; i++)
			ban_edge[i] = 0;
		for(int i = 0; i < n; i++)
			ban_node[i] = 0;
		int j = 0;
		for (int i = 0; i < n; i++) {
			pstart[i] = j;
			while (j < m && etuples[j].u == nodes[i]) {
				edges[j].v = idmp[etuples[j].v];
				edges[j].tag = etuples[j].tag;
				++j;
			}
			std::sort(edges + pstart[i], edges + j, cmp_with_v);
			// 根据每个点边的个数统计点度数。
			degree[i] = j - pstart[i];
		}
		pstart[n] = j;
		// 根据度数对点排序
		degree_rank = new int[n];
		for(int i = 0;i < n; i++)
			degree_rank[i] = i;
		std::sort(degree_rank, degree_rank + n, [&](const int & a, const int & b){ return degree[a] > degree[b]; });

		// 安排从边到结点的映射
		for(int i = 0; i < n; i++){
			for(int j = pstart[i]; j < pstart[i+1]; j++){
				edge_to_nodes[edges[j].tag].first = std::min(i, edges[j].v);
				edge_to_nodes[edges[j].tag].second = std::max(i, edges[j].v);
			}
		}
		
		// 计算图的密度
		if(m != 0)
			R = (double)m /(double)((double)n * ((double)n - 1));
		else
			R = 0;
	}
	// 暂时只处理 dimacs-2 的图
	/*
	void read_binary_graph(const char* filepath) {
		FILE* f = open_file(filepath, "rb");
		if (f == NULL)
			perror("Error opening file");
		int tt;
		fread_wall(&tt, sizeof(int), 1, f);
		fread_wall(&n, sizeof(int), 1, f);
		fread_wall(&m, sizeof(int), 1, f);

		printf("n = %s; m = %s (undirected)\n", integer_to_string(n).c_str(), integer_to_string(m / 2).c_str());

		int* degree = new int[n];
		fread_wall(degree, sizeof(unsigned int), n, f);
		if (pstart == nullptr) pstart = new ept[n + 1];
		if (edges == nullptr) edges = new int[m];

		pstart[0] = 0;
		for (int i = 0; i < n; i++) {
			if (degree[i] > 0) {
				fread_wall(edges + pstart[i], sizeof(int), degree[i], f);

				// remove self loops and parallel edges
				int* buff = edges + pstart[i];
				std::sort(buff, buff + degree[i]);
				int idx = 0;
				for (int j = 0; j < degree[i]; j++) {
					if (buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
					if (buff[j] == i || (j > 0 && buff[j] == buff[j - 1])) continue;
					buff[idx++] = buff[j];
				}
				degree[i] = idx;
			}

			pstart[i + 1] = pstart[i] + degree[i];
		}

		fclose(f);
		delete[] degree;
	}

	int readRawDIM10Text(const char* filepath) {
		std::ifstream infile;
		const int SZBUF = 99999999;
		char* buf = new char[SZBUF];
		std::vector<std::pair<int, int> > epairs;
		std::vector<int> nodes;
		//FILE *f = Utility::open_file(filepath, "r");
		infile.open(filepath, std::ios::in);
		if (!infile.is_open()) {
			fprintf(stderr, "can not find file %s\n", filepath);
			exit(1);
		}

		infile.getline(buf, SZBUF);
		while (buf[0] == '%') infile.getline(buf, SZBUF);

		std::stringstream ss(buf);
		int fmt = 0;
		ss >> n >> m >> fmt;
		if (fmt != 0) {
			printf("Format of %s is not supported yet\n", filepath);
			exit(0);
		}
		m *= 2;
		pstart = new ept[n + 1];
		edges = new int[m];
		int j = 0;
		for (int u = 0; u < n; u++) {
			pstart[u] = j;
			infile.getline(buf, SZBUF);
			std::stringstream ss(buf);
			int nei;
			while (ss >> nei) {
				//printf("%d ", nei);
				if ((nei - 1) != u) {
					edges[j] = nei - 1;
					j++;
					//if (j==745)
					//	printf("pause\n");
				}
			}
			//printf("\n");
			std::sort(edges + pstart[u], edges + j);
		}
		pstart[n] = j;
		assert(j == m);
		printf("n:%u m:%lu\n", n, m / 2);
		return 0;
	}

	void readRawSNAPText(const char* filepath) {
		std::ifstream infile;
		char buf[1024];
		std::vector<std::pair<int, int> > epairs;
		std::vector<int> nodes;
		//FILE *f = Utility::open_file(filepath, "r");
		infile.open(filepath, std::ios::in);
		if (!infile.is_open()) {
			fprintf(stderr, "can not find file %s\n", filepath);
			exit(1);
		}
		int max_id = 0;
		int from, to;
		while (infile.getline(buf, 1024)) {
			char* p = buf;
			while (*p == ' ' && *p != '\0') p++;
			if (*p == '#' || *p == '\0') continue;
			std::stringstream ss(buf);
			ss >> from >> to;
			if (from != to) {
				epairs.push_back(std::make_pair(from, to));
				epairs.push_back(std::make_pair(to, from));
				nodes.push_back(from);
				nodes.push_back(to);
			}
		}
		infile.close();

		sort(nodes.begin(), nodes.end());
		nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

		sort(epairs.begin(), epairs.end());
		epairs.erase(unique(epairs.begin(), epairs.end()), epairs.end());

		int contn = 1;
		std::map<int, int> idmp;
		for (int i = 0; i < nodes.size(); i++) {
			idmp[nodes[i]] = i;
			if (nodes[i] != i) {
				contn = 0;
			}
		}
		if (contn == 0) printf("Node ids are not preserved! \n");

		n = nodes.size();
		m = epairs.size();
		printf("n = %s, (undirected) m = %s\n",
			integer_to_string(n).c_str(),
			integer_to_string(m / 2).c_str());

		pstart = new ept[n + 1];
		edges = new int[m];
		int j = 0;
		for (int i = 0; i < n; i++) {
			pstart[i] = j;
			while (j < m && epairs[j].first == nodes[i]) {
				edges[j] = idmp[epairs[j].second];
				++j;
			}
		}
		pstart[n] = j;
	}

	int writeBinaryGraph_ToBin(const char* filepath) {
		FILE* f = open_file(filepath, "wb");
		int tt = sizeof(int);
		fwrite(&tt, sizeof(int), 1, f); //length of int
		fwrite(&n, sizeof(int), 1, f);
		fwrite(&m, sizeof(int), 1, f);
		int* degree = new int[n];
		for (int i = 0; i < n; i++)
			degree[i] = pstart[i + 1] - pstart[i];
		fwrite(degree, sizeof(int), n, f);
		fwrite(edges, sizeof(int), m, f);
		fclose(f);
		return 0;
	}
	*/
};