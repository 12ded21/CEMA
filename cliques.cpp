#include "cliques.h"

namespace {
int find_edge_tag(int u, int v) {
    if (u < 0 || v < 0 || u >= Graph.n || v >= Graph.n || u == v) {
        return -1;
    }
    if (u > v) {
        std::swap(u, v);
    }

    int l = Graph.pstart[u];
    int r = Graph.pstart[u + 1];
    while (l < r) {
        int mid = l + (r - l) / 2;
        if (Graph.edges[mid].v < v) {
            l = mid + 1;
        } else {
            r = mid;
        }
    }

    if (l < Graph.pstart[u + 1] && Graph.edges[l].v == v) {
        return Graph.edges[l].tag;
    }
    return -1;
}
}

int LB_with_fast_color(int method) {
    int n = Graph.n;
    
    std::vector<int> order(n);

    switch(method) {
        case 1:
            for (int i = 0; i < n; ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return Graph.degree[a] < Graph.degree[b];
            });
            break;
            
        case 2:
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.degree_rank[i];
            }
            break;
            
        case 3:
            if(Graph.vis == nullptr) {
                Graph.vis = new char[Graph.n];
            }
            memset(Graph.vis, 0, sizeof(char) * Graph.n);
            Graph.degen_order();
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.peel_sequence[i];
            }
            break;
            
        case 4:
            if(Graph.vis == nullptr)
                Graph.vis = new char[Graph.n];
            Graph.degen_order();
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.peel_sequence[n - 1 - i]; 
            }
            break;
            
        case 5:
            for (int i = 0; i < n; ++i) order[i] = i;
            std::random_shuffle(order.begin(), order.end());
            break;
            
        default:
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.degree_rank[i];
            }
            break;
    }

    int* used = new int[n];
    for (int i = 0; i < n; ++i)
        used[i] = 0;
    int* covered = new int[n];
    for (int i = 0; i < n; ++i)
        covered[i] = 0;
        
    int cnt = 0, mxclq_size = 0;
    std::vector<int> clqs_size;
    
    for(int i = 0; i < n; i++){
        int u = order[i]; 
        
        if(used[u])
            continue;
            
        std::vector<int> vec;
        vec.clear();
        
        std::vector<int> the_clq;
        the_clq.push_back(u);
        cnt++;
        used[u] = 1;
        
        for(int j = Graph.pstart[u]; j < Graph.pstart[u+1]; j++){
            int to = Graph.edges[j].v;
            if(used[to])
                continue;
            covered[to] = cnt;
            vec.push_back(to);
        }
        
        while(!vec.empty()){
            if ((int)the_clq.size() + (int)vec.size() <= LB)
                break;
                
            int now_w = -1, now_d = -1;
            for (int w : vec)
            {
                int d = 0;
                for (int j = Graph.pstart[w]; j < Graph.pstart[w+1]; j++)
                {
                    int to = Graph.edges[j].v;
                    if (!used[to] && covered[to] == cnt)
                        ++d;
                }
                if (d >= now_d)
                {
                    now_w = w;
                    now_d = d;
                }
            }
            
            cnt++;
            used[now_w] = 1;
            the_clq.push_back(now_w);
            vec.clear();
            
            for (int j = Graph.pstart[now_w]; j < Graph.pstart[now_w+1]; j++)
            {
                int to = Graph.edges[j].v;
                if (!used[to] && covered[to] == cnt - 1)
                {
                    covered[to] = cnt;
                    vec.push_back(to);
                }
            }
        }
        
        if((int)the_clq.size() > LB){
            clqs_size.push_back(the_clq.size());
            mxclq_size = std::max((int)the_clq.size(), mxclq_size);
        }
        else{
            for(auto x : the_clq){
                used[x] = 0;
            }
        }
    }
    
    std::sort(clqs_size.begin(), clqs_size.end(), std::greater<int>());
    
    int l = 0, r = mxclq_size;
    while (l <= r)
    {
        int mid = (l + r) / 2;
        int sum = 0;
        for (int x : clqs_size) {
            if (x > mid){
                sum += func(mid + 1, x);
            }
        }
        if (sum > K)
            l = mid + 1;
        else
            r = mid - 1;
    }
    
    delete[] used;
    delete[] covered;
    
    return l;
}

void Max_clqs_of_edges(std::vector<int>& ban_edges){
    int* omega = new int [Graph.m / 2];
    memset(omega, 0, sizeof(int) * (Graph.m / 2));
    for(int i = 0; i < ban_edges.size(); i++){
        int edge_i = ban_edges[i];
        if(Graph.ban_edge[edge_i]){
            continue;
        }
        std::vector<int> clq_node = neg_clq_edges(Graph.edge_to_nodes[edge_i].first, Graph.edge_to_nodes[edge_i].second);
        clq_node.push_back(Graph.edge_to_nodes[edge_i].first);
        clq_node.push_back(Graph.edge_to_nodes[edge_i].second);

    #ifdef EXTRA
        std::sort(clq_node.begin(), clq_node.end());
    #endif

        std::vector<int> clq_edge;
        std::sort(clq_node.begin(),clq_node.end(),std::less<int>());
        for(int m = 0; m < clq_node.size(); m++){
            for(int k = m + 1; k < clq_node.size(); k++){
                int edge_tag = find_edge_tag(clq_node[m], clq_node[k]);
                if(edge_tag >= 0)
                    clq_edge.push_back(edge_tag);
            }
        }
        if(ceil(sqrt(2 * (int)clq_edge.size())) > LB){
            for(auto x : clq_edge)
                omega[x] = std::max(omega[x], (int)ceil(sqrt(2 * (int)clq_edge.size())));
            Graph.clqs_edge.push_back(clq_edge);
        #ifdef EXTRA
            Graph.clqs_node.push_back(clq_node);
        #endif

        }
        else{
            omega[edge_i] =(int)ceil(sqrt(2 * (int)clq_edge.size()));
        #ifdef CLIQUE
            Graph.ban_edge[edge_i] = 1;
        #endif
        }
    }
}

int common_neighbor_num(int u, int v){
    std::set<int> u_neighbor;
    std::vector<int> common_neighbor;
    for(int i = Graph.pstart[u]; i < Graph.pstart[u+1]; i++){
        if(!Graph.ban_edge[Graph.edges[i].tag])
            u_neighbor.insert(Graph.edges[i].v);
    }
    for(int j = Graph.pstart[v]; j < Graph.pstart[v+1]; j++){
        if(!Graph.ban_edge[Graph.edges[j].tag] && u_neighbor.find(Graph.edges[j].v) != u_neighbor.end()){
            common_neighbor.push_back(Graph.edges[j].v);
        }
    }
    return (int)common_neighbor.size();
}


std::vector<std::pair<int, int>> common_neighbor_vector(){
    std::vector<std::pair<int, int>> temp;
    for(int i = 0; i < Graph.n; i++){
        for(int j = Graph.pstart[i]; j < Graph.pstart[i+1]; j++){
            if(!Graph.ban_edge[Graph.edges[j].tag] && i < Graph.edges[j].v){
                int neighbor_num = common_neighbor_num(i, Graph.edges[j].v);
                temp.push_back(std::make_pair(neighbor_num, Graph.edges[j].tag));
            }
        }
    }
    std::sort(temp.begin(), temp.end());
    return temp;
}
std::vector<std::pair<int, int>> color_vector(){
    if(Graph.m == 0){
        return std::vector<std::pair<int, int>>{(1, 0)};
    }
    std::vector<std::pair<int, int>> temp;
		int mxcol = 0;
		if (Graph.col == nullptr)
			Graph.col = new int[Graph.n];
		if(Graph.vis == nullptr)
			Graph.vis = new char[Graph.n];
		Graph.degen_order();
		for (int i = 0; i < Graph.n; ++i)
			Graph.col[i] = Graph.n;
		for (int i = 0; i < Graph.n; ++i)
		{
			int u = Graph.peel_sequence[Graph.n - i - 1];
			for (int j = Graph.pstart[u]; j < Graph.pstart[u + 1]; j++)
			{
				int c = Graph.col[Graph.edges[j].v];
				if (c != Graph.n)
					Graph.vis[c] = 1;
			}
			int tmp = mxcol;
			for (int j = 0; j < mxcol; j++){
				if (!Graph.vis[j])
				{
					tmp = j;
					break;
				}	
			}
			Graph.col[u] = tmp;
			if (tmp == mxcol)
				Graph.vis[mxcol++] = 0;
			for (int j = Graph.pstart[u]; j < Graph.pstart[u + 1]; j++)
			{
				int c = Graph.col[Graph.edges[j].v];
				if (c != Graph.n)
					Graph.vis[c] = 0;
			}
		}

		for(int u = 0; u < Graph.n; u++){
			std::set<int> u_neighbor;
			for(int i = Graph.pstart[u]; i < Graph.pstart[u + 1]; i++){
				if(!Graph.ban_edge[Graph.edges[i].tag])
					u_neighbor.insert(Graph.edges[i].v);
			}
			for(int j = Graph.pstart[u]; j < Graph.pstart[u + 1]; j++){
				std::set<int> cols;
				std::vector<int> neighbor;
				int v = Graph.edges[j].v;
				int v_tag = Graph.edges[j].tag;
                if(u >= v)
                    continue;
				for(int i = Graph.pstart[v]; i < Graph.pstart[v + 1]; i++){
					if(!Graph.ban_edge[Graph.edges[i].tag] && u_neighbor.find(Graph.edges[i].v) != u_neighbor.end()){
						cols.insert(Graph.col[Graph.edges[i].v]);
					}
				}
                temp.push_back(std::make_pair((int)cols.size()+2, v_tag));
			}
		}
        std::sort(temp.begin(), temp.end());
        return temp;
}

std::vector<int> neg_clq_edges(int u,int v){
    EdgeListGraph subgraph;
    std::vector<std::pair<int,int>> subedges;
    subedges.clear();
    int cnt = 0;
    std::set<int> u_neighbor;
    std::vector<int> id_to_node; 
    int* dict = new int[Graph.n];
    for(int i = 0; i < Graph.n; i++)
        dict[i] = -1;
    std::vector<int> common_neighbor;
    for(int i = Graph.pstart[u]; i < Graph.pstart[u+1]; i++){
        if(!Graph.ban_edge[Graph.edges[i].tag])
            u_neighbor.insert(Graph.edges[i].v);
    }
    for(int j = Graph.pstart[v]; j < Graph.pstart[v+1]; j++){
        if(!Graph.ban_edge[Graph.edges[j].tag] && u_neighbor.find(Graph.edges[j].v) != u_neighbor.end()){
            id_to_node.push_back(Graph.edges[j].v);
            dict[Graph.edges[j].v] = cnt++;
            common_neighbor.push_back(Graph.edges[j].v);
        }
    }
    if(cnt + 2 <= LB){
        delete[] dict;
        return std::vector<int>();      
    }
    for(auto x : common_neighbor){
        assert(dict[x] >= 0);
        int min_node = dict[x];
        for(int i = Graph.pstart[x]; i < Graph.pstart[x+1]; i++){
            if(dict[Graph.edges[i].v] != -1 && min_node < dict[Graph.edges[i].v]){
                subedges.push_back(std::make_pair(min_node, dict[Graph.edges[i].v]));
            }
        }
    }
    if(common_neighbor.empty()){
        delete[] dict;
        return std::vector<int>();
    }
    if(subedges.empty()){
        delete[] dict;
        return std::vector<int>{common_neighbor[0]};
    }

    subgraph.sub_init(subedges, cnt);

    std::vector<int> subgraph_v;
    if(subgraph.n <= 500){
        std::cout << subgraph.n << " " << subgraph.m << std::endl;
        GET* solr = new GET();
        subgraph_v = solr->get_mxclq(subgraph);
        delete solr;
    }
    else{
        subgraph_v = subgraph.max_clq_nodes();
    }
    std::vector<int> maxclq_node;
    for(int x : subgraph_v){
        if(x >= 0 && x < (int)id_to_node.size())
            maxclq_node.push_back(id_to_node[x]);
    }
    delete[] dict;
    return maxclq_node;
}
std::vector<int> final_max_clique_nodes(){
    int cnt = 0;
    int* dict = new int[Graph.n];
    std::vector<int> id_to_node;
    for(int i = 0; i < Graph.n; i++)
        dict[i] = -1;
    for(int i = 0; i < Graph.n; i++){
        bool flag_edge = false;
        for(int j = Graph.pstart[i]; j < Graph.pstart[i+1]; j++){
            if(!Graph.ban_edge[Graph.edges[j].tag]){
                flag_edge = true;
                break;
            }
        }
        if(!flag_edge){
            continue;
        }
        else{
            id_to_node.push_back(i);
            dict[i] = cnt++;
        }
    }
    if (cnt == 0){
        delete[] dict;
        return std::vector<int>();
    }
    EdgeListGraph subgraph;
    std::vector<std::pair<int, int>> subedges;
    subedges.clear();
    for(int i = 0; i < Graph.n; i++){
        int min_node = dict[i];
        for(int j = Graph.pstart[i]; j < Graph.pstart[i+1]; j++){
            if(!Graph.ban_edge[Graph.edges[j].tag] && min_node < dict[Graph.edges[j].v]){
                subedges.push_back(std::make_pair(min_node, dict[Graph.edges[j].v]));
            }
        }
    }

    for (int i = 0; i < Graph.n; i++){
        if(!Graph.ban_edge[i])
            dict[i] = -1;
    }
    subgraph.sub_init(subedges, cnt);
    std::vector<int> tmp;
    if(subgraph.n <= 500){
        std::cout << subgraph.n << " " << subgraph.m << std::endl;
        GET* solr = new GET();
        tmp = solr->get_mxclq(subgraph);
        delete solr;
    }
    else{
        tmp = subgraph.max_clq_nodes();
    }
    
    std::vector<int> fin;
    for (int x : tmp){
        if(x >= 0 && x < (int)id_to_node.size())
            fin.push_back(id_to_node[x]);
    }
    delete[] dict;
    return fin;
}
