#include "cliques.h"

int LB_with_fast_color(int method) {
    int n = Graph.n;
    
    // ================= 开始：5种取点排序方法 =================
    std::vector<int> order(n);

    switch(method) {
        case 1: // 1. 度数正序 (Degree Ascending)
            for (int i = 0; i < n; ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return Graph.degree[a] < Graph.degree[b];
            });
            break;
            
        case 2: // 2. 度数逆序 (Degree Descending) - 【已修复：直接复用原代码的排列】
            // 完美解决 std::sort 的不稳定性导致的 LB 质量下降问题
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.degree_rank[i];
            }
            break;
            
        case 3: // 3. degeneracy order 正序
            // 调用图中自带的方法生成 peel_sequence
            if(Graph.vis == nullptr)
			    Graph.vis = new char[Graph.n];
            Graph.degen_order();
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.peel_sequence[i];
            }
            break;
            
        case 4: // 4. degeneracy order 逆序
            if(Graph.vis == nullptr)
                Graph.vis = new char[Graph.n];
            Graph.degen_order();
            for (int i = 0; i < n; ++i) {
                // 将 peel_sequence 逆向存入
                order[i] = Graph.peel_sequence[n - 1 - i]; 
            }
            break;
            
        case 5: // 5. 随机排序 (Random)
            for (int i = 0; i < n; ++i) order[i] = i;
            std::random_shuffle(order.begin(), order.end());
            break;
            
        default: // 兜底：默认使用 degree_rank
            for (int i = 0; i < n; ++i) {
                order[i] = Graph.degree_rank[i];
            }
            break;
    }
    // ================= 结束：5种取点排序方法 =================

    int* used = new int[n];
    for (int i = 0; i < n; ++i)
        used[i] = 0;
    int* covered = new int[n];
    for (int i = 0; i < n; ++i)
        covered[i] = 0;
        
    int cnt = 0, mxclq_size = 0;
    std::vector<int> clqs_size;
    
    // 开始贪心寻找极大团
    for(int i = 0; i < n; i++){
        // 这里使用我们刚刚根据 method 生成的 order 数组
        int u = order[i]; 
        
        if(used[u])
            continue;
            
        // 存储可以加入极大团的邻居
        std::vector<int> vec;
        vec.clear();
        
        // 存储当前极大团
        std::vector<int> the_clq;
        the_clq.push_back(u);
        cnt++;
        used[u] = 1;
        
        // 遍历点 u 周围的所有点
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
                    // 选择u的邻居中连接其他顶点最多的一个点
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
    
    // 通过二分更新一个可行下界
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
    
    // 清理内存，防止内存泄漏
    delete[] used;
    delete[] covered;
    
    return l;
}

// int LB_with_fast_color(){
//     // 找包含各个点的极大团大小
//     int n = Graph.n;
//     int* used = new int[n];
//     for (int i = 0; i < n; ++i)
//         used[i] = 0;
//     int* covered = new int[n];
//     for (int i = 0; i < n; ++i)
//         covered[i] = 0;
//     int cnt = 0, mxclq_size = 0;
//     std::vector<int> clqs_size;
//     for(int i = 0; i < n; i++){

//         int u = Graph.degree_rank[i];
//         if(used[u])
//             continue;
//         // 存储可以加入极大团的邻居
//         std::vector<int> vec;
//         vec.clear();
//         // 存储当前极大团
//         std::vector<int> the_clq;
//         the_clq.push_back(u);
//         cnt++;
//         used[u] = 1;
//         // 遍历点 u 周围的所有点
//         for(int j = Graph.pstart[u]; j < Graph.pstart[u+1]; j++){
//             int to = Graph.edges[j].v;
//             if(used[to])
//                 continue;
//             covered[to] = cnt;
//             vec.push_back(to);
//         }
//         while(!vec.empty()){
//             // std::cerr << "The size of clq:" << the_clq.size() << " The size of vec: " << (int)vec.size() << std::endl;
//             // std::cerr << "LB:" << LB << std::endl;
//             if ((int)the_clq.size() + (int)vec.size() <= LB)
//                 break;
//             // std::cerr << "debug!\n";
//             int now_w = -1, now_d = -1;
//             for (int w : vec)
//             {
//                 int d = 0;
//                 for (int j = Graph.pstart[w]; j < Graph.pstart[w+1]; j++)
//                 {
//                     int to = Graph.edges[j].v;
//                     // 选择u的邻居中连接其他顶点最多的一个点
//                     if (!used[to] && covered[to] == cnt)
//                         ++d;
//                 }
//                 if (d >= now_d)
//                 {
//                     now_w = w;
//                     now_d = d;
//                 }
//             }
//             cnt++;
//             used[now_w] = 1;
//             the_clq.push_back(now_w);
//             vec.clear();
//             for (int j = Graph.pstart[now_w]; j < Graph.pstart[now_w+1]; j++)
//             {
//                 int to = Graph.edges[j].v;
//                 if (!used[to] && covered[to] == cnt - 1)
//                 {
//                     covered[to] = cnt;
//                     vec.push_back(to);
//                 }
//             }
//         }
//         // std::cerr << "The size of clq:" << the_clq.size() << std::endl;
//         if((int)the_clq.size() > LB){
//             clqs_size.push_back(the_clq.size());
//             mxclq_size = std::max((int)the_clq.size(), mxclq_size);
//         }
//         else{
//             for(auto x : the_clq){
//                 used[x] = 0;
//                 //covered[x] = 0;
//             }
//         }
//     }
//     std::sort(clqs_size.begin(), clqs_size.end(), std::greater<int>());
//     // 通过二分更新一个可行下界。
//     int l = 0, r = mxclq_size;
//     while (l <= r)
//     {
//         // std::cerr << "================================\n";
//         int mid = (l + r) / 2;
//         // std::cerr << "mid: " << mid <<std::endl;
//         int sum = 0;
//         for (int x : clqs_size)
//             if (x > mid){
//                 sum += func(mid + 1, x);
//                 //sum += x - mid;
//                 // std::cerr << x << ": " << func(mid + 1, x) << std::endl;
//             }
//         // std::cerr << "sum: " << sum << std::endl;
//         if (sum > K)
//             l = mid + 1;
//         else
//             r = mid - 1;
//     }

//     return l;
// }


void Max_clqs_of_edges(std::vector<int>& ban_edges){
    int* omega = new int [Graph.m / 2];
    memset(omega, 0, sizeof(int) * (Graph.m / 2));
    for(int i = 0; i < ban_edges.size(); i++){
        int edge_i = ban_edges[i];
        if(Graph.ban_edge[edge_i]){
            continue;
        }
        // clq_node 是包含这条边两个端点的团
        std::vector<int> clq_node = neg_clq_edges(Graph.edge_to_nodes[edge_i].first, Graph.edge_to_nodes[edge_i].second);
        clq_node.push_back(Graph.edge_to_nodes[edge_i].first);
        clq_node.push_back(Graph.edge_to_nodes[edge_i].second);

    #ifdef EXTRA
        std::sort(clq_node.begin(), clq_node.end());
    #endif

        // clq_edge 是由边组成的团
        std::vector<int> clq_edge;
        std::sort(clq_node.begin(),clq_node.end(),std::less<int>());
        for(int m = 0; m < clq_node.size(); m++){
            for(int k = m + 1; k < clq_node.size(); k++){
                int now_node = clq_node[m];
                int l = Graph.pstart[now_node];
                int r = Graph.pstart[now_node + 1] - 1;
                while(l < r){
                    int mid = (l + r) / 2;
                        if(Graph.edges[mid].v < clq_node[k])
                            l = mid + 1;
                        else    
                            r = mid;
                }
                if(Graph.edges[l].v == clq_node[k])
                    clq_edge.push_back(Graph.edges[l].tag);
            }
        }
        if(ceil(sqrt(2 * (int)clq_edge.size())) > LB){
            for(auto x : clq_edge)
                omega[x] = std::max(omega[x], (int)ceil(sqrt(2 * (int)clq_edge.size())));
            Graph.clqs_edge.push_back(clq_edge);

            // 将团的点存储起来在初始阶段
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
    // (两顶点公共邻居数目， 边号)

    //结点 u 的邻居的集合
    std::set<int> u_neighbor;
    // u 和 v的公共邻居
    std::vector<int> common_neighbor;
    for(int i = Graph.pstart[u]; i < Graph.pstart[u+1]; i++){
        // std::cerr << Graph.edges[i].v << " ";
        // 如果通向这个点的边没有被 ban,则标记为u的邻居
        if(!Graph.ban_edge[Graph.edges[i].tag])
            u_neighbor.insert(Graph.edges[i].v);
    }
    // std::cerr << std::endl;
    // 看v的邻居是否在u的邻居中
    for(int j = Graph.pstart[v]; j < Graph.pstart[v+1]; j++){
        // 如果这条边没有被ban，并且同时是u的邻居，那么把这个点加入子图。
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
    // 先进行染色
		int mxcol = 0;
		if (Graph.col == nullptr)
			Graph.col = new int[Graph.n];
		if(Graph.vis == nullptr)
			Graph.vis = new char[Graph.n];
		Graph.degen_order();
		// for(int i = 0; i < Graph.n; i++){
		// 	std::cerr << Graph.peel_sequence[i] << " ";
		// }
		// std::cerr << "\n";
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

		// 求每条边的公共邻居，查看公共邻居中的数目
		for(int u = 0; u < Graph.n; u++){
			//结点 u 的邻居的集合
			std::set<int> u_neighbor;
			// 找到邻居，并且判断是否ban
			for(int i = Graph.pstart[u]; i < Graph.pstart[u + 1]; i++){
				// 如果通向这个点的边没有被 ban,则标记为u的邻居
				if(!Graph.ban_edge[Graph.edges[i].tag])
					u_neighbor.insert(Graph.edges[i].v);
			}
			// 遍历边
			for(int j = Graph.pstart[u]; j < Graph.pstart[u + 1]; j++){
				// 记录所有的颜色
				std::set<int> cols;
				// 找到另外一个点
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
                // 加上u和v的颜色
                temp.push_back(std::make_pair((int)cols.size()+2, v_tag));
			}
		}
        std::sort(temp.begin(), temp.end());
        // std::cerr << temp.size() << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        return temp;
}

std::vector<int> neg_clq_edges(int u,int v){
    // std::cout << "two vertice are:" << "(" << u << ", " << v << ")" << std::endl;
    // 建立一个子图，存储包含边(u，v)的子图
    
    EdgeListGraph subgraph;
    std::vector<std::pair<int,int>> subedges;
    subedges.clear();
    int cnt = 0;
    //结点 u 的邻居的集合
    std::set<int> u_neighbor;
    // 子图中结点编号和原图中结点编号的映射
    std::vector<int> id_to_node; 
    // 原图顶点编号到子图顶点编号的映射
    int* dict = new int[Graph.n];
    for(int i = 0; i < Graph.n; i++)
        dict[i] = -1;
    // u 和 v的公共邻居
    std::vector<int> common_neighbor;
    // 找到点u和点v的所有 公共 邻居
    // 找到 u 的所有邻居，放入集合
    // std::cerr << "u is:" << u << "!\t";
    for(int i = Graph.pstart[u]; i < Graph.pstart[u+1]; i++){
        if(!Graph.ban_edge[Graph.edges[i].tag])
            u_neighbor.insert(Graph.edges[i].v);
    }
    // std::cerr << std::endl;
    // 看v的邻居是否在u的邻居中
    for(int j = Graph.pstart[v]; j < Graph.pstart[v+1]; j++){
        // 如果这条边没有被ban，并且同时是u的邻居，那么把这个点加入子图。
        if(!Graph.ban_edge[Graph.edges[j].tag] && u_neighbor.find(Graph.edges[j].v) != u_neighbor.end()){
            id_to_node.push_back(Graph.edges[j].v);
            dict[Graph.edges[j].v] = cnt++;
            common_neighbor.push_back(Graph.edges[j].v);
        }
    }
    // 如果最大可能的团大小不超过LB，说明答案不可能是这个值，直接消除影响
    if(cnt + 2 <= LB){
        delete[] dict;
        return std::vector<int>();      
    }
    // 把所有边加进去，通过边初始化子图
    for(auto x : common_neighbor){
        //assert(dict[x] >= 0);
        int min_node = dict[x];
        for(int i = Graph.pstart[x]; i < Graph.pstart[x+1]; i++){
            if(dict[Graph.edges[i].v] != -1 && min_node < dict[Graph.edges[i].v]){
                subedges.push_back(std::make_pair(min_node, dict[Graph.edges[i].v]));
            }
        }
    }
    // std::cerr << "subedges: " << subedges.size() << "\n";
    subgraph.sub_init(subedges);

    std::vector<int> subgraph_v;
    // if(subgraph.n <= 500){
    //     GET solr;
    //     subgraph_v = solr.get_mxclq(subgraph);
    // }
    if(subgraph.n <= 500){
        GET* solr = new GET();
        subgraph_v = solr->get_mxclq(subgraph);
        delete solr;  // 千万别忘了释放内存，防止内存泄漏
    }
    else{
        subgraph_v = subgraph.max_clq_nodes();
    }
    // 以点集的形式返回最大团
    std::vector<int> maxclq_node;
    // std::cerr << "id_to_node_size: " << id_to_node.size() << "\n";
    // std::cerr << "subgraph_v: " << subgraph_v.size() << "\n";
    // std::cerr << "common_neighbor: " << common_neighbor.size() << "\n";
    for(int x : subgraph_v){
        //std::cerr << x << ": " << id_to_node[x] << "\n";
        maxclq_node.push_back(id_to_node[x]);
    }
    delete[] dict;
    return maxclq_node;
}
std::vector<int> final_max_clique_nodes(){
    int cnt = 0;
    int* dict = new int[Graph.n];
    std::vector<int> id_to_node;
    // 初始化映射dict全部为-1
    for(int i = 0; i < Graph.n; i++)
        dict[i] = -1;
    // 如果所有边被ban了，那么这个点就无了，不然仍然存在。
    // 进行点的映射。
    for(int i = 0; i < Graph.n; i++){
        bool flag_edge = false;
        for(int j = Graph.pstart[i]; j < Graph.pstart[i+1]; j++){
            // 如果出现了一条没有被 ban 的边
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
    // std::cerr << "debug: 'final_max_clique_nodes' - The num of nodes: " << id_to_node.size() << std::endl;
    if (cnt == 0){
        return std::vector<int>();
    }
    // 初始化子图
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
    // std::cerr << "debug: 'final_max_clique_nodes' - The edges: ";
    // for(auto E : subedges){
    //     std::cerr << "(" << E.first << ", " << E.second << ")" << " ";
    // }
    // std::cerr << std::endl;

    for (int i = 0; i < Graph.n; i++){
        if(!Graph.ban_edge[i])
            dict[i] = -1;
    }
    // std::cerr << "debug: 'final_max_clique_nodes' - The num of edges in subgraph: " << subedges.size() << std::endl;
    subgraph.sub_init(subedges);
    std::vector<int> tmp;
    // 19年代码求解最大团
    // std::vector<int> tmp = subgraph.max_clq_nodes();
    // ILP求解最大团
    // std::vector<int> tmp = max_clq_ILP(subgraph);
    // Hua Jiang code
    // std::cerr << "!!!!!!!!!!!!!!!!!!!!!The num of subgraph: " << subgraph.n <<std::endl;
    // if(subgraph.n <= 500){
    //     GET solr;
    //     tmp = solr.get_mxclq(subgraph);
    // }
    if(subgraph.n <= 500){
        GET* solr = new GET();
        tmp = solr->get_mxclq(subgraph);
        delete solr;  // 千万别忘了释放内存，防止内存泄漏
    }
    else{
        tmp = subgraph.max_clq_nodes();
    }
    
    std::vector<int> fin;
    for (int x : tmp)
        fin.push_back(id_to_node[x]);
    // std::cerr << "debug: 'final_max_clique_nodes' - The size of fin: " << fin.size() << std::endl;
    return fin;
}