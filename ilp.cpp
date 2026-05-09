#include "ilp.h"

// ILP的成员函数

// 定义构造函数
ILP::ILP(int _m, int k, int z_num, int _LB, int _UB){
    // 其中 z_num 为 z 需要的大小。
    m = _m;
    y = _UB - 1;
    k = std::min(k, _m);
    model = IloModel(env);
    x_e = IloBoolVarArray(env, m);
    z = IloBoolVarArray(env, z_num);
    z_cnt = 0;
    s = IloIntVar(env, 0, k);
    // 优化目标是最小化s
    model.add(IloMinimize(env, s));
    // 添加表达式约束
    model.add(IloSum(x_e) == s);
    // 定义求解器
    cplex = IloCplex(model);
    // 定义参数
    cplex.setOut(env.getNullStream());
    cplex.setError(env.getNullStream());
    cplex.setParam(IloCplex::Param::MIP::Display, 2);
    cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);

    cplex.setParam(IloCplex::Param::Threads, 1);
    // cplex.setParam(IloCplex::Param::TimeLimit, 300);
    // 优化
    
    // 注册懒惰回调函数
    cplex.use(new (env) LazyConstraintCallback(env, x_e, z, s, y));
}

// 定义析构函数
ILP::~ILP(){
    env.end();
}

// 懒惰约束回调类的构造函数
ILP::LazyConstraintCallback::LazyConstraintCallback(IloEnv env, IloBoolVarArray x_e, IloBoolVarArray z, IloIntVar s, int& y)
            : IloCplex::LazyConstraintCallbackI(env), _x_e(x_e), _z(z), _s(s), _env(env), _y(y){}

// 懒惰约束回调函数
void ILP::LazyConstraintCallback::main(){
    std::cerr << "================================== Callback function begin ==================================\n";
    // m表示结点变量的数量
    int m = _x_e.getSize();
    // LB = max(LB, (int)getValue(_lb));
    // v 中存储删去的边 e 的编号。
    std::vector<int> del_edges;
    // 与ban边连接的顶点,这些顶点可能是找到的最大团里的顶点。
    // std::vector<int> ban_node;
    for (int i = 0; i < m; ++i){
        if (round(getValue(_x_e[i])) == 1){
            del_edges.push_back(i);
        }
    }
    // std::cerr << "\n";

    // 去重
    // std::sort(ban_node.begin(), ban_node.end());
    // ban_node.erase(std::unique(ban_node.begin(), ban_node.end()), ban_node.end());

    // 设置ban边，用于处理最大团
    for (int x : del_edges)
        Graph.ban_edge[x] = 1;

    // 在剩余图中找到最大团
    clock_t st = clock();
    std::vector<int> clq;
    if(Graph.m / 2 - (int)del_edges.size() > 0)
        clq = final_max_clique_nodes();
    else
        clq.push_back(0);
    clock_t en = clock();
    double duration = (en - st) * 1.0 / CLOCKS_PER_SEC;
    std::cerr << "callback: Time of get_maxclq: " << duration << "\n";
    std::cerr << "callback: size of the max_clq: " << clq.size() << "\n";
    // 最大团中的边
    std::vector<int> mxclq_edges;
    // ans = min(ans, (int)clq.size());
    std::vector<int> neis;
    std::vector<std::vector<int>> clqs;

    // 把顶点搞成边
    std::sort(clq.begin(),clq.end(),std::less<int>());
    for(int m = 0; m < clq.size(); m++){
        for(int k = m + 1; k < clq.size(); k++){
            int now_node = clq[m];
            int l = Graph.pstart[now_node];
            int r = Graph.pstart[now_node+1] - 1;
            // bool search_flag = 0;
            while(l < r){
                int mid = (l + r) / 2;
                    if(Graph.edges[mid].v < clq[k])
                        l = mid + 1;
                    else    
                        r = mid;
            }
            if(Graph.edges[l].v == clq[k])
                mxclq_edges.push_back(Graph.edges[l].tag);
        }
    }

    int temp = _y;
    // while((int)clq.size() <= _y)
    //     _y--;
    if((int)clq.size() <= _y)
        _y = (int)clq.size() - 1;
    now_ans = _y + 1;
    if(_y <= 0 || _y < LB){
        // 如果更新了，那么换上。
        Graph.global_delete_edges = del_edges;
        return ;
    }
        
    if (temp > _y)
    {
        // 如果更新了，那么换上。
        Graph.global_delete_edges = del_edges;

        // 重新加入之前加入的所有约束
        // 加入边约束，对于每个约束，更改右侧不等式的func函数值
        for (auto clq_ : Graph.clqs_edge){
            IloExpr lazyConstraint(_env);
            for(auto edges_ : clq_){
                lazyConstraint += _x_e[edges_];
            }
            add(lazyConstraint >= func(_y+1, ceil(sqrt(2 * clq_.size()))));
            lazyConstraint.end();
            all_constraints++;
        }
        // 加入点约束，对于每个约束，更改右侧不等式的y值
        std::cerr << "The clqs_z: " << Graph.clqs_z.size() << "\n";

        for(auto clq_z : Graph.clqs_z){
            IloExpr lazyConstraint(_env);
            for(auto z_tag_ : clq_z){
                lazyConstraint += _z[z_tag_];
            }
            add(lazyConstraint <= _y);
            lazyConstraint.end();
            // 增加一个点约束
            all_constraints++;
        }
    }
    std::vector<int> edges_clq;
    // 各个团需要添加的边
    std::vector<std::vector<int>> edges_clqs;

    // 拓展点
    if (clq.size() > _y){
        // 添加一个线性表达式
        // 把整张图的最大团加入。
        IloExpr lazyConstraint(_env);
        // std::cerr << "The added testriction: ";
        for (int i : mxclq_edges)
        {
            // std::cerr << i << "\t";
            lazyConstraint += _x_e[i];
        }
        // 添加懒惰约束，最大团
        add(lazyConstraint >= func(_y+1, clq.size()));
        lazyConstraint.end();
        Graph.clqs_edge.push_back(mxclq_edges);
        // 增加一个边约束
        all_constraints++;
    }
    // 检查现在的约束个数
    std::cerr << "The num of constraints: " << getNrows() << std::endl;


    //复原
    for (int x : del_edges)
        Graph.ban_edge[x] = 0;
    std::cerr << "callback: The s is: " << round(getValue(_s)) << std::endl;
    std::cerr << "callback: The y is: " << _y << std::endl;
    // std::cerr << "debug!\n" ;
}

// 添加单个团
void ILP::add_edge_clq(std::vector<int> &clq){
    IloExpr sum_cut_in_clq(env);
    for(int i : clq){
        sum_cut_in_clq += x_e[i];
    }
    model.add(sum_cut_in_clq >= func(y+1, ceil(sqrt(2 * clq.size()))));
    sum_cut_in_clq.end();
    all_constraints++;
}

// 批量添加团
void ILP::add_edge_clqs(std::vector<std::vector<int>> &clqs)
{
    for (int i = 0; i < clqs.size(); ++i)
        add_edge_clq(clqs[i]);
}

// 以点的形式添加额外约束
// 添加单个团, 增加ub变量代表
void ILP::add_node_clq(std::vector<int> &clq){
#ifdef EXTRA
    std::vector<int> z_vec;

    IloExpr sum_z(env);
    for(int i = 0; i < clq.size(); i++){
        int u = clq[i];
        // 遍历该团中所有与 u 相邻的边，二分找到边的另外一点大于u的边并记录下来
        // 在某一序列顺序下，仅添加从前向后的边，表示有向边。
        // std::cout << u << "(" << Graph.degree[u] << ")" << " ";
        IloExpr z_plus_sum_edge(env);
        for(int j = i + 1; j < clq.size(); j++){
            int v = clq[j];
            int l = Graph.pstart[u];
            int r = Graph.pstart[u + 1] - 1;
            while(l < r){
                int mid = (l + r) / 2;
                    if(Graph.edges[mid].v < v)
                        l = mid + 1;
                    else    
                        r = mid;
            }

            // 找到边之后将每条边编号所对应的变量加起来，构造约束。

            if(Graph.edges[l].v == v)
                z_plus_sum_edge += x_e[Graph.edges[l].tag];
            else
                std::cerr << "There is an error when get the node clq\n";  
        }
        z_plus_sum_edge += z[z_cnt];
        model.add(z_plus_sum_edge >= 1);
        z_plus_sum_edge.end();
        // 增加点约束
        all_constraints++;
        sum_z += z[z_cnt];
        z_vec.push_back(z_cnt);
        z_cnt++;
    }
    // std::exit(0);
    model.add(sum_z <= y);
    sum_z.end();
    Graph.clqs_z.push_back(z_vec);
    // 增加点约束
    all_constraints++;
#endif

}

// 批量添加团
void ILP::add_node_clqs(std::vector<std::vector<int>> &clqs)
{
    #ifdef EXTRA
    for (int i = 0; i < clqs.size(); ++i){
        /*
            将clq修改一下，先度数降序排列一个添加约束。
        */
        std::sort(clqs[i].begin(), clqs[i].end(), [](int a, int b){
            return Graph.pstart[a+1] - Graph.pstart[a] > Graph.pstart[b+1] - Graph.pstart[b];
        });

        add_node_clq(clqs[i]);

        /*
            将clq修改一下，后度数升序排列一个添加约束。
        */
        std::sort(clqs[i].begin(), clqs[i].end(), [](int a, int b){
            return Graph.pstart[a+1] - Graph.pstart[a] < Graph.pstart[b+1] - Graph.pstart[b];
        });

        add_node_clq(clqs[i]);
    }
    #endif
}

// 求解器求解
int ILP::solve(){
    int solved = 0;
    cplex.solve();
    
    IloAlgorithm::Status status = cplex.getStatus();
    ilpnode = cplex.getNnodes();
    ilpcut = cplex.getNcuts(IloCplex::CutUser);
    // 断言只能是两种情况，超时或者无解(在试验阶段需要改成其他的)
    // assert(status == IloAlgorithm::Infeasible || 
    //         status == IloAlgorithm::Unknown && cplex.getCplexStatus() == IloCplex::AbortTimeLim);
    
    // 打印日志
    /*
    if(status == IloAlgorithm::Infeasible){
        printf("The CEMA instance is solved successfully!\n");
        solved = 1;
    }
    else{
        IloCplex::Status detailed_status = cplex.getCplexStatus();
        if(detailed_status == IloCplex::AbortTimeLim)
            printf("The time limit is reached!\n");
        solved = 0;
    }
    */

    if(status == IloAlgorithm::Feasible){
        std::cerr << "The answer is Feasible!\n";
        solved = -1;
    }
    else{
        if(status == IloAlgorithm::Optimal){
            std::cerr << "================================== WARNING ==================================\n";
            std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
            std::cerr << "The answer is Optimal!\n";
            std::cerr << "CHECK if LB == UB!\n";
            solved = -2;
        }
        else{
            IloCplex::Status detailed_status = cplex.getCplexStatus();
            if(detailed_status == IloCplex::AbortTimeLim){
                std::cerr << "The time limit is reached!\n";
                solved = 2;
            }
            else{
                if(status == IloAlgorithm::Infeasible){
                    std::cerr << "The CEMA instance is solved successfully!\n";
                    solved = 1;
                }
                else    
                    solved = 0;
            }
        }
    }
    return solved;
}

// 得到最优的目标
int ILP::get_s(){
    return cplex.getValue(s);
}

// 得到答案
int ILP::get_ans(){
    return y;
}

// 判断是不是邻居, 并返回边的编号
int iscon(int u, int w)
{
    int l = Graph.pstart[u], r = Graph.pstart[u + 1];
    if (l >= r)
        return -1;
    while (l + 1 < r)
    {
        int mid = l + (r - l) / 2;
        if (Graph.edges[mid].v > w)
            r = mid;
        else
            l = mid;
    }
    if(Graph.edges[l].v == w)
        return Graph.edges[l].tag;
    else    
        return -1;
}
