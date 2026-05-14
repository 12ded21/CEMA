#include "ilp.h"

ILP::ILP(int _m, int k, int z_num, int _LB, int _UB){
    m = _m;
    y = _UB - 1;
    k = std::min(k, _m);
    model = IloModel(env);
    x_e = IloBoolVarArray(env, m);
    z = IloBoolVarArray(env, z_num);
    z_cnt = 0;
    s = IloIntVar(env, 0, k);
    model.add(IloMinimize(env, s));
    model.add(IloSum(x_e) == s);
    cplex = IloCplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setError(env.getNullStream());
    cplex.setParam(IloCplex::Param::MIP::Display, 2);
    cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.use(new (env) LazyConstraintCallback(env, x_e, z, s, y));
}

ILP::~ILP(){
    env.end();
}

ILP::LazyConstraintCallback::LazyConstraintCallback(IloEnv env, IloBoolVarArray x_e, IloBoolVarArray z, IloIntVar s, int& y)
            : IloCplex::LazyConstraintCallbackI(env), _x_e(x_e), _z(z), _s(s), _env(env), _y(y){}

void ILP::LazyConstraintCallback::main(){
    std::cerr << "================================== Callback function begin ==================================\n";
    int m = _x_e.getSize();
    std::vector<int> del_edges;
    for (int i = 0; i < m; ++i){
        if (round(getValue(_x_e[i])) == 1){
            del_edges.push_back(i);
        }
    }

    for (int x : del_edges)
        Graph.ban_edge[x] = 1;

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
    std::vector<int> mxclq_edges;
    std::vector<int> neis;
    std::vector<std::vector<int>> clqs;

    std::sort(clq.begin(),clq.end(),std::less<int>());
    for(int m = 0; m < clq.size(); m++){
        for(int k = m + 1; k < clq.size(); k++){
            int edge_tag = iscon(clq[m], clq[k]);
            if(edge_tag >= 0)
                mxclq_edges.push_back(edge_tag);
        }
    }

    int temp = _y;
    if((int)clq.size() <= _y)
        _y = (int)clq.size() - 1;
    now_ans = _y + 1;
    if(_y <= 0 || _y < LB){
        Graph.global_delete_edges = del_edges;
        return ;
    }
        
    if (temp > _y)
    {
        Graph.global_delete_edges = del_edges;

        for (auto clq_ : Graph.clqs_edge){
            IloExpr lazyConstraint(_env);
            for(auto edges_ : clq_){
                lazyConstraint += _x_e[edges_];
            }
            add(lazyConstraint >= func(_y+1, ceil(sqrt(2 * clq_.size()))));
            lazyConstraint.end();
            all_constraints++;
        }
        std::cerr << "The clqs_z: " << Graph.clqs_z.size() << "\n";

        for(auto clq_z : Graph.clqs_z){
            IloExpr lazyConstraint(_env);
            for(auto z_tag_ : clq_z){
                lazyConstraint += _z[z_tag_];
            }
            add(lazyConstraint <= _y);
            lazyConstraint.end();
            all_constraints++;
        }
    }
    std::vector<int> edges_clq;
    std::vector<std::vector<int>> edges_clqs;

    if (clq.size() > _y){
        IloExpr lazyConstraint(_env);
        for (int i : mxclq_edges)
        {
            lazyConstraint += _x_e[i];
        }
        add(lazyConstraint >= func(_y+1, clq.size()));
        lazyConstraint.end();
        Graph.clqs_edge.push_back(mxclq_edges);
        all_constraints++;
    }
    std::cerr << "The num of constraints: " << getNrows() << std::endl;


    for (int x : del_edges)
        Graph.ban_edge[x] = 0;
    std::cerr << "callback: The s is: " << round(getValue(_s)) << std::endl;
    std::cerr << "callback: The y is: " << _y << std::endl;
}

void ILP::add_edge_clq(std::vector<int> &clq){
    IloExpr sum_cut_in_clq(env);
    for(int i : clq){
        sum_cut_in_clq += x_e[i];
    }
    model.add(sum_cut_in_clq >= func(y+1, ceil(sqrt(2 * clq.size()))));
    sum_cut_in_clq.end();
    all_constraints++;
}

void ILP::add_edge_clqs(std::vector<std::vector<int>> &clqs)
{
    for (int i = 0; i < clqs.size(); ++i)
        add_edge_clq(clqs[i]);
}

void ILP::add_node_clq(std::vector<int> &clq){
#ifdef EXTRA
    std::vector<int> z_vec;

    IloExpr sum_z(env);
    for(int i = 0; i < clq.size(); i++){
        int u = clq[i];
        IloExpr z_plus_sum_edge(env);
        for(int j = i + 1; j < clq.size(); j++){
            int v = clq[j];
            int edge_tag = iscon(u, v);
            if(edge_tag >= 0)
                z_plus_sum_edge += x_e[edge_tag];
            else
                std::cerr << "There is an error when get the node clq\n";  
        }
        z_plus_sum_edge += z[z_cnt];
        model.add(z_plus_sum_edge >= 1);
        z_plus_sum_edge.end();
        all_constraints++;
        sum_z += z[z_cnt];
        z_vec.push_back(z_cnt);
        z_cnt++;
    }
    model.add(sum_z <= y);
    sum_z.end();
    Graph.clqs_z.push_back(z_vec);
    all_constraints++;
#endif

}

void ILP::add_node_clqs(std::vector<std::vector<int>> &clqs)
{
    #ifdef EXTRA
    for (int i = 0; i < clqs.size(); ++i){
        std::sort(clqs[i].begin(), clqs[i].end(), [](int a, int b){
            return Graph.pstart[a+1] - Graph.pstart[a] > Graph.pstart[b+1] - Graph.pstart[b];
        });

        add_node_clq(clqs[i]);

        std::sort(clqs[i].begin(), clqs[i].end(), [](int a, int b){
            return Graph.pstart[a+1] - Graph.pstart[a] < Graph.pstart[b+1] - Graph.pstart[b];
        });

        add_node_clq(clqs[i]);
    }
    #endif
}

int ILP::solve(){
    int solved = 0;
    cplex.solve();
    
    IloAlgorithm::Status status = cplex.getStatus();
    ilpnode = cplex.getNnodes();
    ilpcut = cplex.getNcuts(IloCplex::CutUser);

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

int ILP::get_s(){
    return cplex.getValue(s);
}

int ILP::get_ans(){
    return y;
}

int iscon(int u, int w)
{
    if (u < 0 || w < 0 || u >= Graph.n || w >= Graph.n || u == w)
        return -1;
    if (u > w)
        std::swap(u, w);

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
