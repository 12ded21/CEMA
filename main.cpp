#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ilcplex/ilocplex.h>
#include <cassert>
#include <time.h>
#include <filesystem>
#include <thread>
#include <cstdlib>
#include <filesystem>
#include "utils.h"
#include "ilp.h"
#include "cliques.h"
#include "get_mxclq.h"

#define BAN_2
#define DEGREE
#define COLOR
// #define CLIQUE is in cliques.h

using namespace std;
namespace fs = std::filesystem;

//The graph
EdgeListGraph Graph;
// the initial number of nodes and edges
//int initial_n, initial_m;
int *degree;
int *dedges;
// time
clock_t start_time;
clock_t end_time;

// The lower bound and the upper bound of the problem.
int LB = -1;
int UB = -1;

// 定义全局变量cishu
int cishu;
extern int Dynamic_Radio;
extern int INIT_CLIQUE;
extern int SHOW_COUNT;

// Parameter storage variables
std::string inputFile;
std::string outputFile;
double Time_Limit = -1.0;
int K = -1;
double Density = -1.0;
bool useStdout = true;

// Container for positional arguments
std::vector<std::string> positionalArgs;

// 输出参数
int initial_n = -1, initial_m = -1, red_n = -1, red_m = -1;
int color_ban_m = -1, clique_ban_m = -1;
int temp_constraints = 0;
int now_ans = -1;
int cut_num = 0;
// int node_num = 0;
// 总共添加的约束个数
long long all_constraints = 0;

double global_time = -1, red_time = -1, ILP_time = -1, LB_time = -1;
double color_ban_time = -1, clique_ban_time = -1;

void monitor_thread(int timeout_sec) {
    std::this_thread::sleep_for(std::chrono::seconds(timeout_sec));
    std::cerr << "Timeout! Forcing termination." << std::endl;
    fs::path pathObj(positionalArgs[0]);
    string graph_name = pathObj.filename().string();
    if(useStdout){
        cerr    << graph_name << "\t" << initial_n << "\t" << initial_m << "\t" << red_n << "\t" << red_m << "\t" << red_time << "\t" << color_ban_m << "\t" 
                << color_ban_time << "\t" << clique_ban_m << "\t" << clique_ban_time << "\t" 
                << all_constraints << "\t" 
                << cut_num << "\t" << ILP_time << "\t" << LB << "\t" << LB_time << "\t" << std::max(LB, now_ans) << "\t" << global_time << "\n";
    }
    else{
        std::ofstream outFile(outputFile, std::ios::app);
        outFile << graph_name << "\t" << initial_n << "\t" << initial_m << "\t" << red_n << "\t" << red_m << "\t" << red_time << "\t" << color_ban_m << "\t" 
                << color_ban_time << "\t" << clique_ban_m << "\t" << clique_ban_time << "\t"
                << all_constraints << "\t" 
                << cut_num << "\t" << ILP_time << "\t" << LB << "\t" << LB_time << "\t" << std::max(LB, now_ans) << "\t" << global_time << "\n";
    }
    std::exit(1); 
}


int main(int argc, char *argv[]){

    // Parse command-line arguments
    int method = 2; 
    int edge_select_method = 5; // 【新增】默认为 5 (随机选择)

    int i = 1;
    while (i < argc) {
        if (argv[i][0] == '-') {  // Option argument
            if (strcmp(argv[i], "-l") == 0) {
                if (i + 1 < argc) {
                    Time_Limit = std::atof(argv[i + 1]);
                    i += 2;  // Skip option and value
                } else {
                    std::cerr << "Error: -l option requires a value" << std::endl;
                    return 1;
                }
            } 
            else if (strcmp(argv[i], "-k") == 0) {
                if (i + 1 < argc) {
                    K = std::atoi(argv[i + 1]);
                    i += 2;
                } else {
                    std::cerr << "Error: -k option requires a value" << std::endl;
                    return 1;
                }
            } 
            else if (strcmp(argv[i], "-d") == 0) {
                if (i + 1 < argc) {
                    Density = std::atof(argv[i + 1]);
                    i += 2;
                } else {
                    std::cerr << "Error: -d option requires a value" << std::endl;
                    return 1;
                }
            } 
            else if (strcmp(argv[i], "-m") == 0) {
                if (i + 1 < argc) {
                    method = std::atoi(argv[i + 1]);
                    i += 2;
                } else {
                    std::cerr << "Error: -m option requires a value (1-5)" << std::endl;
                    return 1;
                }
            }
            // 【新增】解析 -e 参数，用于选边策略
            else if (strcmp(argv[i], "-e") == 0) {
                if (i + 1 < argc) {
                    edge_select_method = std::atoi(argv[i + 1]);
                    if (edge_select_method < 1 || edge_select_method > 5) {
                        std::cerr << "Error: -e option requires a value between 1 and 5" << std::endl;
                        return 1;
                    }
                    i += 2;
                } else {
                    std::cerr << "Error: -e option requires a value (1-5)" << std::endl;
                    return 1;
                }
            }
            else {
                std::cerr << "Error: Unknown option " << argv[i] << std::endl;
                return 1;
            }
        } 
        else {  // Positional argument
            positionalArgs.push_back(argv[i]);
            i++;
        }
    }

    // Process positional arguments...
    if (positionalArgs.empty()) {
        std::cerr << "Error: Missing required input file" << std::endl;
        return 1;
    }
    inputFile = positionalArgs[0];
    
    if (positionalArgs.size() >= 2) {
        outputFile = positionalArgs[1];
        useStdout = false;
    }

    // Validate required parameters
    if (Time_Limit <= 0) {
        std::cerr << "Error: Time limit (-l) must be positive" << std::endl;
        return 1;
    }
    
    // Check mutually exclusive parameters
    if (K >= 0 && Density >= 0) {
        std::cerr << "Error: -k and -d options cannot be used together" << std::endl;
        return 1;
    }
    if (K < 0 && Density < 0) {
        std::cerr << "Error: Must provide either -k or -d parameter" << std::endl;
        return 1;
    }

    // Validate parameter ranges
    if (K >= 0 && K <= 0) {
        std::cerr << "Error: Edge budget (-k) must be positive integer" << std::endl;
        return 1;
    }
    if (Density >= 0 && (Density < 0 || Density > 1)) {
        std::cerr << "Error: k density (-d) must be in range [0,1]" << std::endl;
        return 1;
    }

    // Print configuration
    std::cout << "=== Configuration Parameters ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output destination: " 
              << (useStdout ? "stdout" : outputFile) << std::endl;
    std::cout << "Time limit: " << Time_Limit << "s" << std::endl;
    if (K >= 0) {
        std::cout << "Edge budget (k): " << K << std::endl;
    } else {
        std::cout << "k density (d): " << Density << std::endl;
    }
    // 【新增】打印选取的边排序方法
    std::cout << "Edge select method: " << edge_select_method << " (1:CN Desc, 2:CN Asc, 3:DegSum Desc, 4:DegSum Asc, 5:Random)" << std::endl;

    // Example processing logic
    if (useStdout) {
        // std::cout << "\nProcessing complete. Output to stdout..." << std::endl;
        // Actual processing code
    } else {
        std::ofstream outFile(outputFile, std::ios::app);
        if (outFile) {
            // outFile << "Processing complete. Output to file: " 
            //         << outputFile << std::endl;
            // Actual processing code
        } else {
            std::cerr << "Error: Could not open output file" << std::endl;
            return 1;
        }
    }


    // 开启线程以强制结束程序
    std::thread watchdog(monitor_thread, Time_Limit);

    // 开始程序
    int ordering = -1, _all = FALSE;

    // DIMACS-2, 所有图都以DIMACS-2格式读入
    Graph.readDIMACS2Text(inputFile.c_str());

    // 设置随机种子
    srand(time(0));

    // 开始记录程序运行总时间
    auto global_start = std::chrono::high_resolution_clock::now();

    initial_n = Graph.n;
    initial_m = Graph.m / 2;

    std::vector<int> tep = final_max_clique_nodes();
    std::cerr << tep.size() << std::endl;
    

    if(Density >= 0){
        K = ceil(Density * (double)Graph.m / 2);
    }
    std::cerr << "K: " << K << "\n";
    std::cerr << "================================== Get LB begin ==================================" << "\n";
    auto LB_start = std::chrono::high_resolution_clock::now();
    // 染色，找到包含每条边的极大团，快速找到一个下界，只关心
    // 极大团的数目而不关心这个图案的具体边和点。
    LB = LB_with_fast_color(method);
    // LB = LB_with_fast_color();
    auto LB_end = std::chrono::high_resolution_clock::now();
    LB_time = std::chrono::duration_cast<std::chrono::microseconds>(LB_end - LB_start).count();
    LB_time = LB_time / 1000000;
    std::cerr << "Time of init_LB: " << LB_time << "\n";
    std::cerr << "LB: " << LB << std::endl;

    // fs::path pathObj(positionalArgs[0]);
    // string graph_name = pathObj.filename().string();
    // if(useStdout){
    //     cerr    << graph_name << "\t" << LB << "\t" << LB_time << "\n";
    // }
    // else{
    //     std::ofstream outFile(outputFile, std::ios::app);
    //     outFile << graph_name << "\t" << LB << "\t" << LB_time << "\n";
    // }
    // std::exit(0);

#ifdef DEGREE
    std::cerr << "================================== Degree Reduction begin ==================================" << "\n";
    auto degree_reduction_start = std::chrono::high_resolution_clock::now();


    Graph.degen();

    int temp_n = Graph.n;
    int temp_m = Graph.m;
    Graph.core_shrink_graph();

    auto degree_reduction_end = std::chrono::high_resolution_clock::now();
    red_time = std::chrono::duration_cast<std::chrono::microseconds>(degree_reduction_end - degree_reduction_start).count();
    red_time = red_time / 1000000;

    red_n = Graph.n;
    red_m = Graph.m / 2;
    std::cerr << "Time of reduction: " << red_time << "\n";
    std::cerr << "num of nodes: " << red_n << "(" << "reduced nodes: " << temp_n - Graph.n << ")\n";
    std::cerr << "num of edges: " << red_m << "(" << "reduced edges: " << (temp_m - Graph.m) / 2<< ")\n";
    
    // 发现所有的点都没了,插入一个结点
    if(initial_n != 0 && Graph.n == 0){
        Graph.n = 1;
        Graph.pstart[0] = Graph.pstart[1] = 0;
    }

#endif

// 启发式ban边
#ifdef BAN_2
    std::cerr << "================================== Ban edges ==================================" << "\n";
    std::cerr << "The density of graph: " << Graph.R << "\n";
    // 找到包含每个点的最大团，可以在开始把这些团加入到ILP中
    vector<pair<int, int>> all_edges;
    vector<int> store_edges;
    int color_ban = 0;
    int clique_ban = 0;
    int after_degree = Graph.m / 2;
#ifdef COLOR
    auto color_ban_start = std::chrono::high_resolution_clock::now();
    all_edges = color_vector();
    for(int i = 0; i < all_edges.size(); i++){
        if(all_edges[i].first <= LB){
            Graph.ban_edge[all_edges[i].second] = 1;
            color_ban++;
        }
        else{
            store_edges.push_back(all_edges[i].second);
        }
    }
    auto color_ban_end = std::chrono::high_resolution_clock::now();
    color_ban_time = std::chrono::duration_cast<std::chrono::microseconds>(color_ban_end - color_ban_start).count();
    color_ban_time = color_ban_time / 1000000;
    std::cerr << "Time of color_banning: " << color_ban_time << "\n";
#endif
#ifndef COLOR
    all_edges = color_vector();
    for(int i = 0; i < Graph.m / 2; i++){
        store_edges.push_back(all_edges[i].second);
    }
#endif
    int edges_num = 0;
    // n太大
    if( Graph.m / 2 > 0 ){
        auto clique_ban_start = std::chrono::high_resolution_clock::now();
        // 存放应当判断的边个数
        double temp_R = 0;
        if(Graph.m != 0)
            temp_R = 2 * (double)(Graph.m / 2 - color_ban) /(double)((double)Graph.n * ((double)Graph.n - 1));
        edges_num = (int)((1 - temp_R) / 4 * (double)Graph.m / 2);
        
        edges_num = std::min(edges_num, 600);
        edges_num = std::min(edges_num, (int)store_edges.size());

        // ================= 开始：5种选边排序方法 =================
        // 直接使用由命令行 argv 解析得到的 edge_select_method 变量
        if (edge_select_method == 5) {
            // 方法5：随机打乱顺序
            std::random_shuffle(store_edges.begin(), store_edges.end());
        } else {
            // 预计算边的权重以加速排序
            std::vector<std::pair<int, int>> edge_weights; // <权重, 边编号>
            edge_weights.reserve(store_edges.size());
            
            for (int edge_tag : store_edges) {
                int u = Graph.edge_to_nodes[edge_tag].first;
                int v = Graph.edge_to_nodes[edge_tag].second;
                int weight = 0;
                
                if (edge_select_method == 1 || edge_select_method == 2) {
                    weight = common_neighbor_num(u, v);
                } else if (edge_select_method == 3 || edge_select_method == 4) {
                    weight = Graph.degree[u] + Graph.degree[v];
                }
                edge_weights.push_back({weight, edge_tag});
            }

            // 根据方法进行排序
            if (edge_select_method == 1 || edge_select_method == 3) {
                // 方法1, 3: 从大到小降序排列
                std::sort(edge_weights.begin(), edge_weights.end(), std::greater<std::pair<int, int>>());
            } else {
                // 方法2, 4: 从小到大升序排列
                std::sort(edge_weights.begin(), edge_weights.end(), std::less<std::pair<int, int>>());
            }

            // 将排好序的边编号写回 store_edges
            for (size_t i = 0; i < store_edges.size(); i++) {
                store_edges[i] = edge_weights[i].second;
            }
        }
        // ================= 结束：5种选边排序方法 =================

        vector<int> ban_edges(store_edges.begin(), store_edges.begin() + edges_num);
        Max_clqs_of_edges(ban_edges);

        auto clique_ban_end = std::chrono::high_resolution_clock::now();
        clique_ban_time = std::chrono::duration_cast<std::chrono::microseconds>(clique_ban_end - clique_ban_start).count();
        clique_ban_time = clique_ban_time / 1000000;
    }
    else{
        clique_ban_time = 0;
    }
    std::cerr << "Time of clique_banning: " << clique_ban_time << "\n";
    std::cerr << "number of selected edges: " << edges_num << "\n";

    // 把边ban掉，直接把边删掉
    Graph.ban();

    color_ban_m = after_degree - color_ban;
    clique_ban_m = Graph.m / 2;

    clique_ban = after_degree - Graph.m / 2 - color_ban;

    std::cerr << "color_ban: " << color_ban << "\n";
    std::cerr << "clique_ban: " << clique_ban << "\n";
    std::cerr << "banned edges: " << after_degree - Graph.m / 2 << "\n";
    std::cerr << "reserved edges: " << Graph.m / 2 << "\n";
#endif

// 初始化，逐渐删去边，直到所有的 k 用完，每次从最大团删掉若干条边，并将每个最大团记录下来作为约束。
    std::cerr << "================================== UB initialize ==================================" << "\n";
    // output the graph
    // for(int i = 0; i < Graph.n; i++){
    //     std::cerr << "Node " << i << ": ";
    //     for(int j = Graph.pstart[i]; j < Graph.pstart[i + 1]; j++){
    //         std::cerr << Graph.edges[j].v << "(" << Graph.edges[j].tag << ") ";
    //     }
    //     std::cerr << "\n";
    // }
    // 设定预算。
    int budget = K;
    // 不是真的删边，所以在while之后需要复原。
    std::vector<int> temp_del_edges;
    // UB_clqs 存储每次找到的最大团
    std::vector<std::vector<int>> UB_clqs;
    // 记录UB
    int temp_UB = Graph.n;
    // 记录每次删边的数量
    int del_num = -1;


    while(budget > 0){
        std::vector<int> clq = final_max_clique_nodes();
        std::cerr << "Current max clique size: " << clq.size() << "\n";
        if(clq.size() <= 1)
            break;
        del_num = f3(budget, clq.size(), del_num, temp_UB > (int)clq.size() ? 1 : 0);
        del_num = std::min(del_num, budget);
        budget -= del_num;
        std::cerr << "Deleting " << del_num << " edges from the max clique.\n";

        temp_UB = std::min(temp_UB, (int)clq.size());
        // 把团的边弄成边的形式
        std::vector<int> clq_edges;
        for(int m = 0; m < clq.size(); m++){
            for(int k = m + 1; k < clq.size(); k++){
                int edge_tag = iscon(clq[m], clq[k]);
                if(edge_tag >= 0)
                    clq_edges.push_back(edge_tag);
            }
        }
        UB_clqs.push_back(clq_edges);
        // 随机删掉del_num条边
        std::random_shuffle(clq_edges.begin(), clq_edges.end());
        for(int i = 0; i < del_num; i++){
            Graph.ban_edge[clq_edges[i]] = 1;
            temp_del_edges.push_back(clq_edges[i]);
        }
        // std::cerr << "\n";
    }
    // 复原边
    for(int x : temp_del_edges)
        Graph.ban_edge[x] = 0;


    // 把得到的最大团加入到 ILP 开始进行 ILP 求解。
    auto ILP_start = std::chrono::high_resolution_clock::now();
    int solve_log;

    // 统计z变量开多大
    int zz = 0;

#ifdef EXTRA
    for(int i = 0; i < Graph.clqs_node.size(); i++){
        // // 其中两倍是因为按照两个顺序加了两次
        zz += 2 * Graph.clqs_node[i].size();

        // 单个顺序
        // zz += Graph.clqs_node[i].size();
    }
#endif

    ILP ilp_solver(Graph.m / 2, K, zz, LB, temp_UB);
    ilp_solver.add_edge_clqs(Graph.clqs_edge);
    // 加入上界的约束
    ilp_solver.add_edge_clqs(UB_clqs);

    // 添加额外的约束, 限制密度0.8
    #ifdef EXTRA 
    if(Graph.R >= 0.8){
        std::random_shuffle(Graph.clqs_node.begin(), Graph.clqs_node.end());
        // 随机取20个
        int extra_clq_num = std::min(20, (int)Graph.clqs_node.size());
        for(int i = 0; i < extra_clq_num; i++){
            ilp_solver.add_node_clq(Graph.clqs_node[i]);
        }
    }
    #endif

    solve_log = ilp_solver.solve();
    cerr << "================================== Solving information ==================================\n";
    cerr << "solving result: " << solve_log << endl;
    try{
        int final_s = ilp_solver.get_s();
        cerr << "The s is: " << final_s << "\n";
    } catch(const IloException& e){
        std::cerr << "ILP has no solution s!" << "\n";
    }
    cerr << "The ans is: " << ilp_solver.y + 1 << "\n";

    cerr << "delete "; cerr << Graph.global_delete_edges.size() << " "; cerr << "edges: " << "\n";
    for(auto x : Graph.global_delete_edges)
        cerr << x << " ";
    cerr << std::endl;

    auto ILP_end = std::chrono::high_resolution_clock::now();
    ILP_time = std::chrono::duration_cast<std::chrono::microseconds>(ILP_end - ILP_start).count();
    ILP_time = ILP_time / 1000000;
    cut_num = ilp_solver.ilpcut;
    // node_num = ilp_solver.ilpnode;
    std::cerr << "Time of ILP: " << ILP_time << "\n";

    auto global_end = std::chrono::high_resolution_clock::now();
    global_time = std::chrono::duration_cast<std::chrono::microseconds>(global_end - global_start).count();
    global_time = global_time / 1000000;
    std::cerr << "Time: " << global_time << "\n";
    // 检验最优解
    cerr << "================================== Check the solution ==================================\n";
    for(auto x : Graph.global_delete_edges)
        Graph.ban_edge[x] = 1;
    vector<int> temp;
    if(Graph.m / 2 - (int)Graph.global_delete_edges.size()> 0)
        temp =  final_max_clique_nodes();
    else   
        temp.push_back(0);
    if(temp.size() == ilp_solver.y + 1)
        std::cerr << "The ans is right!!!\n";
    else   
        std::cerr << "The ans is Wrong???\n";

    fs::path pathObj(positionalArgs[0]);
    string graph_name = pathObj.filename().string();
    if(useStdout){
        cerr    << graph_name << "\t" << initial_n << "\t" << initial_m << "\t" << red_n << "\t" << red_m << "\t" << red_time << "\t" << color_ban_m << "\t" 
                << color_ban_time << "\t" << clique_ban_m << "\t" << clique_ban_time << "\t"
                << all_constraints << "\t" 
                << cut_num << "\t" << ILP_time << "\t" << LB << "\t" << LB_time << "\t" << std::max(LB, ilp_solver.y + 1) << "\t" << global_time << "\n";
    }
    else{
        std::ofstream outFile(outputFile, std::ios::app);
        outFile << graph_name << "\t" << initial_n << "\t" << initial_m << "\t" << red_n << "\t" << red_m << "\t" << red_time << "\t" << color_ban_m << "\t" 
                << color_ban_time << "\t" << clique_ban_m << "\t" << clique_ban_time << "\t"
                << all_constraints << "\t" 
                << cut_num << "\t" << ILP_time << "\t" << LB << "\t" << LB_time << "\t" << std::max(LB, ilp_solver.y + 1) << "\t" << global_time << "\n";
    }
    watchdog.detach();
    return 0;
}
