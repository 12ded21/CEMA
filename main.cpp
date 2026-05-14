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

using namespace std;
namespace fs = std::filesystem;

EdgeListGraph Graph;
int *degree;
int *dedges;
clock_t start_time;
clock_t end_time;

int LB = -1;
int UB = -1;

int cishu;
extern int Dynamic_Radio;
extern int INIT_CLIQUE;
extern int SHOW_COUNT;

std::string inputFile;
std::string outputFile;
double Time_Limit = -1.0;
int K = -1;
double Density = -1.0;
bool useStdout = true;

std::vector<std::string> positionalArgs;

int initial_n = -1, initial_m = -1, red_n = -1, red_m = -1;
int color_ban_m = -1, clique_ban_m = -1;
int temp_constraints = 0;
int now_ans = -1;
int cut_num = 0;
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

    int method = 2; 
    int edge_select_method = 5;

    int i = 1;
    while (i < argc) {
        if (argv[i][0] == '-') {
            if (strcmp(argv[i], "-l") == 0) {
                if (i + 1 < argc) {
                    Time_Limit = std::atof(argv[i + 1]);
                    i += 2;
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
        else {
            positionalArgs.push_back(argv[i]);
            i++;
        }
    }

    if (positionalArgs.empty()) {
        std::cerr << "Error: Missing required input file" << std::endl;
        return 1;
    }
    inputFile = positionalArgs[0];
    
    if (positionalArgs.size() >= 2) {
        outputFile = positionalArgs[1];
        useStdout = false;
    }

    if (Time_Limit <= 0) {
        std::cerr << "Error: Time limit (-l) must be positive" << std::endl;
        return 1;
    }
    
    if (K >= 0 && Density >= 0) {
        std::cerr << "Error: -k and -d options cannot be used together" << std::endl;
        return 1;
    }
    if (K < 0 && Density < 0) {
        std::cerr << "Error: Must provide either -k or -d parameter" << std::endl;
        return 1;
    }

    if (K >= 0 && K <= 0) {
        std::cerr << "Error: Edge budget (-k) must be positive integer" << std::endl;
        return 1;
    }
    if (Density >= 0 && (Density < 0 || Density > 1)) {
        std::cerr << "Error: k density (-d) must be in range [0,1]" << std::endl;
        return 1;
    }

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
    std::cout << "Edge select method: " << edge_select_method << " (1:CN Desc, 2:CN Asc, 3:DegSum Desc, 4:DegSum Asc, 5:Random)" << std::endl;

    if (!useStdout) {
        std::ofstream outFile(outputFile, std::ios::app);
        if (!outFile) {
            std::cerr << "Error: Could not open output file" << std::endl;
            return 1;
        }
    }


    std::thread watchdog(monitor_thread, Time_Limit);

    int ordering = -1, _all = FALSE;

    Graph.readDIMACS2Text(inputFile.c_str());

    srand(time(0));

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
    LB = LB_with_fast_color(method);
    auto LB_end = std::chrono::high_resolution_clock::now();
    LB_time = std::chrono::duration_cast<std::chrono::microseconds>(LB_end - LB_start).count();
    LB_time = LB_time / 1000000;
    std::cerr << "Time of init_LB: " << LB_time << "\n";
    std::cerr << "LB: " << LB << std::endl;

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
    
    if(initial_n != 0 && Graph.n == 0){
        Graph.n = 1;
        Graph.pstart[0] = Graph.pstart[1] = 0;
    }

#endif

#ifdef BAN_2
    std::cerr << "================================== Ban edges ==================================" << "\n";
    std::cerr << "The density of graph: " << Graph.R << "\n";
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
    if( Graph.m / 2 > 0 ){
        auto clique_ban_start = std::chrono::high_resolution_clock::now();
        double temp_R = 0;
        if(Graph.m != 0)
            temp_R = 2 * (double)(Graph.m / 2 - color_ban) /(double)((double)Graph.n * ((double)Graph.n - 1));
        edges_num = (int)((1 - temp_R) / 4 * (double)Graph.m / 2);
        
        edges_num = std::min(edges_num, 600);
        edges_num = std::min(edges_num, (int)store_edges.size());

        if (edge_select_method == 5) {
            std::random_shuffle(store_edges.begin(), store_edges.end());
        } else {
            std::vector<std::pair<int, int>> edge_weights;
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

            if (edge_select_method == 1 || edge_select_method == 3) {
                std::sort(edge_weights.begin(), edge_weights.end(), std::greater<std::pair<int, int>>());
            } else {
                std::sort(edge_weights.begin(), edge_weights.end(), std::less<std::pair<int, int>>());
            }

            for (size_t i = 0; i < store_edges.size(); i++) {
                store_edges[i] = edge_weights[i].second;
            }
        }

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

    Graph.ban();

    color_ban_m = after_degree - color_ban;
    clique_ban_m = Graph.m / 2;

    clique_ban = after_degree - Graph.m / 2 - color_ban;

    std::cerr << "color_ban: " << color_ban << "\n";
    std::cerr << "clique_ban: " << clique_ban << "\n";
    std::cerr << "banned edges: " << after_degree - Graph.m / 2 << "\n";
    std::cerr << "reserved edges: " << Graph.m / 2 << "\n";
#endif

    std::cerr << "================================== UB initialize ==================================" << "\n";
    int budget = K;
    std::vector<int> temp_del_edges;
    std::vector<std::vector<int>> UB_clqs;
    int temp_UB = Graph.n;
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
        std::vector<int> clq_edges;
        for(int m = 0; m < clq.size(); m++){
            for(int k = m + 1; k < clq.size(); k++){
                int edge_tag = iscon(clq[m], clq[k]);
                if(edge_tag >= 0)
                    clq_edges.push_back(edge_tag);
            }
        }
        UB_clqs.push_back(clq_edges);
        std::random_shuffle(clq_edges.begin(), clq_edges.end());
        for(int i = 0; i < del_num; i++){
            Graph.ban_edge[clq_edges[i]] = 1;
            temp_del_edges.push_back(clq_edges[i]);
        }
    }
    for(int x : temp_del_edges)
        Graph.ban_edge[x] = 0;


    auto ILP_start = std::chrono::high_resolution_clock::now();
    int solve_log;

    int zz = 0;

#ifdef EXTRA
    for(int i = 0; i < Graph.clqs_node.size(); i++){
        zz += 2 * Graph.clqs_node[i].size();
    }
#endif

    ILP ilp_solver(Graph.m / 2, K, zz, LB, temp_UB);
    ilp_solver.add_edge_clqs(Graph.clqs_edge);
    ilp_solver.add_edge_clqs(UB_clqs);

    #ifdef EXTRA 
    if(Graph.R >= 0.8){
        std::random_shuffle(Graph.clqs_node.begin(), Graph.clqs_node.end());
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
    std::cerr << "Time of ILP: " << ILP_time << "\n";

    auto global_end = std::chrono::high_resolution_clock::now();
    global_time = std::chrono::duration_cast<std::chrono::microseconds>(global_end - global_start).count();
    global_time = global_time / 1000000;
    std::cerr << "Time: " << global_time << "\n";
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
