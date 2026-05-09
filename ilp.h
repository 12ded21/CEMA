#pragma once
#include <ilcplex/ilocplex.h>
#include <vector>
#include <cassert>
#include <cmath>
#include "EdgeListGraph.h"
#include "utils.h"
#include "cliques.h"

extern EdgeListGraph Graph;
extern int LB;
extern int now_ans;
extern long long all_constraints;

class ILP{
public:
    IloEnv env;
    IloModel model;
    IloBoolVarArray x_e;
    IloBoolVarArray z;
    IloIntVar s;
    IloCplex cplex;
    // 定义ILP树上经过的节点数和剪枝次数
    int ilpnode, ilpcut;
    int m;
    int y;
    int z_cnt = 0;
    // // 记录所有加过的约束,记录所有团
    // std::vector<std::vector<int>> all_constraints;
    // 懒惰约束回调函数
    class LazyConstraintCallback: public IloCplex::LazyConstraintCallbackI{
    private:
        IloBoolVarArray _x_e;
        IloBoolVarArray _z;
        IloIntVar _s;
        IloEnv _env;
        int &_y;

    public:
        LazyConstraintCallback(IloEnv, IloBoolVarArray,  IloBoolVarArray, IloIntVar, int&);
        void main() override ;
        IloCplex::CallbackI *duplicateCallback() const override{
            return new (_env) LazyConstraintCallback(_env, _x_e, _z, _s, _y);
        }
    };

public:
    // _m表示变量个数，k表示限制删掉的边数，_LB和_UB来确定变量的范围。
    ILP(int, int, int, int, int);
    ~ILP();
    // 添加单个团，团中存储的是边
    void add_edge_clq(std::vector<int> &);
    // 批量添加团，团中存储的是边
    void add_edge_clqs(std::vector<std::vector<int>> &);
    // 添加单个团，团中存储的是点
    void add_node_clq(std::vector<int> &);
    // 批量添加团，团中存储的是点
    void add_node_clqs(std::vector<std::vector<int>> &);
    // 求解
    int solve();
    // 得到y
    int get_ans();
    // 得到s
    int get_s();
};

// 判断 u 和 w 是不是邻居
int iscon(int u, int w);