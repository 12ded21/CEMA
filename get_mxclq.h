#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>
#include <assert.h>
#include "EdgeListGraph.h"
#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define P_TRUE 2
#define P_FALSE 0
#define NO_REASON -3
#define CONFLICT -1978
#define tab_node_size  510
#define max_expand_depth 1200
#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer
#define is_neibor(i,j) matrice[i][j]
#define candidate_start branch_stack[branch_stack_fill_pointer][0]
#define candidate_end Candidate_Stack_fill_pointer
#define branching_node branch_stack[branch_stack_fill_pointer][1]
#define candidate_is_not_empty() (Candidate_Stack_fill_pointer-candidate_start>0)
#define candidate_is_empty() (Candidate_Stack_fill_pointer==candidate_start)
#define CUR_CLQ_SIZE Clique_Stack_fill_pointer
#define NOT_ROOT (ptr(Cursor_Stack)>1)
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]
#define SSINDEX SIndex_Stack[SIndex_Stack_fill_pointer-1]
#define MIN(a,b) a<=b?a:b
#define time_in_s(clock) (double)(clock)/clock_per_second
#define fix_node(node,iset) ((node > NB_NODE)? fix_newNode_for_iset(node, iset):fix_oldNode_for_iset(node, iset))
#define assign_node(node, value, reason) {\
	node_state[node] = value;\
	node_reason[node] = reason;\
	push(node, FIXED_NODE_STACK);\
}

class GET{
public:

static int static_degree[tab_node_size];
static int none_degree[tab_node_size];

int FORMAT = 1, DENSITY, NB_NODE, ADDED_NODE, NB_EDGE, MAX_CLQ_SIZE,
		MAX_ISET_SIZE, NB_BACK_CLIQUE;
char node_state[2 * tab_node_size];
int node_reason[2 * tab_node_size];
char *static_matrix;
char matrice[tab_node_size][tab_node_size];
int iSET[tab_node_size][tab_node_size];
int iSET_COUNT = 0;
int iSET_Size[tab_node_size];
char iSET_State[tab_node_size];
char iSET_Used[tab_node_size];
char iSET_Tested[tab_node_size];
int iSET_Index[tab_node_size];
int REDUCED_iSET_STACK[tab_node_size * 10];
int REDUCED_iSET_STACK_fill_pointer = 0;
int PASSIVE_iSET_STACK[tab_node_size];
int PASSIVE_iSET_STACK_fill_pointer = 0;
int FIXED_NODE_STACK[2 * tab_node_size];
int FIXED_NODE_STACK_fill_pointer = 0;
int UNIT_STACK[tab_node_size];
int UNIT_STACK_fill_pointer = 0;
int NEW_UNIT_STACK[tab_node_size];
int NEW_UNIT_STACK_fill_pointer = 0;
int *node_neibors[tab_node_size];
int *none_neibors[tab_node_size];
int active_degree[tab_node_size];
int Clique_Stack[tab_node_size];
int Clique_Stack_fill_pointer = 0;
int *INIT_Stack;
int INIT_Stack_fill_pointer = 0;
int MaxCLQ_Stack[tab_node_size];
int Candidate_Stack[tab_node_size * max_expand_depth];
int Candidate_Stack_fill_pointer = 0;
int Cursor_Stack[max_expand_depth];
int Cursor_Stack_fill_pointer = 0;
int Vertex_UB[tab_node_size * max_expand_depth];
int Rollback_Point;
int Branching_Point;
int Tmp_Stack[tab_node_size * 2];
int Tmp_Stack_fill_pointer = 0;
int NEW_OLD[tab_node_size];
int OLD_NEW[tab_node_size];
int NB_CANDIDATE = 0, FIRST_INDEX, REBUILD_MATRIX = FALSE;

int Extra_Node_Stack[1000];
int Last_Idx = 0;
int cut_ver = 0, total_cut_ver = 0;
int cut_inc = 0, total_cut_inc = 0;
int cut_iset = 0, total_cut_iset = 0;
int cut_satz = 0, total_cut_satz = 0;
long long Branches_Nodes[6];
int STATIC_ORDERING = TRUE, LAST_IN, INIT_CLIQUE = 0;
float Dynamic_Radio = 0.6;
float Mean_Dynamic_Radio = 0.0;
int Dynamic_Count = 0;
int LIST_ALL = FALSE;
int Branches[1200];
int MAX_COUNT = 0, SHOW_COUNT = 0;
int *APPEND_STACK;
int APPEND_STACK_SIZE = 0;
int APPEND_STACK_USED = 0;
int CONFLICT_ISET_STACK[tab_node_size * tab_node_size];
int CONFLICT_ISET_STACK_fill_pointer;
int ADDED_NODE_iSET[2 * tab_node_size];
int iSET_Involved[tab_node_size];
int REASON_STACK[tab_node_size];
int REASON_STACK_fill_pointer = 0;
char tested[tab_node_size];
int BRANCHING_COUNT = 0, LAST_BRANCHING_COUNT = 0;
struct rusage lasttime;
int** flag;

public:

GET();
~GET();

static int none_degree_inc(const void *a, const void *b);
static int static_degree_dec(const void *a, const void *b);
float density(int nb_node, int nb_edge);

int further_test_reduced_iset(int start);
int inc_maxsatz_lookahead_by_fl2();
int inc_maxsatz_on_last_iset(int end);
int open_new_iset_old(int i);
int simple_further_test_node(int start);
int test_node_for_failed_nodes(int node, int iset);
int test_by_eliminate_failed_nodes();

int build_simple_graph_instance(EdgeListGraph& graph);
void search_initial_maximum_clique();
void complement_graph();
void sort_by_active_degree();

void my_sort_by_degree_dec();
void my_sort_by_degree_inc();
int iset_smaller_than(int iset1, int iset2);
void sort_isets_and_push_nodes();
int choose_candidate_node();

void print_isets(int i_set);
void print_all_iset();
void printallMaxClique();
void print_help();
void print_version();
void print_branches();

int addIntoIsetTomitaBis(int node);
int re_number(int node);
int fix_newNode_for_iset(int fix_node, int fix_iset);
int fix_oldNode_for_iset(int fix_node, int fix_iset);
int fix_anyNode_for_iset(int fix_node, int fix_iset);
int fix_node_iset(int fix_iset);

int cut_by_iset_last_renumber();
int cut_by_inc_ub(int print_info);
int cut_by_inc_maxsat_eliminate_first();
void rollback_context_for_maxsatz(int start_fixed, int start_passive, int start_reduced);
void reset_context_for_maxsatz();

void init_for_search();
void search_maxclique(int cutoff, int print_info);
void init_for_maxclique(int ordering, int list_all);
int sort_by_maxiset(int mandatory);
void rebuild_matrix(int start);

void identify_conflict_sets(int iset_idx);
void enlarge_conflict_sets();
void push_reason_iset(int iset_idx);
int unit_iset_process();
int unit_iset_process_used_first();

std::vector<int> get_mxclq(EdgeListGraph& graph);
void re_code();
void store_maximal_clique(int node);
void store_maximum_clique(int node, int print_info);

};
