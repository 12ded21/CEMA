// 包含读入读出函数、类型转换函数
// 全局变量的定义
#pragma once
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <set>
#define _BITSET_ //use bit set to represent the adjacency matrix
#define _STATISTIC_
#define _KERNEL_
#define _RECOLOR_
// #define NDEBUG
#include <cassert>

// #define NAIVE
#ifdef _BITSET_
#define set_bit(array, pos) (((array)[(pos)>>3]) |= (1<<((pos)&7)))
#define reverse_bit(array, pos) (((array)[(pos)>>3]) ^= (1<<((pos)&7)))
#define test_bit1(array, pos) (((array)[(pos)>>3])&(1<<((pos)&7)))
#define test_bit(array, pos) ((((array)[(pos)>>3])>>((pos)&7))&1)
#else
#define set_bit(array, pos) (((array)[pos]) = 1)
#define reverse_bit(array, pos) (((array)[pos]) = 1- ((array)[pos]))
#define test_bit(array, pos) ((array)[pos])
#endif

#ifndef NAIVE
#define MY_SOLVER
#endif

using ept = unsigned long;
using ll = long long;
using VertexSet = std::vector<int>;
using Edge = std::pair<int, int>;
using EdgeSet = std::vector<Edge>;

#define pb push_back
#define mp std::make_pair

// type conversion, integer to string
std::string integeter_to_str(ll);
//open the file
FILE* open_file(const char* , const char*);
//
void fread_wall(void* , size_t , size_t , FILE*);
// 计算 f 函数
int func(int, int);
// 初始化UB控制每次的个数
int f1(int budget, int clq_size);
int f2(int budget, int clq_size, int del, bool flag);
int f3(int budget, int clq_size, int del, bool flag);