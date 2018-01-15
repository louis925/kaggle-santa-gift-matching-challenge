#ifndef SANTA_GLOBAL_INCLUDED
#define SANTA_GLOBAL_INCLUDED
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include "santa_par.h"

using namespace std;

// Global variables 
const long int max_child_happiness= (long long int) points_childIN * Nchprf * Nchild;
const long int max_gifts_happiness= (long long int) points_giftsIN * Ngfprf * Ntypes * Nqttys;
const long int child_gifts_happiness_ratio = max_gifts_happiness / max_child_happiness;

const double max_gifts_happiness_3 = pow((double) max_gifts_happiness, 3); // max_gifts_happiness to the power of 3
const long int child_gifts_happiness_ratio_3 = child_gifts_happiness_ratio*child_gifts_happiness_ratio*child_gifts_happiness_ratio;  // child_gifts_happiness_ratio to the power of 3

extern vector<vector<int>> wlist_child;
extern vector<vector<int>> wlist_gifts;
extern vector<vector<int>> states;
extern vector<vector<int>> gift_list;

extern vector<int> unhappy_child;

extern vector<vector<int>> child_happy_map; // Map of happyness given gift for each child. child_happy_map[child][gift] == child_happiness
extern vector<vector<int>> gift_happy_map; // Map of happyness given child for each gift. gift_happy_map[gift][child] == gift_happiness
//extern vector<unordered_map<int, int>> child_happy_map;  // Map of happyness given gift for each child. child_happy_map[child][gift] == child_happiness
//extern vector<unordered_map<int, int>> gift_happy_map; // Map of happyness given child for each gift. gift_happy_map[gift][child] == gift_happiness

// Global functions
inline double ran_generator();
void error(string errinfo, int nnum=0, double num1=0, double num2=0);
void read_csv(const char* filename, vector<vector<int>>& matrix, int header= 0, int index_col= 0);
double init();
inline int cal_single_child(int ichild, int igift);
inline int cal_single_gifts(int ichild, int igift);
inline long long int raw_single_happiness(int ichild, int igift);
long long int cal_child_all();
long long int cal_gift_all();
double cal_all();
void write_conf(double best_happy, bool isrestart= false);

void build_child_happiness_map();
void build_gift_happiness_map();

#endif // SANTA_GLOBAL_INCLUDED
