#ifndef SANTA_GLOBAL_INCLUDED
#define SANTA_GLOBAL_INCLUDED
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include "santa_par.h"

using namespace std;

// Global variables 
const long int max_child_happiness = (long long int) points_childIN * Nchprf * Nchild;
const long int max_gifts_happiness = (long long int) points_giftsIN * Ngfprf * Ntypes * Nqttys;
const long int child_gifts_happiness_ratio = max_gifts_happiness / max_child_happiness;

const double max_gifts_happiness_3 = pow((double) max_gifts_happiness, 3); // max_gifts_happiness to the power of 3
const long int child_gifts_happiness_ratio_3 = child_gifts_happiness_ratio * child_gifts_happiness_ratio * child_gifts_happiness_ratio;  // child_gifts_happiness_ratio to the power of 3

#ifdef ZERO_GIFT_R
const HappyInt gift_r = 0;  // weight of gift happiness point in approx. happiness 
#ifdef MERGED
const HappyInt single_r = 1; // effective happiness ratio r^3 C^2 / G^2, 
#else
const HappyInt single_r = 6;
#endif // MERGED
const HappyInt eff_r = 1;

#else // !ZERO_GIFT_R
#ifdef MERGED
const HappyInt gift_r = 1;
#else
const HappyInt gift_r = 6;
#endif // MERGED
extern HappyInt eff_r;  // effective happiness ratio r^3 C^2 / G^2, 
extern HappyInt single_r;
#endif // ZERO_GIFT_R

#ifdef MERGED
extern HappyInt triplet_r;
extern HappyInt twin_r;
const HappyInt gift_triplet_r = gift_r;
const HappyInt gift_twin_r = gift_r;
#else
extern HappyInt triplet_r;
extern HappyInt twin_r;
const HappyInt gift_triplet_r = gift_r / 3;
const HappyInt gift_twin_r = gift_r / 2;
#endif // MERGED


extern vector<vector<int>> wlist_child;
#ifndef NO_GIFT_MAP
extern vector<vector<int>> wlist_gifts;
#endif // !NO_GIFT_MAP
extern vector<vector<int>> states;
extern vector<int> answ;

extern vector<vector<HappyInt>> happiness;  // happyness given gift for each child. happiness[child][gift] == approx_happiness

// Global functions
void error(string errinfo, int nnum=0, double num1=0, double num2=0);
void read_csv(const char* filename, vector<vector<int>>& matrix, int header= 0, int index_col= 0);
double init();
int triplet_consistent(vector<int> &answ);
int twin_consistent(vector<int> &answ);
int triplet_twin_order();
bool check_gift_count(vector<int> &answ);

void reset_happiness_map(HappyInt non);
void update_eff_r(HappyInt total_child_happiness_points, HappyInt total_gift_happiness_points);

void build_happines_all();
void build_app_triplet_happiness_map();
void build_app_twin_happiness_map();
#ifndef NO_SINGLE_MAP
void build_app_single_happiness_map();
void build_app_gift_happiness_map();
#endif // !NO_SINGLE_MAP

inline long long int raw_single_happiness(int ichild, int igift);
long long int cal_child_all(vector<int> &answ);
long long int cal_gift_all(vector<int> &answ);
double cal_all(vector<int> &answ);

void states_to_answ();
void answ_to_states();

void write_conf(vector<int> &answ, const char* tag);
void write_conf(double best_happy, bool isrestart);

#endif // SANTA_GLOBAL_INCLUDED
