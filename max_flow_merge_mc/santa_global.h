#ifndef SANTA_GLOBAL_INCLUDED
#define SANTA_GLOBAL_INCLUDED
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include "santa_par.h"

using namespace std;

// Global variables 
const long long int max_child_happiness = (long long int) points_childIN * Nchprf * Nchild;
const long long int max_gifts_happiness = (long long int) points_giftsIN * Ngfprf * Ntypes * Nqttys;
const int child_gifts_happiness_ratio = max_gifts_happiness / max_child_happiness;

const double max_gifts_happiness_3 = pow((double) max_gifts_happiness, 3); // max_gifts_happiness to the power of 3
const int child_gifts_happiness_ratio_3 = child_gifts_happiness_ratio * child_gifts_happiness_ratio * child_gifts_happiness_ratio;  // child_gifts_happiness_ratio to the power of 3

#if defined(ZERO_GIFT_R)
const Eff_R eff_r = 1;
const HappyInt gift_r = 0;  // weight of gift happiness point in approx. happiness 
const HappyInt gift_twin_r = 0;
const HappyInt gift_triplet_r = 0;
const HappyInt single_r = 6;
const HappyInt twin_r = 3;
const HappyInt triplet_r = 2;
#elif defined(FIX_EFF_R)
const Eff_R gift_r = 6;  // weight of gift happiness point in approx. happiness 
const HappyInt gift_twin_r = 3;
const HappyInt gift_triplet_r = 2;
const HappyInt single_r = 6 * eff_r;
const HappyInt twin_r = 3 * eff_r;
const HappyInt triplet_r = 2 * eff_r;
#elif defined(NO_FIX)
extern Eff_R eff_r;  // effective happiness ratio r^3 C^2 / G^2, 
const HappyInt gift_r = 1;
const HappyInt gift_twin_r = 1;
const HappyInt gift_triplet_r = 1;
extern HappyInt single_r;
extern HappyInt twin_r;
extern HappyInt triplet_r;
//#elif defined(MERGED_NO_FIX)
#else  // Variable eff_r case
extern Eff_R eff_r;  // effective happiness ratio r^3 C^2 / G^2, 
const HappyInt gift_r = 6;
const HappyInt gift_twin_r = gift_r / 2;
const HappyInt gift_triplet_r = gift_r / 3;
const HappyInt single_r = 6;
const HappyInt twin_r = single_r / 2;
const HappyInt triplet_r = single_r / 3;
#endif


extern vector<vector<short>> wlist_child;
#ifndef NO_GIFT_MAP
extern vector<vector<int>> wlist_gifts;
#endif // !NO_GIFT_MAP
extern vector<vector<int>> states;
extern vector<int> answ;


// Global functions
void init();
int check_states();
int triplet_consistent(const vector<int>& answ);
int twin_consistent(const vector<int>& answ);
int new_gift_inconsistent(const vector<int>& new_gift, const vector<int>& child_swap, vector<int>& inconsistent_child);
int triplet_twin_order();
bool check_gift_count(const vector<int>& answ);

void reset_happiness_map(vector<vector<HappyInt>>& happiness, HappyInt non);
Eff_R update_eff_r(long long int total_child_happiness_points, long long int total_gift_happiness_points);

void build_happines_all(vector<vector<HappyInt>>& happiness_c, vector<vector<HappyInt>>& happiness_g);
void build_app_triplet_happiness_map(vector<vector<HappyInt>>& happiness_c);
void build_app_twin_happiness_map(vector<vector<HappyInt>>& happiness_c);
#ifndef NO_SINGLE_MAP
void build_app_single_happiness_map(vector<vector<HappyInt>>& happiness_c);
void build_app_gift_happiness_map(vector<vector<HappyInt>>& happiness_g);
#endif // !NO_SINGLE_MAP

//inline HappyInt raw_single_happiness(int ichild, int igift);
long long int cal_child_all(const vector<int>& answ, const vector<vector<HappyInt>>& happiness_c);
//long long int cal_child_all(const vector<int>& answ);
long long int cal_gift_all(const vector<int>& answ, const vector<vector<HappyInt>>& happiness_g);
//long long int cal_gift_all(const vector<int>& answ);
//double cal_all(const vector<int>& answ, const vector<vector<HappyInt>>& happiness_c, const vector<vector<HappyInt>>& happiness_g);
//double cal_all(const vector<int>& answ);
double cal_all_from_points(long long int C, long long int G);

//bool verify_app_happy_points(vector<int>& answ, vector<vector<HappyInt>>& happiness,
//	long long int C, long long int G);
int compare_answ_diff(const vector<int>& answ1, const vector<int>& answ2);

void states_to_answ();                   // Update answ from states
void answ_to_states(const vector<int>& answ);  // Update states from answ

void error(string errinfo, int nnum=0, double num1=0., double num2=0.);
//void error(string errinfo, int nnum = 0, long long num1 = 0, long long num2 = 0);

void read_csv(const char* filename, vector<vector<int>>& matrix, int header, int index_col);
void read_csv(const char* filename, vector<vector<short>>& matrix, int header, int index_col);

void write_conf(vector<int>& answ, const char* tag);
void write_conf(double best_happy, bool isrestart);


// Old one
int cal_single_child_ori(int c, int g);
int cal_single_gifts_ori(int c, int g);
long long int cal_child_all_ori(const vector<int>& answ);
long long int cal_gift_all_ori(const vector<int>& answ);
double cal_all_ori(const vector<int>& answ);
#endif // SANTA_GLOBAL_INCLUDED
