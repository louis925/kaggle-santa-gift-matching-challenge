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

const HappyInt gift_r = 6;
const HappyInt gift_twin_r = gift_r / 2;
const HappyInt gift_triplet_r = gift_r / 3;
const HappyInt single_r = 6;
const HappyInt twin_r = single_r / 2;
const HappyInt triplet_r = single_r / 3;
#ifndef FIX_EFF_R
extern Eff_R eff_r;  // effective happiness ratio r^3 C^2 / G^2, 
extern Eff_R happy_min;
#else
const Eff_R happy_min = points_childOUT * single_r * eff_r + points_giftsOUT * gift_r;
#endif


extern vector<vector<short>> wlist_child;
#ifndef NO_GIFT_MAP
extern vector<vector<int>> wlist_gifts;
#endif // !NO_GIFT_MAP
//extern vector<vector<int>> states;
extern vector<short> answ;
extern vector<vector<HappyInt>> happiness_c;
extern vector<vector<HappyInt>> happiness_g;


// Global functions
void init();
void reset_happiness_map(vector<vector<HappyInt>>& happiness, HappyInt non);
Eff_R update_eff_r(long long int total_child_happiness_points, long long int total_gift_happiness_points);

void build_happines_all(vector<vector<HappyInt>>& happiness_c, vector<vector<HappyInt>>& happiness_g);
void build_app_triplet_happiness_map(vector<vector<HappyInt>>& happiness_c);
void build_app_twin_happiness_map(vector<vector<HappyInt>>& happiness_c);
void build_app_single_happiness_map(vector<vector<HappyInt>>& happiness_c);
void build_app_gift_happiness_map(vector<vector<HappyInt>>& happiness_g);

long long int cal_child_all(const vector<short>& answ, const vector<vector<HappyInt>>& happiness_c);
long long int cal_gift_all(const vector<short>& answ, const vector<vector<HappyInt>>& happiness_g);
double cal_all_from_points(long long int C, long long int G);

int cal_single_child_ori(int c, int g);
int cal_single_gifts_ori(int c, int g);
long long int cal_child_all_ori(const vector<short>& answ);
long long int cal_gift_all_ori(const vector<short>& answ);
double cal_all_ori(const vector<short>& answ);

void states_to_answ();                           // Update answ from states
void answ_to_states(const vector<short>& answ);  // Update states from answ

int check_states();
int triplet_twin_order();
int triplet_consistent(const vector<short>& answ);
int twin_consistent(const vector<short>& answ);
int check_negative(const vector<short>& answ);
bool check_gift_count(const vector<short>& answ);
int check_unassigned_gifts(const vector<int>& unassigned_gifts);
int compare_answ_diff(const vector<short>& answ1, const vector<short>& answ2);
void print_vector(vector<int>& vec);

void error(string errinfo, int nnum = 0, double num1 = 0, double num2 = 0);
void read_csv(const char* filename, vector<vector<int>>& matrix, int header, int index_col);
void read_csv(const char* filename, vector<vector<short>>& matrix, int header, int index_col);
void write_conf(vector<short>& answ, const char* tag);
void write_conf(vector<short>& answ, double best_happy, bool isrestart);
#endif // SANTA_GLOBAL_INCLUDED
