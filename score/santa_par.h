#ifndef SANTA_PAR_INCLUDED
#define SANTA_PAR_INCLUDED

// System
const int Nchild=   1000000; // N of children
const int Ntypes=      1000; // N of different types of gifts
const int Nqttys=      1000; // quantity of one gifts
const int Ntripl=      5001; // N of triplet children
const int Ntwins=     45001; // N of twin and triplet children
const int Nchprf=       100; // length of child wish list
const int Ngfprf=      1000; // length of gifts wish list

// Mechanism
const int points_childIN= 2;
const int points_giftsIN= 2;
const int points_childOUT= -1;
const int points_giftsOUT= -1;

//#define ZERO_TEMP
//const double temperature= 0;

//const int weak_search = 1;  // Number of mini_step before forced search in the unhappy_child
//const int non_weak_search = 50; // Number of mini_step before search outside the unhappy_child list
const double output_threshold = 0.95;

// Time
const long long int step_end= 50000;
const long long int step_out= 20000;

const int max_runtime = 600;  // Max run time in sec [28800s = 8hrs]

const char name_wl_child[100] = "../input/child_wishlist_v2.csv";
const char name_wl_gifts[100] = "../input/gift_goodkids_v2.csv";
//const char name_initCF[100] = "../input/sample_submission_random_v2.csv";
const char name_initCF[100] = "../results/max_flow_mc_v11_optimum_01_936301288849572.csv";
//const char name_initCF[100] = "../results/result_v2_4.csv";
const char output_path[100] = "";  // path and prefix for output result

#endif // SANTA_PAR_INCLUDED
