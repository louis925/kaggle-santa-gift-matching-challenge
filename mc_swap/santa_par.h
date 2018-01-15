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

#define ZERO_TEMP
const double temperature= 0;
// const long int eff_r = child_gifts_happiness_ratio; // For initial few runs
const long long int eff_r = 478314120; // effective happiness ratio r^3 C^2 / G^2, need to adjust for later stage
const int match_times = 1000; // maximum number of search times to find another child to match the condition

//#define WEAK_SEARCH_MAIN      // More weak search
#ifdef WEAK_SEARCH_MAIN
const int non_weak_search = 2; // Number of mini_step before search outside the unhappy_child list
#else
const int weak_search = 10;  // Number of mini_step before forced search in the unhappy_child
#endif
const double unhappy_threshold = 0.3*eff_r*Nchprf*points_childIN;  // Below which counts as unhappy
const double output_threshold = 0.91;

const int pref_gift_search = 1;  // Number of mini_step before forced pick gift 2 only from child 1's prefer gift
//#define APP_DIFF  // Use first order approximate for cal_diff()

// Time
const int step_end= 28000;
const int step_out= 1;

const int max_runtime = 30*60;  // Max run time in sec [28800s = 8hrs]

const char name_wl_child[100]=  "../input/child_wishlist_v2.csv";
const char name_wl_gifts[100]=  "../input/gift_goodkids_v2.csv";
//const char name_initCF[100] = "../input/sample_submission_random_v2.csv";  // To start from very beginning
//const char name_initCF[100] = "../others/02_public_subm.csv";
const char name_initCF[100] = "../results/result_v2_6.csv";  // To start from previous saved result
const char output_path[100] = "../results/result_v2_7";  // path and prefix for output result

#endif // SANTA_PAR_INCLUDED
