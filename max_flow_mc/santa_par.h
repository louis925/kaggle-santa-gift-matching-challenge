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

// Precision of happiness points
typedef short HappyInt;
typedef long long int Eff_R;
//typedef long long int HappyInt;
//typedef double HappyInt;

const double keep_r = 0;                 // Ratio of successful swapped gifts using for next run

const int n_mf = 50;                   // Number of Max Flow runs in MC
const int n_gift_swap = 100;             // Number of gifts to be swap per MF run
const int n_child_swap = 100000;            // Number of child to be swap per MF run
										 // Note: This version (v5) is only efficient for small n_gift_swap (< 200) because of memory limit
//  8 sec for  72 gift swap ( 18 avg run/s 1.476 swap/s at max_fix = inf)
// 60 sec for 200 gift swap (134 avg run/s 1.484 swap/s at max_fix = 50)
const int step_output = 1;               // Number of steps before write improvement output
const int step_eq_output = 5;          // Number of steps before write same score output
const int max_eq_fixer_step = 2;         // Maximum step before stop fixing an same score inconsistent answer
const int max_fixer_step = 20;            // Maximum step before stop fixing an improvement inconsistent answer
const int time_limit = 30*60;          // Maximum run time in sec before stop next max flow run (excluding fixer)

//#define NO_FIX
#define FIX_REC                        // v9           - Fixer will fix recursively
//#define FIX_BEST                       // v7.1         - Fixer only fix the best improved triplet or twin and leaving other unfixed
//#define FIX_RAND                       // same as v7   - Fixer only fix one random triplet or twin and leaving other unfixed
//#define FIX_ALL                          // same as v5   - Fixer fix all the inconsistent triplets or twins
//#define MERGED_NO_FIX

//#define UPDATE_EFF_R_DURING_FIX        // same as v7.1 - Update eff_r even in the fixer step
//#define REBUILD_HAPPY // Rebuild happiness map each run to save memory but take more time
//#define NO_SINGLE_MAP // Do not put single child in the map
//#define NO_GIFT_MAP   // Do not include gift wishlist
//#define NO_INIT       // No initialization (load files) for debug only
//#define ZERO_GIFT_R   // whether to set contribution from gift happiness to zero
#define MERGED        // whether to merge twins and triplets or not
						
						// Initial effective ratio between child and gift happiness (r^3 C^2 / G^2) 
						// (Has no effect in this version v5)
						// Maximum eff_r = 10^6 to avoid overflow for int, 
						// Can change HappyInt to long long int if needed.
//#define FIX_EFF_R       // Use a fixed eff_r. Note this has no effect if you are using ZERO_GIFT_R
#if defined(FIX_EFF_R) && !defined(ZERO_EFF_R)
const Eff_R eff_r = 227585;
#endif // FIX_EFF_R

//#define VERBOSE
#define WINDOWS

// Files
const char name_wl_child[100] = "../input/child_wishlist_v2.csv";
const char name_wl_gifts[100] = "../input/gift_goodkids_v2.csv";
//const char name_initCF[100] = "../input/sample_submission_random_v2.csv";  // To start from very beginning
//const char name_initCF[100] = "../others/02_public_subm.csv";
const char name_initCF[100] = "../results/max_flow_merge_mc_optimum_14_936301288849572.csv";  // To start from previous saved result
const char output_path[100] = "../results/max_flow_mc_v11";      // path or output result

#endif // SANTA_PAR_INCLUDED
