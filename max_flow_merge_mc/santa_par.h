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


const int n_mf = 300;                   // Number of Max Flow runs in MC
//#define CONST_N_GIFT_SWAP
#ifdef CONST_N_GIFT_SWAP
const double keep_r = 0;                  // Ratio of successful swapped gifts using for next run
const int n_gift_swap = 400;              // Number of gifts to be swap per MIP run
										   // Note: This version is only efficient for small n_gift_swap (< 200) because of memory limit
#else  // Variable n_gift_swap
// Number of gifts to be swap per run for each group_size (only group_size = 1,2,3 and 6 have been implemented.)
const int n_gift_swap_max = 1000; //460  650
//const int n_gift_swap_list[6] = { n_gift_swap_max, n_gift_swap_max, n_gift_swap_max, 0, 0, n_gift_swap_max };
//const int n_gift_swap_list[6] = { 0.5773 * n_gift_swap_max, 0.8165 * n_gift_swap_max, n_gift_swap_max, 0, 0, n_gift_swap_max };
const int n_gift_swap_list[6] = { 600, 400, 500, 0, 0, 600 };
#endif // CONST_N_GIFT_SWAP
const int n_child_swap = 190000;            // Number of child to be swap per MF run for group_size = 1

//242    sec/run for 400 gift (1)           grouping, 26.1 swap/s
// 85.7  sec/run for 460 gift (1*4 + 2 + 3) grouping, 20.2 swap/s
//199    sec/run for 1000 gift 120000 (1)   grouping, 25.6 swap/s, 10.1GB
const int step_output = 1;               // Number of steps before write improvement output
const int step_eq_output = 6;            // Number of steps before write same score output
//const int max_eq_fixer_step = 2;         // Maximum step before stop fixing an same score inconsistent answer
//const int max_fixer_step = 200;          // Maximum step before stop fixing an improvement inconsistent answer
const int time_limit = 5*60*60;            // Maximum run time in sec before stop next max flow run (excluding fixer)

//#define FIX_GROUPSIZE                  // If you don't want to search in all group_size
#ifdef FIX_GROUPSIZE
const int group_size_fixed = 1;        // the group_size you want to use
#else
const int group_size_list[] = {1,1,1,1,1,1,2,3,6};  // Running order of group size
const int group_size_list_len = 9;
#endif // FIX_GROUPSIZE

#define REBUILD_HAPPY // Rebuild happiness map each run to save memory but take more time
//#define NO_SINGLE_MAP // Do not put single child in the map
//#define NO_GIFT_MAP   // Do not include gift wishlist
//#define NO_INIT       // No initialization (load files) for debug only
//#define ZERO_GIFT_R   // whether to set contribution from gift happiness to zero
						
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
const char name_initCF[100] = "../results/max_flow_mc_v10.csv";  // To start from previous saved result
const char output_path[100] = "../results/max_flow_merge_mc_pc";      // path or output result

#endif // SANTA_PAR_INCLUDED
