#ifndef SANTA_PAR_INCLUDED
#define SANTA_PAR_INCLUDED

// System
const unsigned int Nchild=   1000000; // N of children
const unsigned int Ntypes=      1000; // N of different types of gifts
const unsigned int Nqttys=      1000; // quantity of one gifts
const unsigned int Ntripl=      5001; // N of triplet children
const unsigned int Ntwins=     45001; // N of twin and triplet children
const unsigned int Nchprf=       100; // length of child wish list
const unsigned int NchprfCut= Nchprf/5; // cutoff length of child wish list
const unsigned int Ngfprf=      1000; // length of gifts wish list

// Mechanism
const unsigned int points_childIN= 2;
const unsigned int points_giftsIN= 2;
const int points_childOUT= -1;
const int points_giftsOUT= -1;

typedef long long int HappyInt;  // Precision of happiness points

const HappyInt eff_r_init = 471351034;   // Initial effective ratio between child and 
										 // gift happiness (r^3 C^2 / G^2) 
										 // Maximum eff_r = 10^6 to avoid overflow for int, 
										 // Can change HappyInt to long long int if needed.
const double keep_r = 0;                 // Ratio of successful swapped gifts using for next run
const int n_mip = 10;                   // Number of MIP runs in MC
const int n_gift_swap = 40;              // Number of gifts to be swap per MIP run
// 40 sec for 10 gift swap
const int step_output = 1;              // Number of steps before write output and update eff_r

//#define NO_SINGLE_MAP // Do not put single child in the map
//#define NO_GIFT_MAP   // Do not include gift wishlist
//#define NO_INIT       // No initialization (load files) for debug only
//#define ZERO_GIFT_R   // whether to set contribution from gift happiness to zero
#define MERGED          // whether to merge twins and triplets or not
                        // (non-merge method haven't finished yet)

// Files
const char name_wl_child[100] = "../input/child_wishlist_v2.csv";
const char name_wl_gifts[100] = "../input/gift_goodkids_v2.csv";
//const char name_initCF[100] = "../input/sample_submission_random_v2.csv";  // To start from very beginning
//const char name_initCF[100] = "../others/02_public_subm.csv";
//const char name_initCF[100] = "../results/result_v2_6.csv";  // To start from previous saved result
const char name_initCF[100] = "../results/max_flow_merge_mc_optimum_07_93630118.csv";  // To start from previous saved result
const char output_path[100] = "../results/mip_mc_2";      // path or output result

#endif // SANTA_PAR_INCLUDED
