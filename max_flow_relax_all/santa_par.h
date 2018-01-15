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

// Initial effective ratio between child and gift happiness (r^3 C^2 / G^2) 
const Eff_R eff_r = 483279342;
const Eff_R twin_shift = -1E5;
const Eff_R triplet_shift = -1E5;

#define VERBOSE
#define WINDOWS

// Files
const char name_wl_child[100] = "../input/child_wishlist_v2.csv";
const char name_wl_gifts[100] = "../input/gift_goodkids_v2.csv";
const char output_path[100] = "../results/max_flow_relax_all_-5";      // path or output result

#endif // SANTA_PAR_INCLUDED
