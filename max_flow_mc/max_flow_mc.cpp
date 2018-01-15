#define USE_GLOP
#define _CRT_SECURE_NO_WARNINGS

#include <fstream>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm> 
#include <random>
#include <stack> 

//#include "ortools/linear_solver/linear_solver.h"
//#include "ortools/linear_solver/linear_solver.pb.h"
#include "ortools/graph/min_cost_flow.h"

#include "santa_par.h"
#include "santa_global.h"

using namespace std;

random_device rd;        //Will be used to obtain a seed for the random number engine
mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
						 //uniform_int_distribution<int> distribution(0, 6);
uniform_real_distribution<double> ran_dis(0.0, 1.0);

inline double ran_generator() {
	return ran_dis(generator);
	//return rand() / ((double)RAND_MAX + 1.0);
}
inline void get_twins(int ichild, int& twin) {
	if (ichild % 2 == 0) twin = ichild - 1;  // Note the twin now start at odd number id
	else                 twin = ichild + 1;
}
inline void get_triplets(int ichild, int& triplet1, int& triplet2) {
	switch (ichild % 3)
	{
	case 0:	triplet1 = ichild + 1; triplet2 = ichild + 2; break;
	case 1: triplet1 = ichild - 1; triplet2 = ichild + 1; break;
	case 2: triplet1 = ichild - 2; triplet2 = ichild - 1; break;
	default: break;
	}
}

void pick_gifts(const vector<int>& gift_swap_count, vector<int>& gift_swap, unordered_map<int, int>& gift_to_ig) {
	int n_prev_gift = 0;  // Number of previous successful swapped gifts
	for (int i = 0; i < n_gift_swap; i++) {  // Pick gift that has successful swap in previous step
		if (gift_swap_count[i] > 0) {
			gift_swap[n_prev_gift] = gift_swap[i];  // Collect the previous swapped gift to the front of the list
			n_prev_gift++;
		}
	}
	int max_prev_gift = (int)(n_gift_swap * keep_r);  // Keep at most only the last n_gift_swap * keep_r gift
	if (n_prev_gift > max_prev_gift) {
		for (int i = 0; i < max_prev_gift; i++) {
			gift_swap[i] = gift_swap[n_prev_gift - max_prev_gift + i];
		}
		n_prev_gift = max_prev_gift;
	}
	for (int i = n_prev_gift; i < gift_swap.size(); i++) {  // Pick random gift for the rest
		int g;
		vector<int>::iterator it;
		do {
			g = (int)(Ntypes * ran_generator());
			it = find(gift_swap.begin(), gift_swap.end(), g);  // Not to take the same gift including prevous searched gift
		} while (it != gift_swap.end());
		gift_swap[i] = g;
	}

	// Build map from giftId to index in gift_swap list
	for (int ig = 0; ig < gift_swap.size(); ig++)  gift_to_ig[gift_swap[ig]] = ig;
}

void pick_gifts_rand(vector<int>& gift_swap, unordered_map<int, int>& gift_to_ig) {
	for (int i = 0; i < gift_swap.size(); i++) gift_swap[i] = -1;

	for (int i = 0; i < gift_swap.size(); i++) {  // Pick random gift all the time
		int g;
		vector<int>::iterator it;
		do {
			g = (int)(Ntypes * ran_generator());
			it = find(gift_swap.begin(), gift_swap.end(), g);  // Not to take the same gift
		} while (it != gift_swap.end());
		gift_swap[i] = g;
	}

	// Build map from giftId to index in gift_swap list
	for (int ig = 0; ig < gift_swap.size(); ig++)  gift_to_ig[gift_swap[ig]] = ig;
}

namespace operations_research {
	void RunMaxFlow2(const vector<int>& gift_swap, const vector<int>& gift_qty, const vector<int>& child_swap, vector<int>& new_gift, 
		const vector<vector<HappyInt>>& happiness_c, const vector<vector<HappyInt>>& happiness_g) {
		// Given gift in the gift_swap list each with quantity gift_qty
		// Swap childs given in the child_swap
		// Return the new gift assignment in new_gift
		// child_swap mark as -1 won't be swapped and the new_gift won't be affect

		int n_gift_swap = gift_swap.size();

		// BUILD HAPPINESS MAP
		//vector<vector<HappyInt>> happiness(Nchild, vector<HappyInt>(Ntypes, 0));  // Map of happyness given gift for each child. happiness[child][gift] == approx_happiness 
		//build_happines_all(happiness);

		SimpleMinCostFlow max_flow;

		// ADD ARCS
		cout << "Adding arcs... ";
		for (int ic = 0; ic < child_swap.size(); ic++) {
			int c = child_swap[ic];
			if (c >= 0) {
				for (int ig = 0; ig < gift_swap.size(); ig++) {
					int g = gift_swap[ig];
#ifdef MERGED
					if (c >= Ntwins) {       // single
						max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 1, -happiness_c[c][g] * eff_r - happiness_g[c][g]);
					}
					else if (c >= Ntripl) {  // twin
						max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 2, -happiness_c[c][g] * eff_r - happiness_g[c][g]);  // count 2 times
					}
					else {                   // triplet
						max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 3, -happiness_c[c][g] * eff_r - happiness_g[c][g]);  // count 3 times
					}
#else
					max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 1, -happiness_c[c][g] * eff_r - (long long)happiness_g[c][g]);
#endif // MERGED	
				}
#ifdef MERGED
				// Skip the following twins and triplets (assuming they are in order)
				if (c < Ntripl) { ic += 1; }  // triplet
				if (c < Ntwins) { ic += 1; }  // twin
#endif // MERGED
#ifdef REBUILD_HAPPY  // Save some memory			
#ifdef MERGED
				if (c < Ntripl) { happiness_c[c + 2].~vector(); happiness_g[c + 2].~vector(); }  // triplet
				if (c < Ntwins) { happiness_c[c + 1].~vector(); happiness_g[c + 1].~vector(); }  // twin
#endif // MERGED
				happiness_c[c].~vector();  happiness_g[c].~vector();
#endif
			}
		}
#ifdef REBUILD_HAPPY
		happiness_c.~vector();  // Save some memory
		happiness_g.~vector();  // Save some memory
#endif

		// ASSIGN SUPPLY
		cout << "Assigning supply... " << endl;
		// Child
		for (int ic = 0; ic < child_swap.size(); ic++) {
			int c = child_swap[ic];
			if (c >= 0) {
#ifdef MERGED
				// Skip following twins and triplets
				if (c >= Ntwins) {       // single
					max_flow.SetNodeSupply(c, 1);
				}
				else if (c >= Ntripl) {  // twin
					max_flow.SetNodeSupply(c, 2); ic += 1;
				}
				else {                   // triplet
					max_flow.SetNodeSupply(c, 3); ic += 2;
				}
#else
				max_flow.SetNodeSupply(c, 1);  // every child counts as 1
#endif // MERGED	
			}
		}
		// Gift
		for (int ig = 0; ig < n_gift_swap; ig++) {
			max_flow.SetNodeSupply(Nchild + gift_swap[ig], -gift_qty[ig]);
		}
#ifdef VERBOSE
		cout << "N of arcs:                         " << max_flow.NumArcs() << endl;
		cout << "N of nodes:                        " << max_flow.NumNodes() << endl;
#endif // VERBOSE



		// SOLVE MAX FLOW MIN COST
		cout << "Start solving Max Flow with Min Cost...";
		int result_status = max_flow.SolveMaxFlowWithMinCost();
		cout << " Done!" << endl;
#ifdef VERBOSE
		cout << "Solver status:                     " << result_status << endl;
		cout << "Maximum flow:                      " << max_flow.MaximumFlow() << endl;
		cout << "Optimal cost:                      " << -max_flow.OptimalCost() / ((double)eff_r * max_child_happiness) << endl;
#endif // VERBOSE

		

		// OBTAIN RESULT
		//cout << "Obtaining result..." << endl;
		vector<int> c_to_ic(Nchild, -1);  // Map from childId to index in the child_swap list (non-swap child appears as -1)
		for (int ic = 0; ic < child_swap.size(); ic++) {
			int c = child_swap[ic];
			if (c >= 0) c_to_ic[c] = ic;
		}

#ifdef MERGED
		vector<int> n_gift_tri_twin(Ntwins, 0);  // Triplet and twin filling position shift in new_gift
												 // Note the number will only store at the head of the twin or triplet
#endif // MERGED

		int gift_total = 0;
		for (int i = 0; i < max_flow.NumArcs(); i++) {
			if (max_flow.Flow(i) > 0) {
				int c = max_flow.Tail(i);
				int new_g = max_flow.Head(i) - Nchild;
#ifdef MERGED
				if (c < Ntwins) {              // twin or triplet
					int ic = c_to_ic[c];
					ic += n_gift_tri_twin[c];  // shift the filling by the amount
					int f = max_flow.Flow(i);  // number of gift for the child from this arcs
					n_gift_tri_twin[c] += f;   // update the number of gift that twin or triplet will have
					gift_total += f;           // update total number of gift assigned by the solver
					while (f > 0) {
						new_gift[ic] = new_g;  // fill the new gift
						ic++;
						f--;
					}
				}
				else {                         // single
					new_gift[c_to_ic[c]] = new_g;
					gift_total++;
				}
#else  // NON-MERGED
				new_gift[c_to_ic[c]] = new_g;
				gift_total++;
#endif // MERGED
			}
		}
		cout << "gift_total:                        " << gift_total << endl;

		return;
	}
}  // namespace operations_research

/* Fixer need to fix the happiness two maps */
// Fixer
// If there is inconsistent  twin and triplet, take the best gift for them
void fixer(const vector<int>& incon_child, vector<int>& child_swap, const vector<int>& new_gift,
	const vector<vector<HappyInt>>& happiness_c, const vector<vector<HappyInt>>& happiness_g, vector<int>& answ_this, vector<int>& gift_qty,
	unordered_map<int, int>& gift_to_ig) {
	for (int iic = 0; iic < incon_child.size(); iic++) {
		int ic = incon_child[iic];
		int c = child_swap[ic];
		if (c < Ntripl) {
			// Triplets
			int g1 = new_gift[ic];
			int g2 = new_gift[ic + 1];
			int g3 = new_gift[ic + 2];
			Eff_R h1 = eff_r * happiness_c[c][g1] + happiness_g[c][g1];
			Eff_R h2 = eff_r * happiness_c[c][g2] + happiness_g[c][g2];
			Eff_R h3 = eff_r * happiness_c[c][g3] + happiness_g[c][g3];
			
			int g_select;
			if (h1 > h2) {
				if (h1 > h3) g_select = g1;
				else         g_select = g3;
			}
			else {
				if (h2 > h3) g_select = g2;
				else         g_select = g3;
			}

			answ_this[c] = g_select;
			answ_this[c + 1] = g_select;
			answ_this[c + 2] = g_select;

			child_swap[ic] = -1;
			child_swap[ic + 1] = -1;
			child_swap[ic + 2] = -1;
			gift_qty[gift_to_ig[g_select]] -= 3;  // Note there might not have enough gifts
		}
		else {
			// Twins
			int g1 = new_gift[ic];
			int g2 = new_gift[ic + 1];
			Eff_R h1 = eff_r * happiness_c[c][g1] + happiness_g[c][g1];
			Eff_R h2 = eff_r * happiness_c[c][g2] + happiness_g[c][g2];

			int g_select;
			if (h1 > h2) g_select = g1;
			else         g_select = g2;

			answ_this[c] = g_select;
			answ_this[c + 1] = g_select;

			child_swap[ic] = -1;
			child_swap[ic + 1] = -1;
			gift_qty[gift_to_ig[g_select]] -= 2;  // Note there might not have enough gifts
		}
	}
}

//// Randomly fix only one triplet or twin
//void fixer_rand_one(const vector<int>& incon_child, vector<int>& child_swap, const vector<int>& new_gift,
//	const vector<vector<HappyInt>>& happiness, vector<int>& answ_this, vector<int>& gift_qty,
//	unordered_map<int, int>& gift_to_ig) {
//	//for (int iic = 0; iic < incon_child.size(); iic++) {
//	int iic = incon_child.size() * ran_generator();
//	int ic = incon_child[iic];
//	int c = child_swap[ic];
//	if (c < Ntripl) {
//		// Triplets
//		int g1 = new_gift[ic];
//		int g2 = new_gift[ic + 1];
//		int g3 = new_gift[ic + 2];
//		HappyInt h1 = happiness[c][g1];
//		HappyInt h2 = happiness[c][g2];
//		HappyInt h3 = happiness[c][g3];
//
//		int g_select;
//		if (h1 > h2) {
//			if (h1 > h3) g_select = g1;
//			else         g_select = g3;
//		}
//		else {
//			if (h2 > h3) g_select = g2;
//			else         g_select = g3;
//		}
//
//		answ_this[c] = g_select;
//		answ_this[c + 1] = g_select;
//		answ_this[c + 2] = g_select;
//
//		child_swap[ic] = -1;
//		child_swap[ic + 1] = -1;
//		child_swap[ic + 2] = -1;
//		gift_qty[gift_to_ig[g_select]] -= 3;  // Note there might not have enough gifts
//	}
//	else {
//		// Twins
//		int g1 = new_gift[ic];
//		int g2 = new_gift[ic + 1];
//		HappyInt h1 = happiness[c][g1];
//		HappyInt h2 = happiness[c][g2];
//
//		int g_select;
//		if (h1 > h2) g_select = g1;
//		else         g_select = g2;
//
//		answ_this[c] = g_select;
//		answ_this[c + 1] = g_select;
//
//		child_swap[ic] = -1;
//		child_swap[ic + 1] = -1;
//		gift_qty[gift_to_ig[g_select]] -= 2;  // Note there might not have enough gifts
//	}	
//}
//
//// Fix only the best one triplet or twin
//void fixer_best_one(const vector<int>& incon_child, vector<int>& child_swap, const vector<int>& new_gift,
//	const vector<vector<HappyInt>>& happiness, vector<int>& answ_this, vector<int>& gift_qty,
//	unordered_map<int, int>& gift_to_ig) {
//
//	// find the best difference in each twin or triplet
//	vector<HappyInt> max_h_diff(incon_child.size());
//	for (int iic = 0; iic < incon_child.size(); iic++) {
//		int ic = incon_child[iic];
//		int c = child_swap[ic];
//		if (c < Ntripl) {
//			// Triplets
//			int g1 = new_gift[ic];
//			int g2 = new_gift[ic + 1];
//			int g3 = new_gift[ic + 2];
//			HappyInt h1 = happiness[c][g1];
//			HappyInt h2 = happiness[c][g2];
//			HappyInt h3 = happiness[c][g3];
//			HappyInt hmax = max(h1, h2);
//			hmax = max(hmax, h3);
//			max_h_diff[iic] = 3 * hmax - h1 - h2 - h3;
//		}
//		else {
//			// Twins
//			int g1 = new_gift[ic];
//			int g2 = new_gift[ic + 1];
//			HappyInt h1 = happiness[c][g1];
//			HappyInt h2 = happiness[c][g2];
//			HappyInt hmax = max(h1, h2);
//			max_h_diff[iic] = 2 * hmax - h1 - h2;
//		}
//	}
//
//	// find the best one among incon_child
//	int best_iic = 0; HappyInt best_h_diff = numeric_limits<long long>::min();
//	for (int iic = 0; iic < incon_child.size(); iic++) {
//		if (max_h_diff[iic] > best_h_diff) {
//			best_h_diff = max_h_diff[iic]; best_iic = iic;
//		}
//	}
//	int ic = incon_child[best_iic];
//	int c = child_swap[ic];
//	if (c < Ntripl) {
//		// Triplets
//		int g1 = new_gift[ic];
//		int g2 = new_gift[ic + 1];
//		int g3 = new_gift[ic + 2];
//		HappyInt h1 = happiness[c][g1];
//		HappyInt h2 = happiness[c][g2];
//		HappyInt h3 = happiness[c][g3];
//
//		int g_select;
//		if (h1 > h2) {
//			if (h1 > h3) g_select = g1;
//			else         g_select = g3;
//		}
//		else {
//			if (h2 > h3) g_select = g2;
//			else         g_select = g3;
//		}
//
//		answ_this[c] = g_select;
//		answ_this[c + 1] = g_select;
//		answ_this[c + 2] = g_select;
//
//		child_swap[ic] = -1;
//		child_swap[ic + 1] = -1;
//		child_swap[ic + 2] = -1;
//		gift_qty[gift_to_ig[g_select]] -= 3;  // Note there might not have enough gifts
//	}
//	else {
//		// Twins
//		int g1 = new_gift[ic];
//		int g2 = new_gift[ic + 1];
//		HappyInt h1 = happiness[c][g1];
//		HappyInt h2 = happiness[c][g2];
//
//		int g_select;
//		if (h1 > h2) g_select = g1;
//		else         g_select = g2;
//
//		answ_this[c] = g_select;
//		answ_this[c + 1] = g_select;
//
//		child_swap[ic] = -1;
//		child_swap[ic + 1] = -1;
//		gift_qty[gift_to_ig[g_select]] -= 2;  // Note there might not have enough gifts
//	}
//}

void pick_child_rand(vector<int>& child_swap, const vector<int>& answ_this, unordered_map<int, int>& gift_to_ig, vector<int>& gift_qty) {
	vector<int> child_rand_list(Nchild); for (int i = 0; i < Nchild; i++) { child_rand_list[i] = i; }
	shuffle(child_rand_list.begin(), child_rand_list.end(), generator);
	vector<bool> twin_triplet_exit(Ntwins, false);   // Mark if the twin or triplet exist in the child_swap list
	int fill_i = 0; int take_i = 0;
	while (fill_i < child_swap.size()) {
		int c = child_rand_list[take_i];
		int g = answ_this[c];  // the gift which the child currently has
		if (gift_to_ig.find(g) != gift_to_ig.end()) {  // take the child only in the gift list
			if (c < Ntripl) {       // Triplet
				if (child_swap.size() - fill_i >= 3) {   // if there is enought room for triplet
					int c2, c3;
					get_triplets(c, c2, c3);
					c = min(c, c2); c = min(c, c3);    // find the head of triplet
					if (!twin_triplet_exit[c]) {       // If haven't added this triplet to child_swap
						child_swap[fill_i] = c;
						child_swap[fill_i + 1] = c + 1;
						child_swap[fill_i + 2] = c + 2;
						twin_triplet_exit[c] = true;
						fill_i += 3;
						gift_qty[gift_to_ig[g]] += 3;
					}
				}
			}
			else if (c < Ntwins) {  // Twin
				if (child_swap.size() - fill_i >= 2) {  // if there is enought room for twin
					int c2;
					get_twins(c, c2);
					c = min(c, c2);                    // find the head of twin
					if (!twin_triplet_exit[c]) {       // If haven't added this triplet to child_swap
						child_swap[fill_i] = c;
						child_swap[fill_i + 1] = c + 1;
						twin_triplet_exit[c] = true;
						fill_i += 2;
						gift_qty[gift_to_ig[g]] += 2;
					}
				}
			}
			else {                  // Single
				child_swap[fill_i] = c;
				fill_i += 1;
				gift_qty[gift_to_ig[g]] += 1;
			}
		}
		take_i++;
	}
}

// Pick child without fix. Only work for full number of gift quantity
void pick_child_no_fix(vector<int>& child_swap, const vector<int>& answ_this, unordered_map<int, int>& gift_to_ig, vector<int>& gift_qty) {
	int i_fill = 0;
	for (int c = 0; c < Nchild; c++) {
		int g = answ_this[c];
		if (gift_to_ig.find(g) != gift_to_ig.end()) {
			child_swap[i_fill] = c;
			i_fill++;
		}
	}
	for (int i = 0; i < gift_qty.size(); i++) {
		gift_qty[i] = Ntypes;
	}
}

int main() {
	cout << "======================    MC Max Flow    ======================" << endl;
	int t0cpu = time(0);
#ifndef NO_INIT
	//double init_happy = init();                  // Read wishlists and build happiness maps
	init();                                      // Read wishlists and build happiness maps
	vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
	vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
	build_happines_all(happiness_c, happiness_g);

	long long int C_init = cal_child_all(answ, happiness_c);  // Initial total child happiness points
	long long int G_init = cal_gift_all(answ, happiness_g);   // Initial total gift happiness points
	double init_happy = cal_all_from_points(C_init, G_init);

	vector<int> answ_init(answ);                 // Initial answer
												 // Note answ is the best answ, initialized by init() from input data
	long long int C_best = C_init;               // Total child happiness points
	long long int G_best = G_init;               // Total gift happiness points
	double best_happy = init_happy;              // Previous best happiness

	printf( "Initial score:                     %.15f\n", init_happy);
	cout << "Child / Gifts happiness points:    " << C_init << " / " << G_init << endl;
	cout << "Child / Gifts ori happiness points:" << cal_child_all_ori(answ) << " / " << cal_gift_all_ori(answ) << endl;
	//cout << "Initial effective happiness ratio: " << eff_r_init;
	Eff_R eff_r_init = update_eff_r(C_init, G_init);  // Initial effective happiness ratio
	//cout << "Use effective happiness ratio:     " << eff_r;
#ifdef FIX_EFF_R
	cout << "Using fixed eff_r:                 " << eff_r << endl;
#endif // FIX_EFF_R
#ifdef ZERO_GIFT_R
	cout << "Initial gift_r = 0" << endl;
#endif // ZERO_GIFT_R

#ifdef REBUILD_HAPPY
	happiness_c.~vector();
	happiness_g.~vector();
#endif // !REBUILD_HAPPY
	
	cout << "Finished initialization." << endl;

	cout << "### Parameters" << endl;
	cout << "Max fix steps:                     " << max_fixer_step << endl;
	cout << "Time limit:                        " << time_limit << "s" << endl;
	cout << "Keep ratio:                        " << keep_r << endl;
#endif // !NO_INIT

	// ====== MC on Max Flow ======
	int t0mc = time(0);
	vector<int> gift_swap(n_gift_swap, -1);       // Gift you want to swap their assigned children
	vector<int> gift_swap_count(n_gift_swap, 0);  // N of child been swapped with gift in the gift_swap
	for (int step_mf = 1; step_mf <= n_mf; step_mf++) {
		cout << endl;
		cout << "#####  MF run " << step_mf << " / " << n_mf << "  #####"<< endl;
		cout << "n_gift_swap = " << n_gift_swap << endl;
		double time_run = -time(NULL);  // Runtime for the large loop run


		// Copy best answ to this answ
		vector<int> answ_this(answ);
		long long int C_this = C_best;
		long long int G_this = G_best;
		double this_happy = best_happy;

		
#ifndef UPDATE_EFF_R_DURING_FIX
		// Update eff_r and rebuild happiness map if they are different
		Eff_R eff_r_before = eff_r;
		update_eff_r(C_this, G_this);
#ifdef REBUILD_HAPPY
		vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
		vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
		build_happines_all(happiness_c, happiness_g);
//#else
//		if (eff_r_before != eff_r) {
//			reset_happiness_map(happiness_c, points_childOUT * single_r);
//			reset_happiness_map(happiness_g, points_giftsOUT * gift_r);
//			build_happines_all(happiness_c, happiness_g);
//		}
#endif // REBUILD_HAPPY
#endif // !UPDATE_EFF_R_DURING_FIX


		// Pick gifts for exchange
		unordered_map<int, int> gift_to_ig;                  // Map from giftId to index in gift_swap list
		pick_gifts_rand(gift_swap, gift_to_ig);              // Build gift_swap list and gift_to_ig map
		//pick_gifts(gift_swap_count, gift_swap, gift_to_ig);  // Build gift_swap list and gift_to_ig map
#ifdef VERBOSE
		cout << "gift_swap: ";
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			cout << gift_swap[ig] << " ";
		}
		cout << endl;
#endif // VERBOSE


		// Pick children for exchange
		vector<int> gift_qty(gift_swap.size(), 0);       // Default quantity of each gift
		vector<int> child_swap(n_child_swap, -1);        // Create child list for swapping
#if defined(NO_FIX) || defined(MERGED_NO_FIX)
		pick_child_no_fix(child_swap, answ_this, gift_to_ig, gift_qty);
#else
		pick_child_rand(child_swap, answ_this, gift_to_ig, gift_qty);
#endif
#ifdef VERBOSE
		for (int iic = 0; iic < child_swap.size(); iic++) {
			if (child_swap[iic] < 0) {
				cout << "no child swap mark at " << iic << endl;
			}
		}
#endif // VERBOSE
		// ************************************* //


		vector<int> new_gift(child_swap.size(), -1);    // The solver will return new gift assignment to this



		// FIXER LOOPS -------------------------------------------------------
		int n_inconsistent = 0;
		vector<int> incon_child;  // first index of the inconsistent twins or triplet in the child_swap
#ifdef FIX_REC
		stack<vector<int>> child_swap_stack;
		stack<vector<int>> answ_this_stack;
		stack<vector<int>> gift_qty_stack;
		stack<int> n_incon_stack;
		int prev_n_incon = numeric_limits<int>::max();
#endif // FIX_REC

		int fixer_step = 0;
		do {
#if !defined(NO_FIX) && !defined(MERGED_NO_FIX)
			cout << "=== Fixer step " << fixer_step << endl;

#ifdef UPDATE_EFF_R_DURING_FIX
			// Update eff_r and rebuild happiness map if they are different
			HappyInt eff_r_before = eff_r;
			update_eff_r(C_this, G_this);

#ifdef REBUILD_HAPPY
			vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
			vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
			build_happines_all(happiness_c, happiness_g);
//#else
//			if (eff_r_before != eff_r) {
//				reset_happiness_map(happiness_c, points_childOUT * single_r);
//				reset_happiness_map(happiness_g, points_giftsOUT * gift_r);
//				build_happines_all(happiness_c, happiness_g);
//			}
#endif // REBUILD_HAPPY
#endif // UPDATE_EFF_R_DURING_FIX

			// Fixer
			// If there is inconsistent  twin and triplet, take the best gift for them
#ifdef FIX_REC
				//void fixer_recursive(answ_this, gift_qty, new_gift);
				// push incon_child to incon_child_stack
				// push answ_this to the stack
				// push 
			cout << "stack size: " << child_swap_stack.size() << endl;
			if (!child_swap_stack.empty()) {
				child_swap = child_swap_stack.top();
				answ_this = answ_this_stack.top();
				gift_qty = gift_qty_stack.top();
				prev_n_incon = n_incon_stack.top();
				child_swap_stack.pop();
				answ_this_stack.pop();
				gift_qty_stack.pop();
				n_incon_stack.pop();
			}
#else
			if (incon_child.size() > 0) {
#ifdef FIX_BEST
				fixer_best_one(incon_child, child_swap, new_gift, happiness, answ_this, gift_qty, gift_to_ig);
#else 
#ifdef FIX_RAND
				fixer_rand_one(incon_child, child_swap, new_gift, happiness, answ_this, gift_qty, gift_to_ig);
#else
#ifdef FIX_ALL
				fixer(incon_child, child_swap, new_gift, happiness_c, happiness_g, answ_this, gift_qty, gift_to_ig);
#else
				cout << "No fixer!" << endl;
#endif
#endif // FIX_RAND
#endif // FIX_BEST
			}
#endif // FIX_REC
#endif // NO_FIX

			// Solve the Max Flow Min Cost for gift in the gift_swap list
			double time_solved = -time(NULL);  // Runtime used by max flow solver
			operations_research::RunMaxFlow2(gift_swap, gift_qty, child_swap, new_gift, happiness_c, happiness_g);                   // THE ACTUAL SOLVER
			time_solved += time(NULL);
			cout << "Max flow runtime:                  " << time_solved << " s" << endl;
			
			
#if !defined(NO_FIX) && !defined(MERGED_NO_FIX)
			// Check if twin and triplet are consistent
			n_inconsistent = new_gift_inconsistent(new_gift, child_swap, incon_child);
			cout << "N of inconsistent twins+triplets:  " << n_inconsistent << endl;			
#endif // !NO_FIX			

			// Write the new_gift to answ_this
			for (int ic = 0; ic < child_swap.size(); ic++) {
				int c = child_swap[ic];
				if (c >= 0)  answ_this[c] = new_gift[ic];
			}


			// Compute this total child and gift happiness points
			C_this = cal_child_all(answ_this, happiness_c);
			G_this = cal_gift_all(answ_this, happiness_g);
			printf("Child happiness points:            %17lld    [ %+lld ]\n", C_this, C_this - C_best);
			printf("Gift  happiness points:            %17lld    [ %+lld ]\n", G_this, G_this - G_best);
			this_happy = cal_all_from_points(C_this, G_this);
			printf("This score:                        %.15f    [ %+.15f ]\n", this_happy, this_happy - best_happy);

			
#ifdef VERBOSE  //DEBUG
			if (this_happy - best_happy < 0 && n_inconsistent == 0 && fixer_step == 0) {
				//verify_app_happy_points(answ_this, happiness, C_this, G_this);
				cout << "Why Worse?" << endl;
				check_gift_count(answ_this);				
			}
#endif // VERBOSE

#if defined(NO_FIX) || defined(MERGED_NO_FIX)
			break;  // No Fixer loop
#elif defined(FIX_REC)
			if (n_inconsistent != 0 && this_happy > best_happy && n_inconsistent <= prev_n_incon) {
				for (int iic = 0; iic < incon_child.size(); iic++) {
					int ic = incon_child[iic];
					int c = child_swap[ic];
					vector<int> temp_child_swap(child_swap);
					if (c < Ntripl) {
						// Triplets
						temp_child_swap[ic] = -1;
						temp_child_swap[ic + 1] = -1;
						temp_child_swap[ic + 2] = -1;

						for (int t = 0; t < 3; t++) {
							int g = new_gift[ic + t];
							vector<int> temp_answ_this(answ_this);
							temp_answ_this[c] = g;
							temp_answ_this[c + 1] = g;
							temp_answ_this[c + 2] = g;

							vector<int> temp_gift_qty(gift_qty);
							temp_gift_qty[gift_to_ig[g]] -= 3;  // Note there might not have enough gifts

							child_swap_stack.push(temp_child_swap);
							answ_this_stack.push(temp_answ_this);
							gift_qty_stack.push(temp_gift_qty);
							n_incon_stack.push(n_inconsistent);
						}
					}
					else {
						// Twins
						temp_child_swap[ic] = -1;
						temp_child_swap[ic + 1] = -1;

						for (int t = 0; t < 2; t++) {
							int g = new_gift[ic + t];
							vector<int> temp_answ_this(answ_this);
							temp_answ_this[c] = g;
							temp_answ_this[c + 1] = g;

							vector<int> temp_gift_qty(gift_qty);
							temp_gift_qty[gift_to_ig[g]] -= 2;  // Note there might not have enough gifts

							child_swap_stack.push(temp_child_swap);
							answ_this_stack.push(temp_answ_this);
							gift_qty_stack.push(temp_gift_qty);
							n_incon_stack.push(n_inconsistent);
						}
					}
				}
			}
			if (child_swap_stack.empty() || //fixer_step >= max_fixer_step || 
				(n_inconsistent == 0 && this_happy > best_happy)) break;
#else
			if (n_inconsistent == 0 ||
				(this_happy == best_happy && fixer_step >= max_eq_fixer_step) ||
				this_happy < best_happy || 
				fixer_step >= max_fixer_step) {
				break;
			}
#endif // FIX_REC
			fixer_step += 1;
		} while (true);
		// END of FIXER LOOP -------------------------------------------------------------
		time_run += time(NULL);
		cout << endl;


		// Check difference and efficiency, and return N of child in swapped gift
		for (int ig = 0; ig < gift_swap_count.size(); ig++)  gift_swap_count[ig] = 0;  // Clear previous swap count
		int n_diff = 0;  // N of swap child
		for (int c = 0; c < Nchild; c++) {
			if (answ_this[c] != answ[c]) {
				gift_swap_count[gift_to_ig[answ[c]]]++;  // Add swap count of gift
				n_diff++;
			}
		}
		printf( "N of swap child:                   %d    [ %.6f/s ]\n", n_diff, n_diff / time_run);

		// Update if there is differnce
		if (n_diff > 0) {
			// Check result
			cout << "N of inconsistent triplets, twins: " << triplet_consistent(answ_this) << ", "
				<< twin_consistent(answ_this) << endl;
			check_gift_count(answ_this);
			//double this_happy = cal_all(answ_this);
			//C_this = cal_child_all(answ_this);
			//G_this = cal_gift_all(answ_this);
			printf("Score:                             %.15f    [ %+.15f ]\n", this_happy, this_happy - best_happy);
			cout << "Child / Gifts happiness points:    " << C_this << " / " << G_this << endl;

			// Check if improved and is consistent
			if (this_happy >= best_happy && n_inconsistent == 0) {
				// Update best answ and scores
				C_best = C_this; G_best = G_this;
				for (int c = 0; c < answ.size(); c++)  answ[c] = answ_this[c];
				answ_to_states(answ_this);  // update best states from answ

				// Output optimization result
				if (this_happy > best_happy) {
					best_happy = this_happy;
					if (step_mf % step_output == 0) {
						cout << "Output improved result..." << endl;
						write_conf(this_happy, false);
					}
				}
				else if (step_mf % step_eq_output == 0) {
					cout << "Output record result..." << endl;
					write_conf(this_happy, false);
				}
			}
			else {
				// Not updating the best answ
				if (n_inconsistent != 0) cout << "Inconsistent, not using this result..." << endl;
				else                     cout << "Worse, not using this result..." << endl;
			}
		}  // n_diff > 0

		if (time(0) - t0mc > time_limit) break;
	}
	int tfmc = time(0);
	

	// Check Final Result
	update_eff_r(C_best, G_best);  // Update eff_r to the best one

	cout << endl;
	cout << "### Summary" << endl;
	printf("Final score:                       %.15f    [ %+.15f ]\n", best_happy, best_happy - init_happy);
	printf("Child happiness points:            %17lld    [ %+lld ]\n", C_best, C_best - C_init);
	printf("Gift  happiness points:            %17lld    [ %+lld ]\n", G_best, G_best - G_init);
	printf("Recommend effective ratio:         %17lld    [ %+lld ]\n", eff_r, eff_r - eff_r_init);
	int n_diff = compare_answ_diff(answ, answ_init);  // N of swap child
	printf("N of swap child:                   %17d    [ %.6f/s ]\n", n_diff, ((double) n_diff) / (tfmc - t0mc + 1e-3));

	// Output
	cout << "Writing result...";
	if (n_mf % step_output != 0) {
		write_conf(best_happy, false);
	}
	write_conf(best_happy, true);
	cout << " Done!" << endl;

	// Remind you the parameters
	cout << "### Remind you the parameters" << endl;
	cout << "n_gift_swap:                       " << n_gift_swap << endl;
	cout << "Max fix steps:                     " << max_fixer_step << endl;
	cout << "Time limit:                        " << time_limit << "s" << endl;
	cout << "Keep ratio:                        " << keep_r << endl;
	cout << "Total runtime:                     " << time(0) - t0cpu << "s" << endl;
	cout << "===================== End of MC Max Flow ======================" << endl;
#ifdef WINDOWS
	cin.get();
#endif // WINDOWS
	return 0;
}