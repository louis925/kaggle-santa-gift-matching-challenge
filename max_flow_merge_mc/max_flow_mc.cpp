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
#include <numeric>

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

inline bool is_single(int ichild) { return ichild >= Ntwins; }
inline bool is_not_single(int ichild) { return ichild < Ntwins; }
inline bool is_twin(int ichild) { return ichild < Ntwins && ichild >= Ntripl; }
inline bool is_not_twin(int ichild) { return ichild >= Ntwins || ichild < Ntripl; }
inline bool is_triplet(int ichild) { return ichild < Ntripl; }
inline bool is_not_triplet(int ichild) { return ichild >= Ntripl; }

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

#ifdef CONST_N_GIFT_SWAP
void pick_gifts(const vector<int>& gift_swap_count, vector<int>& gift_swap, unordered_map<int, int>& gift_to_ig) {
	int n_prev_gift = 0;  // Number of previous successful swapped gifts
	for (int i = 0; i < gift_swap.size(); i++) {  // Pick gift that has successful swap in previous step
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
#endif // CONST_N_GIFT_SWAP

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

// sort vector of vector only by the first element from large to small
bool happy_sort(const vector<Eff_R>& a, const vector<Eff_R>& b) { return a[0] > b[0]; }

void build_child_swap_sort(int group_size, const vector<int>& gift_swap, 
	vector<int>& gift_qty,	vector<vector<int>>& child_swap, 
	vector<vector<HappyInt>>& happiness_c, vector<vector<HappyInt>>& happiness_g) {

	if (group_size == 6) {
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			int n_child_swap_init = child_swap.size();
			vector<vector<int>> triplets;  // list of triplets who have gift g
			vector<vector<int>> twins;     // list of twins who have gift g
			vector<int> singles;           // list of singles who has gift g
			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // extract triplets
				int c = states[g][ic];
				vector<int> triplet = { c, c + 1, c + 2 };
				triplets.push_back(triplet);
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // extract twins
				int c = states[g][ic];
				vector<int> pair = { c, c + 1 };
				twins.push_back(pair);
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				singles.push_back(c);
				ic += 1;
			}

			// shuffle
			shuffle(triplets.begin(), triplets.end(), generator);
			shuffle(twins.begin(), twins.end(), generator);
			shuffle(singles.begin(), singles.end(), generator);

			int i_single = 0; int i_tri;
			for (i_tri = 0; i_tri + 1 < triplets.size(); i_tri += 2) {       // pair all the triplets
				vector<int> row = { triplets[i_tri][0], triplets[i_tri][1], triplets[i_tri][2],
					triplets[i_tri + 1][0], triplets[i_tri + 1][1], triplets[i_tri + 1][2] };
				child_swap.push_back(row);
			}
			if (i_tri < triplets.size()) {  // leftover triplet
				vector<int> row = { triplets[i_tri][0], triplets[i_tri][1], triplets[i_tri][2] };
				while (row.size() < group_size) {
					row.push_back(singles[i_single]);
					i_single++;
				}
				child_swap.push_back(row);
			}

			int i_twin;
			for (i_twin = 0; i_twin + 2 < twins.size(); i_twin += 3) {    // 3 twins into 1 group
				vector<int> row = { twins[i_twin][0], twins[i_twin][1], twins[i_twin + 1][0],
					twins[i_twin + 1][1], twins[i_twin + 2][0], twins[i_twin + 2][1] };
				child_swap.push_back(row);
			}
			if (i_twin < twins.size()) {  // leftover twins
				vector<int> row;
				while (i_twin < twins.size()) {
					row.push_back(twins[i_twin][0]);
					row.push_back(twins[i_twin][1]);
					i_twin++;
				}
				while (row.size() < group_size) {
					row.push_back(singles[i_single]);
					i_single++;
				}
				child_swap.push_back(row);
			}

			// put remaining 6 singles in 1 group
			for (; i_single + 5 < singles.size(); i_single += 6) {
				vector<int> gp_c = { singles[i_single + 0], singles[i_single + 1], singles[i_single + 2],
					singles[i_single + 3], singles[i_single + 4], singles[i_single + 5] };
				child_swap.push_back(gp_c);
			}

			gift_qty[ig] = group_size * (child_swap.size() - n_child_swap_init);
		}
	}
	else if (group_size == 3) {
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			vector<vector<int>> twins(0);  // list of twins who have gift g
			vector<int> singles;        // list of singles who has gift g
			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // push triplets
				int c = states[g][ic];
				vector<int> gp_c = { c, c + 1, c + 2 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 3;
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // extract twins
				int c = states[g][ic];
				vector<int> pair = { c, c + 1 };
				twins.push_back(pair);
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				singles.push_back(c);
				ic += 1;
			}

			// shuffle
			shuffle(twins.begin(), twins.end(), generator);
			shuffle(singles.begin(), singles.end(), generator);

			int n_single_twin_gp = min(singles.size(), twins.size());  // number of group with 1 twin pairs and 1 single
																	   // put 1 twin pair with 1 single in 1 group
			for (int i = 0; i < n_single_twin_gp; i++) {
				int c1 = twins[i][0];
				int c2 = twins[i][1];
				int c3 = singles[i];
				vector<int> row = { c1, c2, c3 };
				child_swap.push_back(row);
			}
			gift_qty[ig] += 3 * n_single_twin_gp;

			// put 3 singles in 1 group
			for (int i = n_single_twin_gp; i + 2 < singles.size(); i += 3) {
				int c1 = singles[i];
				int c2 = singles[i + 1];
				int c3 = singles[i + 2];
				vector<int> gp_c = { c1, c2, c3 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 3;
				// there might be left over so N of 3 singles != singles.size() / 3
			}
		}
	}
	else if (group_size == 2) {
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];

			vector<vector<Eff_R>> singles_h;        // list of singles who has gift g
			int ic = 0;
			while (ic + 2 < Nqttys && is_triplet(states[g][ic])) {  // skip triplets
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // push twins
				int c = states[g][ic];
				vector<int> gp_c = { c, c + 1 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 2;
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * happiness_c[c][g] + happiness_g[c][g], c };
				singles_h.push_back(c_h);
				ic += 1;
			}

			// Sort singles by happiness from large to small
			sort(singles_h.begin(), singles_h.end(), happy_sort);

			// put 2 singles in 1 group
			for (int i = 0; i + 1 < singles_h.size(); i += 2) {
				int c1 = singles_h[i][1];
				int c2 = singles_h[i + 1][1];
				vector<int> gp_c = { c1, c2 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 2;
				// there might be left over so N of 3 singles != singles.size() / 3
			}
		}
	}
	else { // group == 1
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];

			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // skip triplets
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // skip twins
				ic += 2;
			}
			while (ic < Nqttys) {                               // push singles
				int c = states[g][ic];
				vector<int> gp_c = { c };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 1;
				ic += 1;
			}
		}
	}
}

// Same as build_child_swap_sort but improve the score on g1
void build_child_swap_sort_enhanced(int group_size, const vector<int>& gift_swap,
	vector<int>& gift_qty, vector<vector<int>>& child_swap,
	vector<vector<HappyInt>>& happiness_c, vector<vector<HappyInt>>& happiness_g, int ig1_shift) {

	static int include_twin = 0;
	static int include_twin_6 = 0;
	
	if (group_size == 6) {
		cout << "ig1_shift: " << ig1_shift << endl;
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			int g1 = gift_swap[(ig + ig1_shift) % gift_swap.size()];

			int n_child_swap_init = child_swap.size();
			vector<vector<Eff_R>> triplets_h;  // list of triplets who has gift g
			vector<vector<Eff_R>> twins_h;     // list of twins who has gift g
			vector<vector<Eff_R>> singles_h;   // list of singles who has gift g
			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // extract triplets
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g])
					+ (happiness_g[c][g1] - happiness_g[c][g]), c, c + 1, c + 2 };
				triplets_h.push_back(c_h);
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // extract twins
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g])
					+ (happiness_g[c][g1] - happiness_g[c][g]), c, c + 1 };
				twins_h.push_back(c_h);
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g])
					+ (happiness_g[c][g1] - happiness_g[c][g]), c };
				singles_h.push_back(c_h);
				ic += 1;
			}

			// Sort all by happiness from large to small
			sort(triplets_h.begin(), triplets_h.end(), happy_sort);
			sort(twins_h.begin(), twins_h.end(), happy_sort);
			sort(singles_h.begin(), singles_h.end(), happy_sort);

			int i_single = 0; int i_tri;
			for (i_tri = 0; i_tri + 1 < triplets_h.size(); i_tri += 2) {       // pair all the triplets
				vector<int> row = { (int)triplets_h[i_tri][1], (int)triplets_h[i_tri][2], (int)triplets_h[i_tri][3],
					(int)triplets_h[i_tri + 1][1], (int)triplets_h[i_tri + 1][2], (int)triplets_h[i_tri + 1][3] };
				child_swap.push_back(row);
			}
			//if (include_twin_6 % 2 == 1) {
			//	if (i_tri < triplets_h.size()) {  // leftover triplet
			//		vector<int> row = { (int)triplets_h[i_tri][1], (int)triplets_h[i_tri][2], (int)triplets_h[i_tri][3] };
			//		while (row.size() < group_size) {
			//			row.push_back(singles_h[i_single][1]);
			//			i_single++;
			//		}
			//		child_swap.push_back(row);
			//	}
			//}

			int i_twin;
			for (i_twin = 0; i_twin + 2 < twins_h.size(); i_twin += 3) {    // 3 twins into 1 group
				vector<int> row = { (int)twins_h[i_twin][1], (int)twins_h[i_twin][2], (int)twins_h[i_twin + 1][1],
					(int)twins_h[i_twin + 1][2], (int)twins_h[i_twin + 2][1], (int)twins_h[i_twin + 2][2] };
				child_swap.push_back(row);
			}
			if (i_twin < twins_h.size()) {  // leftover twins
				vector<int> row;
				while (i_twin < twins_h.size()) {
					row.push_back(twins_h[i_twin][1]);
					row.push_back(twins_h[i_twin][2]);
					i_twin++;
				}
				while (row.size() < group_size) {
					row.push_back(singles_h[i_single][1]);
					i_single++;
				}
				child_swap.push_back(row);
			}

			// put remaining 6 singles in 1 group
			for (; i_single + 5 < singles_h.size(); i_single += 6) {
				vector<int> gp_c = { (int)singles_h[i_single + 0][1], (int)singles_h[i_single + 1][1], (int)singles_h[i_single + 2][1],
					(int)singles_h[i_single + 3][1], (int)singles_h[i_single + 4][1], (int)singles_h[i_single + 5][1] };
				child_swap.push_back(gp_c);
			}

			gift_qty[ig] = group_size * (child_swap.size() - n_child_swap_init);
			include_twin_6++;
		}
	}
	else if (group_size == 3) {
		cout << "ig1_shift: " << ig1_shift << "  include twin: " << include_twin << endl;
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			int g1 = gift_swap[(ig + ig1_shift) % gift_swap.size()];

			vector<vector<Eff_R>> singles_h;        // list of singles who has gift g
			vector<vector<Eff_R>> twins_h;        // list of singles who has gift g
			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // push triplets
				int c = states[g][ic];
				vector<int> gp_c = { c, c + 1, c + 2 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 3;
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // extract twins
				if (include_twin % 2 == 1) {
					int c = states[g][ic];
					vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g])
						+ (happiness_g[c][g1] - happiness_g[c][g]), c, c + 1 };
					twins_h.push_back(c_h);
				}
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g])
					+ (happiness_g[c][g1] - happiness_g[c][g]), c };
				singles_h.push_back(c_h);
				ic += 1;
			}

			// Sort singles by happiness from large to small
			if (include_twin % 2 == 1) {
				sort(twins_h.begin(), twins_h.end(), happy_sort);
			}
			sort(singles_h.begin(), singles_h.end(), happy_sort);

			
			int n_single_twin_gp = 0;
			if (include_twin % 2 == 1) {
				n_single_twin_gp = min(singles_h.size(), twins_h.size());  // number of group with 1 twin pairs and 1 single
																			   // put 1 twin pair with 1 single in 1 group
				for (int i = 0; i < n_single_twin_gp; i++) {
					int c1 = twins_h[i][1];
					int c2 = twins_h[i][2];
					int c3 = singles_h[i][1];
					vector<int> row = { c1, c2, c3 };
					child_swap.push_back(row);
				}
				gift_qty[ig] += 3 * n_single_twin_gp;
			}
			// put 3 singles in 1 group
			for (int i = n_single_twin_gp; i + 2 < singles_h.size(); i += 3) {
				int c1 = singles_h[i][1];
				int c2 = singles_h[i + 1][1];
				int c3 = singles_h[i + 2][1];
				vector<int> gp_c = { c1, c2, c3 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 3;
				// there might be left over so N of 3 singles != singles.size() / 3
			}

			include_twin++;
		}
	}
	else if (group_size == 2) {
		cout << "ig1_shift: " << ig1_shift << endl;
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			int g1 = gift_swap[(ig + ig1_shift) % gift_swap.size()];

			vector<vector<Eff_R>> singles_h;        // list of singles who has gift g
			int ic = 0;
			while (ic + 2 < Nqttys && is_triplet(states[g][ic])) {  // skip triplets
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // push twins
				int c = states[g][ic];
				vector<int> gp_c = { c, c + 1 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 2;
				ic += 2;
			}
			while (ic < Nqttys) {                               // extract singles
				int c = states[g][ic];
				vector<Eff_R> c_h = { eff_r * (happiness_c[c][g1] - happiness_c[c][g]) 
					+ (happiness_g[c][g1] - happiness_g[c][g]), c };
				singles_h.push_back(c_h);
				ic += 1;
			}

			// Sort singles by happiness from large to small
			sort(singles_h.begin(), singles_h.end(), happy_sort);

			// put 2 singles in 1 group
			for (int i = 0; i + 1 < singles_h.size(); i += 2) {
				int c1 = singles_h[i][1];
				int c2 = singles_h[i + 1][1];
				vector<int> gp_c = { c1, c2 };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 2;
				// there might be left over so N of 3 singles != singles.size() / 3
			}
		}
	}
	else { // group == 1
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];

			int ic = 0;
			while (ic < Nqttys && is_triplet(states[g][ic])) {  // skip triplets
				ic += 3;
			}
			while (ic < Nqttys && is_twin(states[g][ic])) {     // skip twins
				ic += 2;
			}
			while (ic < Nqttys) {                               // push singles
				int c = states[g][ic];
				vector<int> gp_c = { c };
				child_swap.push_back(gp_c);
				gift_qty[ig] += 1;
				ic += 1;
			}
		}
	}
}

// Bhild child swap list with only single and random number in each group
void build_child_swap_single_rand(int group_size, const vector<int>& gift_swap,
	vector<int>& gift_qty, vector<vector<int>>& child_swap, 
	const vector<int>& answ_this, unordered_map<int, int>& gift_to_ig) {
	if (group_size != 1) error("Wrong build child swap on", 1, group_size);
	
	// Make shuffled child index list for singles
	vector<int> child_rand_list(Nchild - Ntwins); 
	for (int i = 0; i < child_rand_list.size(); i++) { child_rand_list[i] = i + Ntwins; }
	shuffle(child_rand_list.begin(), child_rand_list.end(), generator);
	
	int take_i = 0;
	while (child_swap.size() < n_child_swap && take_i < child_rand_list.size()) {
		int c = child_rand_list[take_i];
		int g = answ_this[c];  // the gift which the child currently has
		if (gift_to_ig.find(g) != gift_to_ig.end()) {  // take the child only in the gift list
			                // Single
			vector<int> gp_c = { c };
			child_swap.push_back(gp_c);
			gift_qty[gift_to_ig[g]] += 1;			
		}
		take_i++;
	}
}

// Not using this in merge_mc
void pick_child_rand(vector<int>& child_swap, const vector<int>& answ_this, unordered_map<int, int>& gift_to_ig, vector<int>& gift_qty) {
	vector<int> child_rand_list(Nchild); for (int i = 0; i < Nchild; i++) { child_rand_list[i] = i; }
	shuffle(child_rand_list.begin(), child_rand_list.end(), generator);
	vector<bool> twin_triplet_exist(Ntwins, false);   // Mark if the twin or triplet exist in the child_swap list
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
					if (!twin_triplet_exist[c]) {       // If haven't added this triplet to child_swap
						child_swap[fill_i] = c;
						child_swap[fill_i + 1] = c + 1;
						child_swap[fill_i + 2] = c + 2;
						twin_triplet_exist[c] = true;
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
					if (!twin_triplet_exist[c]) {       // If haven't added this triplet to child_swap
						child_swap[fill_i] = c;
						child_swap[fill_i + 1] = c + 1;
						twin_triplet_exist[c] = true;
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

namespace operations_research {
	// Given gift in the gift_swap list each with quantity gift_qty
	// Swap childs given in the child_swap
	// Return the new gift assignment in new_gift
	// child_swap: list of group of children [2D vector]
	// group     : number of children in 1 group (support 1, 2, 3)
	void RunMaxFlow3(const vector<int>& gift_swap, const vector<int>& gift_qty, 
		vector<vector<int>>& child_swap, vector<int>& new_gift, 
		const vector<vector<HappyInt>>& happiness_c, 
		const vector<vector<HappyInt>>& happiness_g, int group_size) {

		int n_gift_swap = gift_swap.size();
		int n_child_swap = child_swap.size();  // number of groups

		SimpleMinCostFlow max_flow;

		// ADD ARCS
		cout << "Adding arcs... ";
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			int g = gift_swap[ig];
			for (int igp = 0; igp < child_swap.size(); igp++) {
				Eff_R h = 0; //Eff_R hc = 0; Eff_R hg = 0;
				for (int ic = 0; ic < group_size; ic++) {
					int c = child_swap[igp][ic];
					//hc += happiness_c[c][g];
					//hg += happiness_g[c][g];  // combine happiness of children in each group
					h += eff_r * happiness_c[c][g] + happiness_g[c][g];  // combine happiness of children in each group
				}
				//h = eff_r * hc + hg;
				max_flow.AddArcWithCapacityAndUnitCost(igp, ig + n_child_swap, 1, -h);
				// AddArcWithCapacityAndUnitCost(tail, head, capicity, cost)
				// tail = index of group in the child_swap list
				// head = index of gift in the gift_swap list + N of child group
			}
		}
#ifdef REBUILD_HAPPY
		happiness_c.~vector();  // Save some memory
		happiness_g.~vector();  // Save some memory
#endif

							  // ASSIGN SUPPLY
		cout << "Assigning supply... " << endl;
		// Child group
		for (int igp = 0; igp < child_swap.size(); igp++) {
#ifdef MERGED
			// Skip following twins and triplets
			if (c >= Ntwins) {       // single
				max_flow.SetNodeSupply(c, 1);
			}
			else if (c >= Ntripl) {  // twin
				max_flow.SetNodeSupply(c, 2);
			}
			else {                   // triplet
				max_flow.SetNodeSupply(c, 3);
				}
#else
			max_flow.SetNodeSupply(igp, 1);  // every child counts as 1
#endif // MERGED				
		}

		// Gift
		for (int ig = 0; ig < n_gift_swap; ig++) {
			int g_q = gift_qty[ig];
#ifdef VERBOSE
			if (g_q % group_size != 0) error("gift_qty and group size don't match", (double)g_q, (double)group_size);
			cout << g_q << " ";
#endif // VERBOSE
			max_flow.SetNodeSupply(ig + n_child_swap, -(g_q / group_size));
		}
#ifdef VERBOSE
		cout << endl;
		cout << "N of arcs:                         " << max_flow.NumArcs() << endl;
		cout << "N of nodes:                        " << max_flow.NumNodes() << endl;
#endif // VERBOSE


		// SOLVE MAX FLOW MIN COST
		cout << "Start solving Max Flow with Min Cost...";
		int result_status = max_flow.Solve();
		cout << " Done!" << endl;
#ifdef VERBOSE
		cout << "Solver status:                     " << result_status << endl;
		cout << "Maximum flow:                      " << max_flow.MaximumFlow() << endl;
		cout << "Optimal cost:                      " << -max_flow.OptimalCost() / ((double)eff_r * max_child_happiness) << endl;
#endif // VERBOSE


		// OBTAIN RESULT
#ifdef VERBOSE
		cout << "Obtaining result..." << endl;
#endif // VERBOSE
		int gift_total = 0;
		for (int i = 0; i < max_flow.NumArcs(); i++) {
			if (max_flow.Flow(i) > 0) {
				int igp = max_flow.Tail(i);
				int new_g = gift_swap[max_flow.Head(i) - n_child_swap];
				// MERGE WON'T WOKR HERE
				new_gift[igp] = new_g;
				gift_total += group_size;
			}
		}
		cout << "gift_total:                        " << gift_total << endl;

		return;
	}

	// Fast for single
	void RunMaxFlow1(const vector<int>& gift_swap, const vector<int>& gift_qty, 
		vector<vector<int>>& child_swap, vector<int>& new_gift,
		const vector<vector<HappyInt>>& happiness_c,
		const vector<vector<HappyInt>>& happiness_g, int group_size) {

		if (group_size != 1) error("Using max flow 1 on", 1, (double)group_size);

		int n_gift_swap = gift_swap.size();
		int n_child_swap = child_swap.size();  // number of groups

		SimpleMinCostFlow max_flow;

		// ADD ARCS
		cout << "Adding arcs... ";
		for (int igp = 0; igp < child_swap.size(); igp++) {
			int c = child_swap[igp][0];
			for (int ig = 0; ig < gift_swap.size(); ig++) {
				int g = gift_swap[ig];
				max_flow.AddArcWithCapacityAndUnitCost(igp, ig + n_child_swap, 1, 
					-1LL * (eff_r * happiness_c[c][g] + happiness_g[c][g]));
				// AddArcWithCapacityAndUnitCost(tail, head, capicity, cost)
				// tail = index of group in the child_swap list
				// head = index of gift in the gift_swap list + N of child group
			}
#ifdef REBUILD_HAPPY
			happiness_c[c].~vector();  // Save some memory
			happiness_g[c].~vector();  // Save some memory
#endif
		}
#ifdef REBUILD_HAPPY
		happiness_c.~vector();  // Save some memory
		happiness_g.~vector();  // Save some memory
#endif

							  // ASSIGN SUPPLY
		cout << "Assigning supply... " << endl;
		// Child group
		for (int igp = 0; igp < child_swap.size(); igp++) {
			max_flow.SetNodeSupply(igp, 1);  // every child counts as 1
			}

		// Gift
		for (int ig = 0; ig < n_gift_swap; ig++) {
			int g_q = gift_qty[ig];
#ifdef VERBOSE
			if (g_q % group_size != 0) error("gift_qty and group size don't match", (double)g_q, (double)group_size);
			cout << g_q << " ";
#endif // VERBOSE
			max_flow.SetNodeSupply(ig + n_child_swap, -(g_q / group_size));
		}
#ifdef VERBOSE
		cout << endl;
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
#ifdef VERBOSE
		cout << "Obtaining result..." << endl;
#endif // VERBOSE
		int gift_total = 0;
		for (int i = 0; i < max_flow.NumArcs(); i++) {
			if (max_flow.Flow(i) > 0) {
				int igp = max_flow.Tail(i);
				int new_g = gift_swap[max_flow.Head(i) - n_child_swap];
				// MERGE WON'T WOKR HERE
				new_gift[igp] = new_g;
				gift_total += group_size;
			}
		}
		cout << "gift_total:                        " << gift_total << endl;

		return;
	}
}  // namespace operations_research


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
	cout << "Time limit:                        " << time_limit << "s" << endl;
#endif // !NO_INIT

	// ====== MC on Max Flow ======
	int t0mc = time(0);	
#ifdef CONST_N_GIFT_SWAP
	vector<int> gift_swap(n_gift_swap, -1);       // Gift you want to swap their assigned children
	vector<int> gift_swap_count(n_gift_swap, 0);  // N of child been swapped with gift in the gift_swap
#endif // CONST_N_GIFT_SWAP

	for (int step_mf = 1; step_mf <= n_mf; step_mf++) {
		cout << endl;
		cout << "#####  MF run " << step_mf << " / " << n_mf << "  #####"<< endl;
		double time_run = -time(NULL);  // Runtime for the large loop run


		// Copy best answ to this answ
		vector<int> answ_this(answ);
		long long int C_this = C_best;
		long long int G_this = G_best;
		double this_happy = best_happy;

		
		// Update eff_r and rebuild happiness map if they are different
		//Eff_R eff_r_before = eff_r;
		//update_eff_r(C_this, G_this);
#ifdef REBUILD_HAPPY
		vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
		vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
		build_happines_all(happiness_c, happiness_g);
#endif // REBUILD_HAPPY

		// ==============================================================

		// Group size
#ifdef FIX_GROUPSIZE
		int group_size = group_size_fixed;
#else
		int group_size = group_size_list[step_mf % group_size_list_len];

		//if (step_mf % 6 < 3) { group_size = 1; }
		//else if (step_mf % 6 == 3) { group_size = 2; }
		//else if (step_mf % 6 == 4) { group_size = 3; }
		//else { group_size = 6; }  // step_mf % 6 == 5
#endif // FIX_GROUPSIZE
		cout << "Group size: " << group_size << endl;



		// Pick gifts for exchange
		unordered_map<int, int> gift_to_ig;                  // Map from giftId to index in gift_swap list
#ifdef CONST_N_GIFT_SWAP
		pick_gifts(gift_swap_count, gift_swap, gift_to_ig);  // Build gift_swap list and gift_to_ig map
#else
		// Variable n_gift_swap depending on group size
		int n_gift_swap = n_gift_swap_list[group_size - 1];
		vector<int> gift_swap(n_gift_swap, -1);       // Gift you want to swap their assigned children
		vector<int> gift_swap_count(n_gift_swap, 0);  // N of child been swapped with gift in the gift_swap
		pick_gifts_rand(gift_swap, gift_to_ig);
#endif
		cout << "n_gift_swap = " << n_gift_swap << endl;

#ifdef VERBOSE
		cout << "gift_swap: ";
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			cout << gift_swap[ig] << " ";
		}
		cout << endl;
#endif // VERBOSE



		// Pick children for exchange
		vector<vector<int>> child_swap;                                      // Child list for swap
		vector<int> gift_qty(gift_swap.size(), 0);       // Default quantity of each gift
		if (group_size == 1) {
			build_child_swap_single_rand(group_size, gift_swap, gift_qty, child_swap, answ_this, gift_to_ig);
		}
		else {  // group_size > 1
			//build_child_swap_rand(group_size, gift_swap, gift_qty, child_swap);  // Build child_swap list and gift_qty
			//build_child_swap_sort(group_size, gift_swap, gift_qty, child_swap, happiness_c, happiness_g);
			int ig1_shift = (n_gift_swap - 1) * ran_generator() + 1;  // group happiness matching gift shift (should be 1 ~ N_gift_swap-1)
			build_child_swap_sort_enhanced(group_size, gift_swap, gift_qty, child_swap, happiness_c, happiness_g, ig1_shift);
		}
		cout << "N of children: " << child_swap.size() << endl;
		vector<int> new_gift(child_swap.size(), -1);                         // The solver will return new gift assignment to this
#ifdef VERBOSE
		int gift_total_before = accumulate(gift_qty.begin(), gift_qty.end(), 0);
		cout << "total gift number before = " << gift_total_before << endl;
#endif // VERBOSE
		
		//int n_inconsistent = 0;
		//vector<int> incon_child;  // first index of the inconsistent twins or triplet in the child_swap

		// ==============================================================

		// Solve the Max Flow Min Cost for gift in the gift_swap list
		double time_solved = -time(NULL);  // Runtime used by max flow solver
		if (group_size == 1) {
			// SOLVER for group_size == 1
			operations_research::RunMaxFlow1(gift_swap, gift_qty, child_swap, new_gift, happiness_c, happiness_g, group_size);
		}
		else {
			// SOLVER for group_size == 2 and 3
			operations_research::RunMaxFlow3(gift_swap, gift_qty, child_swap, new_gift, happiness_c, happiness_g, group_size);
		}
		time_solved += time(NULL);
		cout << "Max flow runtime:                  " << time_solved << " s" << endl;

		// Write the new_gift to answ_this
		for (int igp = 0; igp < child_swap.size(); igp++) {
			for (int ic = 0; ic < group_size; ic++) {
				int c = child_swap[igp][ic];
				answ_this[c] = new_gift[igp];
			}
		}
		// ============================================================== //
		

		// Check difference and efficiency, and return N of child in swapped gift
		for (int ig = 0; ig < gift_swap_count.size(); ig++)  gift_swap_count[ig] = 0;  // Clear previous swap count
		int n_diff = 0;  // N of swap child
		for (int c = 0; c < Nchild; c++) {
			if (answ_this[c] != answ[c]) {
				gift_swap_count[gift_to_ig[answ[c]]]++;  // Add swap count of gift
				n_diff++;
			}
		}
		time_run += time(NULL);
		printf( "N of swap child:                   %d    [ %.6f/s ]\n", n_diff, n_diff / time_run);

		// Update if there is differnce
		if (n_diff > 0) {
			// Check result
			int n_triplet_incon = triplet_consistent(answ_this);
			int n_twin_incon = twin_consistent(answ_this);
			int n_inconsistent = n_triplet_incon + n_twin_incon;
			cout << "N of inconsistent triplets, twins: " << n_triplet_incon << ", " << n_twin_incon << endl;
			check_gift_count(answ_this);
			
			
			// Compute this total child and gift happiness points
#ifdef REBUILD_HAPPY
			vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
			vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
			build_happines_all(happiness_c, happiness_g);
#endif // REBUILD_HAPPY
			C_this = cal_child_all(answ_this, happiness_c);
			G_this = cal_gift_all(answ_this, happiness_g);
			printf("Child happiness points:            %17lld    [ %+lld ]\n", C_this, C_this - C_best);
			printf("Gift  happiness points:            %17lld    [ %+lld ]\n", G_this, G_this - G_best);
			this_happy = cal_all_from_points(C_this, G_this);
			printf("Score:                             %.15f    [ %+.15f ]\n", this_happy, this_happy - best_happy);
						
			// Check if improved and is consistent
			if (this_happy >= best_happy && n_inconsistent == 0) {
				// Update best answ (even if it is an equivalent one)
				for (int c = 0; c < answ.size(); c++)  answ[c] = answ_this[c];
				answ_to_states(answ_this);  // update best states from answ
				
				// Update happiness map and points
				if (C_this != C_best || G_this != G_best) {
					// Update eff_r and Rebuild happiness map
					update_eff_r(C_this, G_this);
					// Update best scores
					C_best = C_this; G_best = G_this;
				}
				
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
#ifdef CONST_N_GIFT_SWAP
	cout << "N gift swap:                       " << n_gift_swap << endl;
	cout << "Keep ratio:                        " << keep_r << endl;
#else
	cout << "N gift swap max:                   " << n_gift_swap_max << endl;
#endif // CONST_N_GIFT_SWAP
cout << "Time limit:                        " << time_limit << "s" << endl;
	cout << "Total runtime:                     " << time(0) - t0cpu << "s" << endl;
	cout << "===================== End of MC Max Flow ======================" << endl;
#ifdef WINDOWS
	cin.get();
#endif // WINDOWS
	return 0;
}