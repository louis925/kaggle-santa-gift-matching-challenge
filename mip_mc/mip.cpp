#define USE_CBC
#define _CRT_SECURE_NO_WARNINGS

#include <fstream>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm> 

#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"
//#include "ortools/graph/min_cost_flow.h"

#include "santa_par.h"
#include "santa_global.h"

using namespace std;

inline double ran_generator() {
	return rand() / ((double)RAND_MAX + 1.0);
}

namespace operations_research {
	int RunMIP(vector<int> &gift_swap, vector<int> &gift_swap_count) {
		int gift_n = gift_swap.size();
		
		vector<int> answ_before(answ);  // Original answer

		MPSolver solver("MIP", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);

		// 100*100*1000 : 3GB (300 B per variable)
		
		// Variables
		cout << "Creat variables... " << gift_n*Nqttys << "x" << gift_n << endl;
		vector<vector<MPVariable*>> variables(gift_n*Nqttys, vector<MPVariable*>(gift_n));  // pointer to the variables in the solver [Nqttys * ig + ic][jg]
		for (int ig = 0; ig < gift_n; ig++) {  // Loop in the original gift of swap chlid
			int g = gift_swap[ig];  // Original gift of swap child
			int ic = 0;  // Loop in the swap child list with gift g
			while (ic < Nqttys) {
				int c = states[g][ic];
				for (int jg = 0; jg < gift_n; jg++) {  // Loop in gift assignment switch
					variables[Nqttys * ig + ic][jg] = solver.MakeBoolVar("");
				}

				// Skip following twin and triplets
				if (c >= Ntwins) ic += 1;       // single
				else if (c >= Ntripl) ic += 2;  // twin
				else ic += 3;                   // triplet
			}
		}
		
		// Objective: maximize happiness
		cout << "Set objective";
		MPObjective* const objective = solver.MutableObjective();
		for (int ig = 0; ig < gift_n; ig++) {  // Loop in the original gift of swap chlid
			int g = gift_swap[ig];  // Original gift of swap child
			int ic = 0;  // Loop in the swap child list with gift g
			while (ic < Nqttys) {
				int c = states[g][ic];
				for (int jg = 0; jg < gift_n; jg++) {  // Loop in gift assignment switch
					int g_swap = gift_swap[jg];  // Target new gift
					objective->SetCoefficient(variables[Nqttys * ig + ic][jg], happiness[c][g_swap]);
				}

				// Skip following twin and triplets
				if (c >= Ntwins) ic += 1;       // single
				else if (c >= Ntripl) ic += 2;  // twin
				else ic += 3;                   // triplet
			}
		}
		objective->SetMaximization();

		// Constraint: for each child, gift_n_total = 1
		cout << " child constraint";
		for (int ig = 0; ig < gift_n; ig++) {          // Loop in the original gift of swap chlid
			int g = gift_swap[ig];  // Original gift of swap child
			int ic = 0;                                // Loop in the swap child list with gift g
			while (ic < Nqttys) {
				int c = states[g][ic];
				MPConstraint* const constraint = solver.MakeRowConstraint(1, 1);  // for each child, gift_n_total = 1
				for (int jg = 0; jg < gift_n; jg++) {  // Loop in gift assignment switch
					constraint->SetCoefficient(variables[Nqttys * ig + ic][jg], 1);
				}

				// Skip following twin and triplets
				if (c >= Ntwins) ic += 1;       // single
				else if (c >= Ntripl) ic += 2;  // twin
				else ic += 3;                   // triplet
			}
		}

		// Constraint: for each gift, gift_n_total = 1000
		//             twin and triplets should have weight factor.
		cout << " gift total constraint";
		for (int jg = 0; jg < gift_n; jg++) {          // Loop in gift assignment switch
			MPConstraint* const constraint = solver.MakeRowConstraint(Nqttys, Nqttys);  // for each gift, gift_quantity = 1000
			for (int ig = 0; ig < gift_n; ig++) {      // Loop in the original gift of swap chlid
				int g = gift_swap[ig];  // Original gift of swap child
				int ic = 0;                            // Loop in the swap child list with gift g
				while (ic < Nqttys) {
					int c = states[g][ic];

					// Skip following twins and triplets
					if (c >= Ntwins) {       // single
						constraint->SetCoefficient(variables[Nqttys * ig + ic][jg], 1);
						ic += 1;
					}
					else if (c >= Ntripl) {  // twin
						constraint->SetCoefficient(variables[Nqttys * ig + ic][jg], 2);  // count 2 times
						ic += 2;
					}
					else {                   // triplet
						constraint->SetCoefficient(variables[Nqttys * ig + ic][jg], 3);  // count 3 times
						ic += 3;
					}
				}
			}
		}
		cout << endl;

		// Solve
		cout << "Solving MIP...";
		int result_status = solver.Solve();
		double time_solved = solver.wall_time() / 1000.0;
		printf(" Done!    %.2f s\n", time_solved);
		cout << "MIP result status: " << result_status << endl;

		// Obtain solution
		//cout << "Obtain solution to answ" << endl;
		for (int ig = 0; ig < gift_n; ig++) {          // Loop in the original gift of swap chlid
			int g = gift_swap[ig];  // Original gift of swap child
			int ic = 0;                                // Loop in the swap child list with gift g
			while (ic < Nqttys) {
				int c = states[g][ic];

				int new_g = -1;
				int n_sol = 0;  // num of non zero solution
				for (int jg = 0; jg < gift_n; jg++) {  // Loop in gift assignment switch
					if (variables[Nqttys * ig + ic][jg]->solution_value() == 1) {
						new_g = gift_swap[jg];
						n_sol++;
					}
				}

				// debug 
				if (new_g == -1) {
					printf("didn't find g for %d", c);
				}

				if (n_sol > 1) { printf("c = %d has %d solutions\n", c, n_sol); }  // no problem

				// Skip following twin and triplets
				if (c >= Ntwins) {       // single
					answ[c] = new_g;
					ic += 1;       
				}
				else if (c >= Ntripl) {  // twin
					answ[c] = new_g;
					answ[c + 1] = new_g;
					ic += 2;
				}
				else {                   // triplet
					answ[c] = new_g;
					answ[c + 1] = new_g;
					answ[c + 2] = new_g;
					ic += 3;                   
				}
			}
		}
		
		// The objective value of the solution.
		printf("Optimal objective value = %.9e\n", objective->Value());

		// Check difference and efficiency, and return N of child in swapped gift
		unordered_map<int, int> gift_to_ig;                                            // Map from giftId to index in gift_swap list
		for (int ig = 0; ig < gift_swap.size(); ig++)  gift_to_ig[gift_swap[ig]] = ig;
		for (int ig = 0; ig < gift_swap_count.size(); ig++)  gift_swap_count[ig] = 0;  // Clear previous record
		int n_diff = 0;  // N of swap child
		for (int c = 0; c < Nchild; c++) {
			if (answ[c] != answ_before[c]) {
				gift_swap_count[gift_to_ig[answ[c]]]++;  // Add swap count of gift
				n_diff++;
			}
		}
		printf("N of swap child: %d    [%.6f /s]\n", n_diff, n_diff / time_solved);

		return n_diff;
	}
}  // namespace operations_research


int main() {
	cout << "========================    MC  MIP    ========================" << endl;
	int t0cpu = time(0);
#ifndef NO_INIT
	double init_happy = init();  // Read wishlists and build happiness maps
	double best_happy = init_happy;  // Previous best happiness
	vector<int> answ_init(answ);  // Initial answer
	long long int total_child_happiness_points = cal_child_all(answ);
	long long int total_gift_happiness_points = cal_gift_all(answ);
	long long int total_child_happiness_points_init = total_child_happiness_points;
	long long int total_gift_happiness_points_init = total_gift_happiness_points;
	
	printf("Initial score: %.10f\n", init_happy);
	cout << "Child / Gifts happiness points = " << total_child_happiness_points
		<< " / " << total_gift_happiness_points << endl;
	cout << "Initial effective happiness ratio = " << eff_r_init;
#ifdef ZERO_GIFT_R
	cout << " with gift_r = 0";
#endif // ZERO_GIFT_R
	cout << endl;

	//wlist_child.~vector();
#ifndef NO_GIFT_MAP
	//wlist_gifts.~vector();
#endif // !NO_GIFT_MAP
	cout << "Finishing initialization." << endl;

#endif // !NO_INIT


	// ====== MC on MIP ======
	int t0mc = time(0);
	vector<int> gift_swap(n_gift_swap, -1);       // Gift you want to swap their assigned children
	vector<int> gift_swap_count(n_gift_swap, 0);  // N of child been swapped with gift in the gift_swap
	for (int step_mip = 1; step_mip <= n_mip; step_mip++) {
		cout << endl;
		cout << "### MIP run " << step_mip << " / " << n_mip << endl;

		// Gifts for exchange
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
				g = Ntypes * ran_generator();
				it = find(gift_swap.begin(), gift_swap.end(), g);  // Not to take the same gift including prevous searched gift
			} while (it != gift_swap.end());
			gift_swap[i] = g;
		}
		cout << "gift_swap: ";
		for (int ig = 0; ig < gift_swap.size(); ig++) {
			cout << gift_swap[ig] << " ";
		}
		cout << endl;

		// Solve the Mixed Integer Problem for gift in the gift_swap list
		int n_diff = operations_research::RunMIP(gift_swap, gift_swap_count);
		answ_to_states();  // Update states

		// Check result
		cout << "N of inconsistent triplets, twins: " << triplet_consistent(answ) << ", "
			<< twin_consistent(answ) << endl;
		check_gift_count(answ);
		double this_happy = cal_all(answ);
		total_child_happiness_points = cal_child_all(answ);
		total_gift_happiness_points = cal_gift_all(answ);
		printf("Score: %.10f    improved by %+.10f\n", this_happy, this_happy - best_happy);
		cout << "Child / Gifts happiness points = " << total_child_happiness_points
			<< " / " << total_gift_happiness_points << endl;
		update_eff_r(total_child_happiness_points, total_gift_happiness_points);
		best_happy = this_happy;

		if (step_mip % step_output == 0) {
			if (n_diff > 0) {
				cout << "Output optimization result" << endl;
				write_conf(best_happy, false);
#ifndef NO_GIFT_R
				build_happines_all();  // Update happiness map with new eff_r
#endif // !NO_GIFT_R
			}
		}
	}
	int tfmc = time(0);
	

	// Check Final Result
	cout << endl;
	cout << "### Summary" << endl;
	printf("Final score: %.10f    [%+.10f]\n", best_happy, best_happy - init_happy);
	printf("Child happiness points: %d    [%+d]\n", total_child_happiness_points, 
		total_child_happiness_points - total_child_happiness_points_init);
	printf("Gift  happiness points: %d    [%+d]\n", total_gift_happiness_points,
		total_gift_happiness_points - total_gift_happiness_points_init);
	printf("Recommend effective ratio: %d    [%+d]\n",
		eff_r, eff_r - eff_r_init);
	int n_diff = 0;  // N of swap child
	for (int c = 0; c < Nchild; c++) {
		if (answ[c] != answ_init[c]) n_diff++;
	}
	printf("N of swap child: %d    (%.2e /s)\n", n_diff, ((double) n_diff) / tfmc);


	// Output
	cout << "Writing result...";
	if (n_mip % step_output != 0) {
		write_conf(best_happy, false);
	}
	write_conf(best_happy, true);
	cout << " Done!" << endl;

	cout << "Total runtime: " << time(0) - t0cpu << "s" << endl;
	cout << "======================== End of MC MIP ========================" << endl;

	cin.get();
	return 0;
}