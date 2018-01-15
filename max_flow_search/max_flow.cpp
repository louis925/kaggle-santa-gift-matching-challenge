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

// Return scaled happiness of child c with gift g
// Twin and triplet are averaged
inline Eff_R happiness_f(int c, int g, 
	const vector<vector<HappyInt>>& happiness_c, 
	const vector<vector<HappyInt>>& happiness_g) {
	return happiness_c[c][g] * eff_r + happiness_g[c][g];
}

namespace operations_research {
	void RunMaxFlowAll(vector<short>& answ_this, 
		vector<int>& unassigned_gifts_single, vector<int>& unassigned_gifts_twin, vector<int>& unassigned_gifts_triplet,
		vector<int>& n_trips, vector<int>& n_twins) {
		cout << "### Run Max Flow" << endl;
		
		// BUILD HAPPINESS MAP
		vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
		vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
		build_happines_all(happiness_c, happiness_g);

		// BUILD GRAPH
		SimpleMinCostFlow max_flow;

		// < Node list >
		// Children          : NodeIndex = 0                 to Nchild            - 1 (same triplet or twin groups are skipped)
		// Gifts for singles : NodeIndex = Nchild            to Nchild +   Ntypes - 1 (NodeIndex = Nchild + g1)
		// Gifts for twins   : NodeIndex = Nchild +   Ntypes to Nchild + 2*Ntypes - 1 (NodeIndex = Nchild +   Ntypes + g2)
		// Gifts for triplets: NodeIndex = Nchild + 2*Ntypes to Nchild + 3*Ntypes - 1 (NodeIndex = Nchild + 2*Ntypes + g3)
		// Supply            : NodeIndex = Nchild + 3*Ntypes
		//
		// Flow: Flow comes from supply_node, goes through each child, and ends at gift nodes.
		//       Flow can also go directly from supply_node to gift nodes. This flow represents the unassigned gift.
		//       The unassigned gifts will be given to unhappy children later in other part of the code.
		//       Before that, the gift for unhappy children will be denoted as -1.
		// Note: Demands are added at gift nodes. But no supply at child nodes. Instead supply are given at supply_node.
		const NodeIndex supply_node = Nchild + 3*Ntypes;  // index of supply node

		// ADD ARCS
		// Add from children to gifts
		// with only non-minimum happiness for child
		cout << "Adding arcs... ";
		cout << "triplets ";
		for (int c = 0; c < Ntripl; c+=3){
			for (int g = 0; g < Ntypes; g++) {
				Eff_R h = happiness_f(c, g, happiness_c, happiness_g);
				if (h > happy_min) max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + 2 * Ntypes + g, 1, -3 * h);  // count 1 times
			}
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, c, 1, 0); // dummy flow			
		}
		cout << "twins ";
		for (int c = Ntripl; c < Ntwins; c += 2) {
			for (int g = 0; g < Ntypes; g++) {
				Eff_R h = happiness_f(c, g, happiness_c, happiness_g);
				if (h > happy_min) max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + Ntypes + g, 1, -2 * h);  // count 1 times				
			}
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, c, 1, 0); // dummy flow			
		}

		const int mem_free_batch = 400000;  // Size of batch memory releasing from happiness_c and happiness_g
											// release memory for building graph
											// range: Ntwins ~ Nchild

		cout << "singles ";
		for (int c = Ntwins; c < Nchild; c += 1) {
			if (c % 100000 == 0) cout << c << " ";
			for (int g = 0; g < Ntypes; g++) {
				Eff_R h = happiness_f(c, g, happiness_c, happiness_g);
				if (h > happy_min) max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 1, -h);
			}
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, c, 1, 0); // dummy flow
			
			// Release memory batch by batch
			if (c % mem_free_batch == 0) {
				cout << "re ";
				for (int cr = c - mem_free_batch; cr < c; cr++) {
					happiness_c[cr].~vector();
					happiness_g[cr].~vector();
				}
			}
		}
		happiness_c.~vector();  // Save some memory
		happiness_g.~vector();  // Save some memory
		
		// Adding arcs from supply_node to gifts (unassigned gifts)
		cout << "non-assigned gifts ";
		for (int g = 0; g < Ntypes; g++) {
			// Since each unassigned gift means there is one kid who got gift not from his or gift's lists,
			// the happiness of this edge will be happy_min.

			// Gifts for triplets
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, Nchild + 2 * Ntypes + g, n_trips[g], -3 * happy_min);

			// Gifts for twins
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, Nchild + Ntypes + g, n_twins[g], -2 * happy_min);

			// Gifts for singles
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, Nchild + g, Nqttys - 2 * n_twins[g] - 3 * n_trips[g], -happy_min);
		}
		cout << endl;


		// ASSIGN SUPPLY
		cout << "Assigning supply... " << endl;
		max_flow.SetNodeSupply(supply_node, Nchild - 2*Ntripl - (Ntwins - Ntripl));           // Supply node		
		// note each twin or triplet flow only counts as one
		for (int g = 0; g < Ntypes; g++) {
			max_flow.SetNodeSupply(Nchild + 2 * Ntypes + g, -n_trips[g]);                     // Gifts for triplets
			max_flow.SetNodeSupply(Nchild + Ntypes + g, -n_twins[g]);                         // Gifts for twins
			max_flow.SetNodeSupply(Nchild + g, -(Nqttys - 2 * n_twins[g] - 3 * n_trips[g]));  // Gifts for singles
		}

#ifdef VERBOSE
		const int n_child_nodes = Nchild - Ntwins + (Ntwins - Ntripl) / 2 + Ntripl / 3;
		cout << "N of arcs:                         " << max_flow.NumArcs()
			<< " (expect ~" << n_child_nodes * (Nchprf + 1) + 3 * Ntypes << ")" << endl;
		cout << "N of nodes:                        " << max_flow.NumNodes() 
			<< " (expect " << n_child_nodes + 3 * Ntypes + 1 << " or " << Nchild + 3 * Ntypes + 1 << ")" << endl;
#endif // VERBOSE
		

		// SOLVE MAX FLOW MIN COST
		cout << "Solving Max Flow with Min Cost... ";
		int result_status = max_flow.SolveMaxFlowWithMinCost();
		cout << "Done!" << endl;
#ifdef VERBOSE
		cout << "Solver status:                     " << result_status << endl;
		cout << "Maximum flow:                      " << max_flow.MaximumFlow() 
			<< " (expect " << n_child_nodes << ")" << endl;
		cout << "Optimal cost:                      " << -max_flow.OptimalCost() / ((double)6 * eff_r * max_child_happiness) << endl;
#endif // VERBOSE


		// OBTAIN RESULT
		int gift_total = 0;
		int happy_assignment = 0;
		int unhappy_assignment = 0;
		int incorrect_assignment = 0;
		for (int i = 0; i < max_flow.NumArcs(); i++) {
			if (max_flow.Flow(i) > 0) {
				int f = max_flow.Flow(i);
				int c = max_flow.Tail(i);
				gift_total += f;
				if (c < Ntripl) {           // Triplets to gifts assignment
					int g = max_flow.Head(i) - Nchild - 2 * Ntypes;
					answ_this[c    ] = g;
					answ_this[c + 1] = g;
					answ_this[c + 2] = g;
					happy_assignment += 3;
				}
				else if (c < Ntwins) {      // Twins to gifts assignment
					int g = max_flow.Head(i) - Nchild - 1 * Ntypes;
					answ_this[c    ] = g;
					answ_this[c + 1] = g;
					happy_assignment += 2;
				}
				else if (c < Nchild) {      // Singles to gifts assignment
					int g = max_flow.Head(i) - Nchild;
					answ_this[c] = g;
					happy_assignment += 1;					
				}
				else if (c == supply_node) {
					int g = max_flow.Head(i) - Nchild;
					if (g >= 2*Ntypes) {    // Non-assigned gifts for triplets
						g -= 2 * Ntypes;
						unassigned_gifts_triplet[g]++;
						unhappy_assignment += 3;
					}
					else if (g >= Ntypes) { // Non-assigned gifts for twins
						g -= Ntypes;
						unassigned_gifts_twin[g]++;
						unhappy_assignment += 2;
					}
					else if (g >= 0) {      // Non-assigned gifts for singles
						unassigned_gifts_single[g]++;
						unhappy_assignment += 1;
					}
					// ignore negative g (supply_node to children)
				}
				else {
					cout << "Incorrect assignment?" << endl;
					incorrect_assignment += f;
				}
			}
		}

		cout << "gift_total:                        " << gift_total << endl;
		cout << "Happy assignment:                  " << happy_assignment << endl;
		cout << "Unhappy assignment:                " << unhappy_assignment << endl;
		cout << "Incorrect assignment:              " << incorrect_assignment << endl;

		return;
	}
}  // namespace operations_research


int main() {
	cout << "======================    MC Max Flow    ======================" << endl;
	int t0cpu = time(0);
#ifndef NO_INIT
	init();                                      // Read wishlists and build happiness maps
	
	long long int C_init = cal_child_all_ori(answ);  // Initial total child happiness points
	long long int G_init = cal_gift_all_ori(answ);   // Initial total gift happiness points
	double init_happy = cal_all_from_points(C_init, G_init);
	
	vector<short> answ_init(answ);               // Initial answer
												 // Note answ is the best answ, initialized by init() from input data
	long long int C_best = C_init;               // Total child happiness points
	long long int G_best = G_init;               // Total gift happiness points
	double best_happy = init_happy;              // Previous best happiness

	printf( "Initial score:                     %.15f\n", init_happy);
	cout << "Child / Gifts happiness points:    " << C_init << " / " << G_init << endl;
	Eff_R eff_r_init = update_eff_r(C_init, G_init);  // Initial effective happiness ratio
	
	// Compute minimum of approximate happiness points
	happy_min = points_childOUT * single_r * eff_r + points_giftsOUT * gift_r;
	
#ifdef FIX_EFF_R
	cout << "Using fixed eff_r:                 " << eff_r << endl;
#endif // FIX_EFF_R
#ifdef ZERO_GIFT_R
	cout << "Initial gift_r = 0" << endl;
#endif // ZERO_GIFT_R

	
	cout << "Finished initialization." << endl;	
#endif // !NO_INIT

	// ====== MC on Max Flow ======
	int t0mc = time(0);

	// Copy best answ to this answ
	long long int C_this = C_best;
	long long int G_this = G_best;
	double this_happy = best_happy;


	// COMPUTE INIT TWIN AND TRIPLET SIZES
	vector<int> n_trip_answ(Ntypes, 0);
	vector<int> n_twin_answ(Ntypes, 0);
	for (int c = 0; c < Ntripl; c++) { n_trip_answ[answ[c]] ++; }
	for (int c = Ntripl; c < Ntwins; c++) { n_twin_answ[answ[c]] ++; }

	// Round the number of twins and triplets in each gifts
	vector<int> n_trips(Ntypes, 0);
	vector<int> n_twins(Ntypes, 0);
	int n_trips_total = 0, n_twins_total = 0;
	for (int g = 0; g < Ntypes - 1; g++) { 
		n_trips[g] = lround(n_trip_answ[g] / 3.0);
		n_twins[g] = lround(n_twin_answ[g] / 2.0);
		n_trips_total += n_trips[g];
		n_twins_total += n_twins[g];
	}
	// last gift type
	n_trips[Ntypes - 1] = Ntripl/3 - n_trips_total;
	n_twins[Ntypes - 1] = (Ntwins - Ntripl) / 2 - n_twins_total;
	
	cout << "n_trips: "; print_vector(n_trips);
	cout << "n_twins: "; print_vector(n_twins);

	cin.get();

	// SOLVE THE MIN COST MAX FLOW PROBLEM
	cout << endl;
	double time_solved = -time(NULL);
	vector<int> unassigned_gifts_single(Ntypes, 0);
	vector<int> unassigned_gifts_twin(Ntypes, 0);
	vector<int> unassigned_gifts_triplet(Ntypes, 0);
	vector<short> answ_this(Nchild, -1);
	
	// the actual solver
	operations_research::RunMaxFlowAll(answ_this, 
		unassigned_gifts_single, unassigned_gifts_twin, unassigned_gifts_triplet,
		n_trips, n_twins);
	
	time_solved += time(NULL);
	cout << "Max flow runtime:                  " << time_solved << " s" << endl;
	

	// OUTPUT UNFILLED RESULT
	cout << "Write unfilled result...";
	write_conf(answ_this, "unfilled");
	cout << " Done!" << endl;
	
	// FILL UNASSIGNED CHILDREN WITH UNASSIGNED GIFTS
	cout << "Filling unassigned children with unassigned gifts... ";
	cout << "triplets ";
	int gu = 0;    // index of next unassigned gift
	for (int c = 0; c < Ntripl; c += 3) {
		int g = answ_this[c];
		if (g < 0) {
			cout << c << " ";
			while (unassigned_gifts_triplet[gu] == 0) { gu++; }
			answ_this[c    ] = gu;
			answ_this[c + 1] = gu;
			answ_this[c + 2] = gu;
			unassigned_gifts_triplet[gu] -= 1;
		}
	}
	cout << "twins ";
	gu = 0;
	for (int c = Ntripl; c < Ntwins; c += 2) {
		int g = answ_this[c];
		if (g < 0) {
			cout << c << " ";
			while (unassigned_gifts_twin[gu] == 0) { gu++; }
			answ_this[c] = gu;
			answ_this[c + 1] = gu;
			unassigned_gifts_twin[gu] -= 1;
		}
	}
	cout << "singles ";
	gu = 0;
	for (int c = Ntwins; c < Nchild; c += 1) {
		int g = answ_this[c];
		if (g < 0) {
			cout << c << " ";
			while (unassigned_gifts_single[gu] == 0) { gu++; }
			answ_this[c] = gu;
			unassigned_gifts_single[gu] -= 1;
		}
	}	
	cout << "done" << endl;
	cout << endl;


	// CHECK RESULT
	cout << "Number of still unassigned children: " << check_negative(answ_this) << endl;
	cout << "Triplet ";
	check_unassigned_gifts(unassigned_gifts_triplet);
	cout << "Twin ";
	check_unassigned_gifts(unassigned_gifts_twin);
	cout << "Single ";
	check_unassigned_gifts(unassigned_gifts_single);

	vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
	vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
	build_happines_all(happiness_c, happiness_g);
	
	// Compute this total child and gift happiness points
	C_this = cal_child_all(answ_this, happiness_c);
	G_this = cal_gift_all(answ_this, happiness_g);
	printf("Child happiness points:            %17lld    [ %+lld ]\n", C_this, C_this - C_best);
	printf("Gift  happiness points:            %17lld    [ %+lld ]\n", G_this, G_this - G_best);
	this_happy = cal_all_from_points(C_this, G_this);
	printf("This score:                        %.15f    [ %+.15f ]\n", this_happy, this_happy - best_happy);
	
	// Output result
	write_conf(answ_this, this_happy, false);
	
	// Check difference and efficiency, and return N of child in swapped gift
	int n_diff = compare_answ_diff(answ, answ_this);  // N of swap child
	printf( "N of swap child:                   %d    [ %.6f/s ]\n", n_diff, n_diff / time_solved);

	// Update if there is differnce
	if (n_diff > 0) {
		// Check result
		cout << "N of inconsistent triplets, twins: " << triplet_consistent(answ_this) << ", "
			<< twin_consistent(answ_this) << endl;
		check_gift_count(answ_this);
		
		printf("Score:                             %.15f    [ %+.15f ]\n", this_happy, this_happy - best_happy);
		cout << "Child / Gifts happiness points:    " << C_this << " / " << G_this << endl;

		// Check if improved and is consistent
		if (this_happy >= best_happy) {
			// Update best answ and scores
			C_best = C_this; G_best = G_this;
			for (int c = 0; c < answ.size(); c++)  answ[c] = answ_this[c];
			
			// Output optimization result
			if (this_happy > best_happy) {
				best_happy = this_happy;
				cout << "Output improved result..." << endl;
				write_conf(answ_this, this_happy, false);
			}
		}
		else {
			cout << "Worse, not using this result..." << endl;
		}
	}  // n_diff > 0
	int tfmc = time(0);
	

	// CHECK FINAL RESULT
	update_eff_r(C_best, G_best);  // Update eff_r to the best one

	cout << endl;
	cout << "### Summary" << endl;
	printf("Final score:                       %.15f    [ %+.15f ]\n", best_happy, best_happy - init_happy);
	printf("Child happiness points:            %17lld    [ %+lld ]\n", C_best, C_best - C_init);
	printf("Gift  happiness points:            %17lld    [ %+lld ]\n", G_best, G_best - G_init);
	printf("Recommend effective ratio:         %17lld    [ %+lld ]\n", eff_r, eff_r - eff_r_init);
	n_diff = compare_answ_diff(answ, answ_init);  // N of swap child
	printf("N of swap child:                   %17d    [ %.6f/s ]\n", n_diff, ((double) n_diff) / (tfmc - t0mc + 1e-3));


	// OUTPUT
	cout << "Writing result...";
	write_conf(answ, best_happy, false);
	write_conf(answ, best_happy, true);
	cout << " Done!" << endl;


	// Remind you the parameters
	cout << "### Remind you the parameters" << endl;
	cout << "Total runtime:                     " << time(0) - t0cpu << "s" << endl;
	cout << "===================== End of MC Max Flow ======================" << endl;
#ifdef WINDOWS
	cin.get();
#endif // WINDOWS
	return 0;
}