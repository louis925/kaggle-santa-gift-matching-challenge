#define USE_GLOP
#define _CRT_SECURE_NO_WARNINGS

#include <fstream>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cmath>

#include "ortools/graph/min_cost_flow.h"
#include "santa_par.h"
#include "santa_global.h"

using namespace std;

// Return scaled happiness of child c with gift g
// Twin and triplet happiness are averaged.
inline Eff_R happiness_f(int c, int g, 
	const vector<vector<HappyInt>>& happiness_c, 
	const vector<vector<HappyInt>>& happiness_g) {
	return happiness_c[c][g] * eff_r + happiness_g[c][g];
}

namespace operations_research {
	void RunMaxFlowAll(vector<short>& answ_this, vector<int>& unassigned_gifts) {
		cout << "### Run Max Flow" << endl;

		// BUILD HAPPINESS MAP
		vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
		vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
		build_happines_all(happiness_c, happiness_g);

		// BUILD GRAPH
		SimpleMinCostFlow max_flow;

		// Node list:
		// Children : NodeIndex = 0 to Nchild - 1
		// Gifts    : NodeIndex = Nchild to Nchild + Ntypes - 1 
		// Supply   : NodeIndex = Nchild + Ntypes
		// Note: Demand are added at gift nodes. But no supply at children. Flow from supply_node directly to gift 
		//       represents the unassigned gift. The unassigned gifts will give to unhappy children later.
		const NodeIndex supply_node = Nchild + Ntypes;  // index of supply node

		// ADD ARCS
		// Add from children to gifts
		// Adding arcs with only non-minimum happiness for child
		cout << "Adding arcs... ";
		cout << "triplets ";
		for (int c = 0; c < Ntripl; c+=3){
			for (int g = 0; g < Ntypes; g++) {
				Eff_R h = happiness_f(c, g, happiness_c, happiness_g) + triplet_shift;
				if(h > happy_min) max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 3, -h);  // count 3 times
			}
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, c, 3, 0); // dummy flow			
		}
		cout << "twins ";
		for (int c = Ntripl; c < Ntwins; c += 2) {
			for (int g = 0; g < Ntypes; g++) {
				Eff_R h = happiness_f(c, g, happiness_c, happiness_g) + twin_shift;
				if (h > happy_min) max_flow.AddArcWithCapacityAndUnitCost(c, Nchild + g, 2, -h);  // count 2 times				
			}
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, c, 2, 0); // dummy flow			
		}

		int mem_free_batch = 400000;  // Size of batch memory releasing from happiness_c and happiness_g
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
			max_flow.AddArcWithCapacityAndUnitCost(supply_node, Nchild + g, Nqttys, -happy_min);
		}
		cout << endl;

		// ASSIGN SUPPLY
		cout << "Assigning supply... " << endl;
		max_flow.SetNodeSupply(supply_node, Nchild);      // Supply node		
		for (int g = 0; g < Ntypes; g++) {
			max_flow.SetNodeSupply(Nchild + g, -Nqttys);  // Gifts
		}

#ifdef VERBOSE
		cout << "N of arcs:                         " << max_flow.NumArcs() << endl;
		cout << "N of nodes:                        " << max_flow.NumNodes() << endl;
#endif // VERBOSE
		
		// SOLVE MAX FLOW MIN COST
		cout << "Solving Max Flow with Min Cost...";
		int result_status = max_flow.SolveMaxFlowWithMinCost();
		cout << " Done!" << endl;
#ifdef VERBOSE
		cout << "Solver status:                     " << result_status << endl;
		cout << "Maximum flow:                      " << max_flow.MaximumFlow() << endl;
		cout << "Optimal cost:                      " << -max_flow.OptimalCost() / ((double)eff_r * max_child_happiness) << endl;
#endif // VERBOSE


		// OBTAIN RESULT
		int gift_total = 0;
		int happy_assignment = 0;
		int unhappy_assignment = 0;
		int incorrect_assignment = 0;
		vector<int> n_gift_tri_twin(Ntwins, 0);  // Number of gifts the twin or triplet alread have
												 // This is the triplet and twin filling position shift.
												 // Note the number will only be stored at the head of the twin or triplet
		for (int i = 0; i < max_flow.NumArcs(); i++) {
			if (max_flow.Flow(i) > 0) {
				int f = max_flow.Flow(i);
				int c = max_flow.Tail(i);
				int g = max_flow.Head(i) - Nchild;
				gift_total += f;
				if (c < Nchild) {  
					// Children to gifts assignment
					happy_assignment += f;
					if (c < Ntwins) {    // twins or triplets
						while (f > 0) {
							answ_this[c + n_gift_tri_twin[c]] = g;  // fill the new gift
							n_gift_tri_twin[c]++;  // update the number of gift that twin or triplet will have
							f--;
						}
					}
					else {               // singles
						answ_this[c] = g;
					}
				}
				else if (c == supply_node) {
					if (g >= 0) {  // ignore dummy flow to children
						// Non-assigned gifts
						unassigned_gifts[g] += f;
						unhappy_assignment += f;
					}
				}
				else {
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
	cout << "=======================    Max Flow     =======================" << endl;
	int t0cpu = (int)time(0);

	// READ WISHLISTS
	init();
	vector<int> unassigned_gifts(Ntypes, 0);
	vector<short> answ(Nchild, -1);
	cout << "Finished initialization." << endl;
	cout << endl;
		

	// SOLVE THE MIN COST MAX FLOW PROBLEM
	int t0mf = (int)time(0);
	
	operations_research::RunMaxFlowAll(answ, unassigned_gifts);                   // THE ACTUAL SOLVER
	
	cout << "Max flow runtime:                  " << time(0) - t0mf << "s" << endl;
	cout << endl;
	

	// OUTPUT UNFILLED RESULT
	cout << "Write unfilled result...";
	write_conf(answ, "unfilled");
	cout << " Done!" << endl;
	

	// FILL UNASSIGNED CHILDREN WITH UNASSIGNED GIFTS
	cout << "Filling unassigned children with unassigned gifts... ";
	int gs = 0;    // index of next unassigned gift
	for (int c = 0; c < Nchild; c += 1) {
		int g = answ[c];
		if (g < 0) {
#ifdef VERBOSE
			cout << c << " ";
#endif // VERBOSE
			while (unassigned_gifts[gs] == 0) { gs++; }
			answ[c] = gs;
			unassigned_gifts[gs] -= 1;
		}
	}
	cout << endl;
	cout << endl;


	// CHECK FINAL RESULT
	cout << "N of still unassigned children:    " << check_negative(answ) << endl;
	check_unassigned_gifts(unassigned_gifts);

	// Rebuild happines maps
	vector<vector<HappyInt>> happiness_c(Nchild, vector<HappyInt>(Ntypes, points_childOUT * single_r));
	vector<vector<HappyInt>> happiness_g(Nchild, vector<HappyInt>(Ntypes, points_giftsOUT * gift_r));
	build_happines_all(happiness_c, happiness_g);
	
	// Compute total child and gift happiness points
	Eff_R C = cal_child_all(answ, happiness_c);
	Eff_R G = cal_gift_all(answ, happiness_g);
	double happiness = cal_all_from_points(C, G);

	printf("Child happiness points:            %lld\n", C);
	printf("Gift  happiness points:            %lld\n", G);
	printf("This score:                        %.15f\n", happiness);
	update_eff_r(C, G);
	cout << "N of inconsistent triplets:        " << triplet_consistent(answ) << endl;
	cout << "N of inconsistent twins:           " << twin_consistent(answ) << endl;
	check_gift_count(answ);


	// OUTPUT FINAL RESULT
	cout << "Writing final result to " << output_path << " ...";
	write_conf(answ, happiness, true);
	cout << " Done!" << endl;


	cout << "Total runtime:                     " << time(0) - t0cpu << "s" << endl;
	cout << "======================= End of Max Flow =======================" << endl;
#ifdef WINDOWS
	cin.get();
#endif // WINDOWS
	return 0;
}