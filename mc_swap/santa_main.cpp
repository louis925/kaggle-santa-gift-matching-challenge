/* Original by Sky Huang */
/* Modified by Louis Yang 2017.12.18 PM 09:02 */
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include "santa_par.h"
#include "santa_global.h"

using namespace std;

inline double ran_generator() {
	return rand() / ((double)RAND_MAX + 1.0);
}

inline int cal_single_chlid(int ichild, int igift) {
	return child_happy_map[ichild][igift];
}


inline int cal_single_gifts(int ichild, int igift) {
	return gift_happy_map[ichild][igift];
}

inline double cal_diff(long long int C, long long int G, int dC, int dG) {
#ifdef APP_DIFF
	return 3.0 * ((double)G * G * dG + (double)child_gifts_happiness_ratio_3 * C * C * dC
		) / max_gifts_happiness_3;  // approximate happiness difference
#else
	return ((double) dG * dG * dG + (double)child_gifts_happiness_ratio_3 * dC * dC * dC
		    + 3.0 * G * dG * dG + 3.0 * child_gifts_happiness_ratio_3 * C * dC * dC
		    + 3.0 * G * G * dG + 3.0 * child_gifts_happiness_ratio_3 * C * C * dC
		   ) / max_gifts_happiness_3;  // exact happiness difference
#endif
}

inline bool is_single(int ichild) { return ichild >= Ntwins; }
inline bool is_not_single(int ichild) { return ichild < Ntwins; }
inline bool is_twin(int ichild) { return ichild < Ntwins && ichild >= Ntripl; }
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

//************************************************************************************************************
int main() 
{
    cout << "********** Metropolis MC simulation on Santa problem **********\n" << endl;
	cout << "########## Initializing System... ##########" << endl;
	long long int timestep = 0;     // number of improvement step that has been made
	long long int unhappy_found = 0;// number of improvement made from unhappy child search
	double happiness = init();      // current happiness during the simulation
	double best_happy = happiness;
	double init_happy = happiness;  // initial happiness before the simulation
	long long int total_child_happiness_points = cal_child_all(); 
	long long int total_gift_happiness_points = cal_gift_all();

	vector<int> answ_init(Nchild, -1);  // initial answer
	for (int c = 0; c < Nchild; c++) {
		answ_init[c] = gift_list[c][0];
	}

    cout << "Initialization completed. The initial happiness= " << happiness << endl;	
	cout << "Effective happiness ratio = " << eff_r << endl;
#ifdef WEAK_SEARCH_MAIN
	cout << "non_weak_search = " << non_weak_search << endl;
#else
	cout << "weak_search = " << weak_search << endl;
#endif
	cout << "unhappy_threshold = " << unhappy_threshold << endl;

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	printf("Steps: %lld, happiness: %.10f\n", timestep, happiness);
	
	// TEMP variables in the calculation ====================================
	bool is_improved, found_all_child; 
	int g1, g2; // gift 1 and 2 id
	int c11, c12, c13; // children who were assigned to g1
	int c21, c22, c23; // children who were assigned to g2
	int gc11;  // index of child 1 in gift list 1
	int gc12;  // index of child 1 pair in gift list 1
	int gc13;  // index of child 1 pair 2 in gift list 1
	int gc21;  // index of child 2 in gift list 2
	int gc22;  // index of child 2 pair in gift list 2
	int gc23;  // index of child 2 pair 2 in gift list 2

	int child_happiness_points_diff;  // difference in child happiness points
	int gift_happiness_points_diff;   // difference in gift happiness points

	int hc11_before, hc12_before, hc13_before, hc21_before, hc22_before, hc23_before;  // child happiness points before swap
	int hg11_before, hg12_before, hg13_before, hg21_before, hg22_before, hg23_before;  // gift happiness points before swap
	int hc11_after, hc12_after, hc13_after, hc21_after, hc22_after, hc23_after;        // child happiness points after swap
	int hg11_after, hg12_after, hg13_after, hg21_after, hg22_after, hg23_after;        // gift happiness points after swap

	long long int h11_before, h21_before, h12_before, h22_before, h13_before, h23_before;     // approx. happiness points before swap
	long long int h11_after, h21_after, h12_after, h22_after, h13_after, h23_after;           // approx. happiness points after swap

	double ediff; double hdiff;  // difference in happiness

	int mini_steps = 1; // counter for the main while loop
	int i;

	int t0cpu = time(0);

	// ============================= Main search loop =============================
    while(timestep < step_end && time(0) - t0cpu < max_runtime){
        is_improved = false;
		found_all_child = false;

		// PICK GIFT 1, CHILD 1
#ifdef WEAK_SEARCH_MAIN
		if (mini_steps % non_weak_search != 0 && unhappy_child.size() > 0)  // More weak search
#else
		if (mini_steps % weak_search == 0 && unhappy_child.size() > 0)  // More non-weak search
#endif
		{
			// Force to pick from unhappy_child (weak search)
			i = (int)(ran_generator()*unhappy_child.size());
			c11 = unhappy_child[i];  // Pick an unhappy child
			g1 = gift_list[c11][0];
			gc11 = gift_list[c11][1];
		}
		else
		{
			g1 = (int)(ran_generator()*Ntypes); // gift 1
			gc11 = (int)(ran_generator()*Nqttys); // PICK CHILD IN GIFT 1
			c11 = states[g1][gc11];
		}
        
		// PICK GIFT 2
		if (mini_steps % pref_gift_search == 0)
		{
			do {
				g2 = wlist_child[c11][(int)(ran_generator()*Nchprf)]; // Pick gift 2 from c11's wishlist
			} while (g2 == g1);
		}
		else
		{
			g2 = (int)(ran_generator()*(Ntypes - 1));
			if (g2 >= g1) g2 += 1;  // Not to take g1
		}

		// ==================================== 1v1 SWAP ====================================
		if (is_single(c11)) {
			// GET CHILDREN FOR GIFT 2
			for (int j=0; j < match_times; j++)	{  // get a single for c21
				gc21 = (int)(ran_generator()*Nqttys);
				c21 = states[g2][gc21];
				if (is_single(c21))	{  // get only single
					found_all_child = true;
					break;
				}
			}

			// COMPUTE DIFFERENCE AND SWAP
			if (found_all_child)
			{
				// compute child and gift happiness points
				hc11_before = child_happy_map[c11][g1]; hg11_before = gift_happy_map[c11][g1];
				hc21_before = child_happy_map[c21][g2]; hg21_before = gift_happy_map[c21][g2];
				hc11_after = child_happy_map[c11][g2]; hg11_after = gift_happy_map[c11][g2];
				hc21_after = child_happy_map[c21][g1]; hg21_after = gift_happy_map[c21][g1];
			
				child_happiness_points_diff = hc11_after  + hc21_after
											-(hc11_before + hc21_before);
				gift_happiness_points_diff  = hg11_after  + hg21_after
											-(hg11_before + hg21_before);

				ediff = cal_diff(total_child_happiness_points, total_gift_happiness_points,
					             child_happiness_points_diff, gift_happiness_points_diff);

				//#ifdef ZERO_TEMP
				if (ediff >= 0) { // perform swap
								  //#else
								  //          if(ediff >= 0 || ran_generator() < exp(ediff/temperature)){ // perform swap
								  //#endif
					states[g2][gc21] = c11;
					states[g1][gc11] = c21;

					gift_list[c11][0] = g2; gift_list[c11][1] = gc21;
					gift_list[c21][0] = g1; gift_list[c21][1] = gc11;

					total_child_happiness_points += child_happiness_points_diff;
					total_gift_happiness_points += gift_happiness_points_diff;

					if (ediff > 0) {
						happiness += ediff;
						is_improved = true;
					}

					// Update unhappy_child list -------------------
					h11_before = eff_r * hc11_before + hg11_before;  // approx. happiness
					h21_before = eff_r * hc21_before + hg21_before;
					h11_after = eff_r * hc11_after + hg11_after;
					h21_after = eff_r * hc21_after + hg21_after;

					if (h11_before < unhappy_threshold && h11_after >= unhappy_threshold) {
						int u11 = 0;
						while (unhappy_child[u11] != c11) u11++;
						unhappy_child.erase(unhappy_child.begin() + u11);
					}
					if (h21_before < unhappy_threshold && h21_after >= unhappy_threshold) {
						int u21 = 0;
						while (unhappy_child[u21] != c21) u21++;
						unhappy_child.erase(unhappy_child.begin() + u21);
					}
					if (h11_before >= unhappy_threshold && h11_after < unhappy_threshold) {
						unhappy_child.push_back(c11);
					}
					if (h21_before >= unhappy_threshold && h21_after < unhappy_threshold) {
						unhappy_child.push_back(c21);
					}
				}
			}
		}
		// ==================================== 2v2 SWAP ====================================
		else if (is_not_triplet(c11)) {  // c11 == twin; 2v11 or 2v2
			get_twins(c11, c12);
			gc12 = gift_list[c12][1];
						
			// GET CHILDREN FOR GIFT 2
			bool found_c21 = false;

			// GET C21
			for (int j = 0; j < match_times; j++) {  // get a twin or single for c21
				gc21 = (int)(ran_generator()*Nqttys);
				c21 = states[g2][gc21];
				if (is_not_triplet(c21)) {  // get only twin or single
					found_c21 = true;
					break;
				}
			}

			// GET C22
			if (found_c21) {
				if (is_single(c21)) { // 2v11
					for (int j = 0; j < match_times; j++) {  // get a single for c22
						gc22 = (int)(ran_generator()*(Nqttys - 1));
						if (gc22 >= gc21) gc22 += 1;  // not to take gc21
						c22 = states[g2][gc22];
						if (is_single(c22)) {  // get only single
							found_all_child = true;
							break;
						}
					}
				}
				else { // c21 == twin; 2v2
					get_twins(c21, c22);  // get twin of c21 as c22
					gc22 = gift_list[c22][1];
					found_all_child = true;
				}
			}

			// COMPUTE DIFFERENCE AND SWAP
			if (found_all_child)
			{
				hc11_before = child_happy_map[c11][g1]; hg11_before = gift_happy_map[c11][g1];
				hc12_before = child_happy_map[c12][g1]; hg12_before = gift_happy_map[c12][g1];
				hc21_before = child_happy_map[c21][g2]; hg21_before = gift_happy_map[c21][g2];
				hc22_before = child_happy_map[c22][g2]; hg22_before = gift_happy_map[c22][g2];
				hc11_after = child_happy_map[c11][g2]; hg11_after = gift_happy_map[c11][g2];
				hc12_after = child_happy_map[c12][g2]; hg12_after = gift_happy_map[c12][g2];
				hc21_after = child_happy_map[c21][g1]; hg21_after = gift_happy_map[c21][g1];
				hc22_after = child_happy_map[c22][g1]; hg22_after = gift_happy_map[c22][g1];
			
				child_happiness_points_diff = hc11_after  + hc12_after  + hc21_after  + hc22_after
											-(hc11_before + hc12_before + hc21_before + hc22_before);
				gift_happiness_points_diff  = hg11_after  + hg12_after  + hg21_after  + hg22_after
											-(hg11_before + hg12_before + hg21_before + hg22_before);

				ediff = cal_diff(total_child_happiness_points, total_gift_happiness_points,
					             child_happiness_points_diff, gift_happiness_points_diff);

				//#ifdef ZERO_TEMP
				if (ediff >= 0) { // perform swap
								  //#else
								  //            if(ediff > 0 || ran_generator() < exp(ediff/temperature)){ // perform swap
								  //#endif
					states[g2][gc21] = c11;
					states[g2][gc22] = c12;
					states[g1][gc11] = c21;
					states[g1][gc12] = c22;

					gift_list[c11][0] = g2; gift_list[c11][1] = gc21;
					gift_list[c12][0] = g2; gift_list[c12][1] = gc22;
					gift_list[c21][0] = g1; gift_list[c21][1] = gc11;
					gift_list[c22][0] = g1; gift_list[c22][1] = gc12;

					total_child_happiness_points += child_happiness_points_diff;
					total_gift_happiness_points += gift_happiness_points_diff;

					if (ediff > 0) {
						happiness += ediff;
						is_improved = true;
					}

					// Update unhappy_child list -------------------
					h11_before = eff_r * hc11_before + hg11_before;  // approx. happiness
					h12_before = eff_r * hc12_before + hg12_before;
					h21_before = eff_r * hc21_before + hg21_before;
					h22_before = eff_r * hc22_before + hg22_before;
					h11_after = eff_r * hc11_after + hg11_after;
					h12_after = eff_r * hc12_after + hg12_after;
					h21_after = eff_r * hc21_after + hg21_after;
					h22_after = eff_r * hc22_after + hg22_after;

					if (h11_before < unhappy_threshold && h11_after >= unhappy_threshold) {
						int u11 = 0;
						while (unhappy_child[u11] != c11) u11++;
						unhappy_child.erase(unhappy_child.begin() + u11);
					}
					if (h12_before < unhappy_threshold && h12_after >= unhappy_threshold) {
						int u12 = 0;
						while (unhappy_child[u12] != c12) u12++;
						unhappy_child.erase(unhappy_child.begin() + u12);
					}
					if (h21_before < unhappy_threshold && h21_after >= unhappy_threshold) {
						int u21 = 0;
						while (unhappy_child[u21] != c21) u21++;
						unhappy_child.erase(unhappy_child.begin() + u21);
					}
					if (h22_before < unhappy_threshold && h22_after >= unhappy_threshold) {
						int u22 = 0;
						while (unhappy_child[u22] != c22) u22++;
						unhappy_child.erase(unhappy_child.begin() + u22);
					}
					if (h11_before >= unhappy_threshold && h11_after < unhappy_threshold) unhappy_child.push_back(c11);
					if (h12_before >= unhappy_threshold && h12_after < unhappy_threshold) unhappy_child.push_back(c12);
					if (h21_before >= unhappy_threshold && h21_after < unhappy_threshold) unhappy_child.push_back(c21);
					if (h22_before >= unhappy_threshold && h22_after < unhappy_threshold) unhappy_child.push_back(c22);
				}
			}
		}
		// ==================================== 3v3 SWAP ====================================
		else if (is_triplet(c11)) {  // c11 == triplet; 3v111, 3v112, or 3v3
			get_triplets(c11, c12, c13);
			gc12 = gift_list[c12][1];
			gc13 = gift_list[c13][1];

			// GET CHILDREN FOR GIFT 2 (C21)
			gc21 = (int)(ran_generator()*Nqttys);  // get any c21
			c21 = states[g2][gc21];

			// GET C22 AND C23
			if (is_single(c21))	{
				bool found_c22 = false;
				
				// GET C22
				for (int j = 0; j < match_times; j++) {  // get a twin or single for c22
					gc22 = (int)(ran_generator()*(Nqttys - 1));
					if (gc22 >= gc21) gc22 += 1;  // not to take gc21
					c22 = states[g2][gc22];
					if (is_not_triplet(c22)) {  // get only twin or single
						found_c22 = true;
						break;
					}
				}
				
				// GET C23
				if (found_c22)	{
					if (is_single(c22)) {  //3v(1+1+1)
						for (int j = 0; j < match_times; j++) {  // get a single for c23
							gc23 = (int)(ran_generator()*Nqttys);
							c23 = states[g2][gc23];
							if (is_single(c23) && gc23 != gc21 && gc23 != gc22) {
								found_all_child = true;
								break;
							}
						}
					}
					else {  // c22 == twin; 3v(1+2)
						get_twins(c22, c23);  // get twin of c22 as c23
						gc23 = gift_list[c23][1];
						found_all_child = true;
					}
				}
			}
			else if (is_twin(c21)) { // c21 == twin; 3v(2+1)
				get_twins(c21, c22);  // get twin of c21 as c22
				gc22 = gift_list[c22][1];

				// GET C23
				for (int j = 0; j < match_times; j++) {  // get a single for c23
					gc23 = (int)(ran_generator()*Nqttys);
					c23 = states[g2][gc23];
					if (is_single(c23) && gc23 != gc21 && gc23 != gc22) {
						found_all_child = true;
						break;
					}
				}
			}
			else { // c21 == triplet; 3v3
				get_triplets(c21, c22, c23);  // get triplet of c21 as c22, c23
				gc22 = gift_list[c22][1];
				gc23 = gift_list[c23][1];
				found_all_child = true;
			}

			// COMPUTE DIFFERENCE AND SWAP
			if (found_all_child)
			{
				// check diff and perform swap if nec
				hc11_before = child_happy_map[c11][g1]; hg11_before = gift_happy_map[c11][g1];
				hc12_before = child_happy_map[c12][g1]; hg12_before = gift_happy_map[c12][g1];
				hc13_before = child_happy_map[c13][g1]; hg13_before = gift_happy_map[c13][g1];
				hc21_before = child_happy_map[c21][g2]; hg21_before = gift_happy_map[c21][g2];
				hc22_before = child_happy_map[c22][g2]; hg22_before = gift_happy_map[c22][g2];
				hc23_before = child_happy_map[c23][g2]; hg23_before = gift_happy_map[c23][g2];
				hc11_after = child_happy_map[c11][g2]; hg11_after = gift_happy_map[c11][g2];
				hc12_after = child_happy_map[c12][g2]; hg12_after = gift_happy_map[c12][g2];
				hc13_after = child_happy_map[c13][g2]; hg13_after = gift_happy_map[c13][g2];
				hc21_after = child_happy_map[c21][g1]; hg21_after = gift_happy_map[c21][g1];
				hc22_after = child_happy_map[c22][g1]; hg22_after = gift_happy_map[c22][g1];
				hc23_after = child_happy_map[c23][g1]; hg23_after = gift_happy_map[c23][g1];

				child_happiness_points_diff = hc11_after  + hc12_after  + hc13_after  + hc21_after  + hc22_after  + hc23_after
											-(hc11_before + hc12_before + hc13_before + hc21_before + hc22_before + hc23_before);
				gift_happiness_points_diff  = hg11_after  + hg12_after  + hg13_after  + hg21_after  + hg22_after  + hg23_after
											-(hg11_before + hg12_before + hg13_before + hg21_before + hg22_before + hg23_before);
			
				ediff = cal_diff(total_child_happiness_points, total_gift_happiness_points,
					             child_happiness_points_diff, gift_happiness_points_diff);
			
				//#ifdef ZERO_TEMP
				if (ediff >= 0) { // perform swap
								  //#else
								  //            if(ediff > 0 || ran_generator() < exp(ediff/temperature)){ // perform swap
								  //#endif
					states[g2][gc21] = c11;
					states[g2][gc22] = c12;
					states[g2][gc23] = c13;
					states[g1][gc11] = c21;
					states[g1][gc12] = c22;
					states[g1][gc13] = c23;

					gift_list[c11][0] = g2; gift_list[c11][1] = gc21;
					gift_list[c12][0] = g2; gift_list[c12][1] = gc22;
					gift_list[c13][0] = g2; gift_list[c13][1] = gc23;
					gift_list[c21][0] = g1; gift_list[c21][1] = gc11;
					gift_list[c22][0] = g1; gift_list[c22][1] = gc12;
					gift_list[c23][0] = g1; gift_list[c23][1] = gc13;

					total_child_happiness_points += child_happiness_points_diff;
					total_gift_happiness_points += gift_happiness_points_diff;

					if (ediff > 0) {
						happiness += ediff;
						is_improved = true;
					}

					// Update unhappy_child list -------------------
					h11_before = eff_r * hc11_before + hg11_before;  // approx. happiness
					h12_before = eff_r * hc12_before + hg12_before;
					h13_before = eff_r * hc13_before + hg13_before;
					h21_before = eff_r * hc21_before + hg21_before;
					h22_before = eff_r * hc22_before + hg22_before;
					h23_before = eff_r * hc23_before + hg23_before;
					h11_after = eff_r * hc11_after + hg11_after;
					h12_after = eff_r * hc12_after + hg12_after;
					h13_after = eff_r * hc13_after + hg13_after;
					h21_after = eff_r * hc21_after + hg21_after;
					h22_after = eff_r * hc22_after + hg22_after;
					h23_after = eff_r * hc23_after + hg23_after;

					if (h11_before < unhappy_threshold && h11_after >= unhappy_threshold) {
						int u11 = 0;
						while (unhappy_child[u11] != c11) u11++;
						unhappy_child.erase(unhappy_child.begin() + u11);
					}
					if (h12_before < unhappy_threshold && h12_after >= unhappy_threshold) {
						int u12 = 0;
						while (unhappy_child[u12] != c12) u12++;
						unhappy_child.erase(unhappy_child.begin() + u12);
					}
					if (h13_before < unhappy_threshold && h13_after >= unhappy_threshold) {
						int u13 = 0;
						while (unhappy_child[u13] != c13) u13++;
						unhappy_child.erase(unhappy_child.begin() + u13);
					}
					if (h21_before < unhappy_threshold && h21_after >= unhappy_threshold) {
						int u21 = 0;
						while (unhappy_child[u21] != c21) u21++;
						unhappy_child.erase(unhappy_child.begin() + u21);
					}
					if (h22_before < unhappy_threshold && h22_after >= unhappy_threshold) {
						int u22 = 0;
						while (unhappy_child[u22] != c22) u22++;
						unhappy_child.erase(unhappy_child.begin() + u22);
					}
					if (h23_before < unhappy_threshold && h23_after >= unhappy_threshold) {
						int u23 = 0;
						while (unhappy_child[u23] != c23) u23++;
						unhappy_child.erase(unhappy_child.begin() + u23);
					}
					if (h11_before >= unhappy_threshold && h11_after < unhappy_threshold) unhappy_child.push_back(c11);
					if (h12_before >= unhappy_threshold && h12_after < unhappy_threshold) unhappy_child.push_back(c12);
					if (h13_before >= unhappy_threshold && h13_after < unhappy_threshold) unhappy_child.push_back(c13);
					if (h21_before >= unhappy_threshold && h21_after < unhappy_threshold) unhappy_child.push_back(c21);
					if (h22_before >= unhappy_threshold && h22_after < unhappy_threshold) unhappy_child.push_back(c22);
					if (h23_before >= unhappy_threshold && h23_after < unhappy_threshold) unhappy_child.push_back(c23);
				}
			}
		}
		else {
			printf("error in c11= %d", c11);  // For debug only
		}

		// SHOW OR SAVE UPDATE RESULT
        if(is_improved) {
			timestep++;

			// CHECK IF GOT BETTER HAPPINESS
			if (happiness > output_threshold) {
				hdiff = happiness - best_happy;
				if (hdiff > 1e-20) {
					printf("%6ds - ", time(0) - t0cpu);
					printf("Improved by %.1e, ", hdiff);
					printf("now at %.10f", happiness);
					best_happy = happiness;
					
					if (h11_before < unhappy_threshold)	{
						unhappy_found++;
						printf(" by unhappy child search.");
					}

					printf("\n");
				}

				// OUTPUT DATA
				if (0 == timestep % step_out) {
					happiness = cal_all();
					printf("\nSteps: %lld, happiness: %.10f\n", timestep, happiness);
					cout << "Child, Gifts happiness points = " << total_child_happiness_points << ", " << total_gift_happiness_points << endl;
					write_conf(happiness);
					fflush(stdout);
				}
			}
		}

		mini_steps++;
	}

	// FINALIZING
	int tfscpu = time(0);
	int tcpu_total = tfscpu - t0cpu;
	cout << "\n### CPU time: " << tcpu_total << " secs ****" << endl;
	
	happiness = cal_all();
	printf("Total steps: %lld (%f / s)\n", timestep, 1.0 * timestep / tcpu_total);
	printf("Final happiness: %.10f\n", happiness);
	printf("Overall improvement: %.1e\n", happiness - init_happy);
	cout << "Child / Gifts happiness points = " << total_child_happiness_points 
		<< " / " << total_gift_happiness_points << endl;
	printf("Recommend effective ratio: %.5e\n", 
		child_gifts_happiness_ratio_3 * pow((double) total_child_happiness_points / total_gift_happiness_points, 2));
#ifdef WEAK_SEARCH_MAIN
	printf("Weak search rate: %.5f", 1.0 - 1.0 / non_weak_search);
#else
	printf("Weak search rate: %.5f,  ", 1.0 / weak_search);
#endif // WEAK_SEARCH_MAIN
	printf("Unhappy found rate: %.5f\n", (double)unhappy_found / timestep);
	cout << "unhappy_child length = " << unhappy_child.size() << endl;

	int n_diff = 0;
	for (int c = 0; c < Nchild; c++) {
		if (gift_list[c][0] != answ_init[c]) n_diff++;
	}
	printf("N of swap children: %d    [ %f/s ]", n_diff, ((double)n_diff) / tcpu_total);

	write_conf(happiness, false);
    write_conf(happiness, true);
	
	int tfcpu = time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;
		
	cin.get();
	return 0;
}


//Debug
/*if (states[g1][gc11] != c11) printf("c11 not consistent\n");
if (states[g1][gc12] != c12) printf("c12 not consistent\n");
if (states[g1][gc13] != c13) printf("c13 not consistent\n");
if (states[g2][gc21] != c21) printf("c21 not consistent\n");
if (states[g2][gc22] != c22) printf("c22 not consistent\n");
if (states[g2][gc23] != c23) printf("c23 not consistent\n");

if (gift_list[c11][0] != g1) printf("gl c11 not consistent\n");
if (gift_list[c12][0] != g1) printf("gl c12 not consistent\n");
if (gift_list[c13][0] != g1) printf("gl c13 not consistent\n");
if (gift_list[c21][0] != g2) printf("gl c21 not consistent\n");
if (gift_list[c22][0] != g2) printf("gl c22 not consistent\n");
if (gift_list[c23][0] != g2) printf("gl c23 not consistent\n");*/
