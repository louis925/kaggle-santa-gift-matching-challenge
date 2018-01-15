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
	//static bool first= true;
	//if(first){
	//	srand(time(NULL));
	//	first= false;
	//}
	return rand() / ((double)RAND_MAX + 1.0);
}

inline int cal_single_child(int ichild, int igift) {
	if (child_happy_map[ichild].find(igift) != child_happy_map[ichild].end())
		return child_happy_map[ichild][igift];
	else
		return points_childOUT;
}

inline int cal_single_gifts(int ichild, int igift) {
	if (gift_happy_map[igift].find(ichild) != gift_happy_map[igift].end())
		return gift_happy_map[igift][ichild];
	else
		return points_giftsOUT;
}

// Sum of single child-gift happiness without normalization
inline int raw_single_happiness(int ichild, int igift) {
	return cal_single_child(ichild, igift) * child_gifts_happiness_ratio + cal_single_gifts(ichild, igift);
}

inline bool is_single(int ichild) { return ichild >= Ntwins; }
inline bool is_not_single(int ichild) { return ichild < Ntwins; }
inline bool is_twin(int ichild) { return ichild < Ntwins && ichild >= Ntripl; }
inline bool is_triplet(int ichild) { return ichild < Ntripl; }
inline bool is_not_triplet(int ichild) { return ichild >= Ntripl; }
inline void get_twins(int ichild, int* twin) {
	if (ichild % 2 == 0) *twin = ichild - 1;  // Note the twin now start at odd number id
	else                 *twin = ichild + 1;
}
inline void get_triplets(int ichild, int* triplet1, int* triplet2) {
	switch (ichild % 3)
	{
	case 0:	*triplet1 = ichild + 1; *triplet2 = ichild + 2; break;
	case 1: *triplet1 = ichild - 1; *triplet2 = ichild + 1; break;
	case 2: *triplet1 = ichild - 2; *triplet2 = ichild - 1; break;
	default: break;
	}
}

//************************************************************************************************************
int main(){
    long long int timestep= 0;
    double happiness;
	int t0cpu = time(0);

    cout << "********** Metropolis MC simulation on Santa problem **********\n" << endl;
	cout << "########## Initializing System... ##########" << endl;
	double best_happy=  init();
    happiness= best_happy;
    cout << "Initialization completed. The initial old happiness = " << happiness << endl;
	//cout << "weak_search = " << weak_search << endl;
	//cout << "non_weak_search = " << non_weak_search << endl;

	//cout << "\n########## The Simulation Begins !! ##########" << endl;
	//printf("Steps: %lld, happiness: %.10f\n", timestep, happiness);
	

    // FINALIZING
	double final_new_happiness = cal_new_all();  // new score for v2 of the competition
	
	printf("Final new happiness: %.15f\n", final_new_happiness);

	long long int total_child_happiness = cal_child_all();
	long long int total_gift_happiness = cal_gift_all();
	printf("total child happiness points: %lld\n", total_child_happiness);
	printf("total gift  happiness points: %lld\n", total_gift_happiness);
	double effective_ratio = pow(child_gifts_happiness_ratio, 3) * pow((double)total_child_happiness, 2) 
		/ pow((double)total_gift_happiness, 2) ;
	printf("effective first order ratio: %.1f\n", effective_ratio);

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	cin.get();

	return 0;
}
