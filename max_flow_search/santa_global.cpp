#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <unordered_map>
#include <algorithm> 

#include "santa_par.h"
#include "santa_global.h"

using namespace std;

////////// GLOBAL VARIABLES //////////
vector<vector<short>> wlist_child;  // Child's wishlist of gifts
vector<vector<int>> wlist_gifts;  // Gift's wishlist of children
vector<vector<int>> states(Ntypes, vector<int>(Nqttys, -1));  // Child of each gift assigned to [Took 4MB for int]
vector<short> answ(Nchild, -1); // gift of each child now have (best case). [Took 4MB for int]

#ifndef FIX_EFF_R
Eff_R eff_r;     // Effective happiness ratio r^3 C^2 / G^2,
Eff_R happy_min; // Minimum of approximate happiness points
#endif // !FIX_EFF_R


////////// GLOBAL FUNCTIONS //////////

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

void init(){
	srand(time(NULL));

    // READ WISH LISTS
	cout << "Reading " << name_wl_child;
	read_csv(name_wl_child, wlist_child, 0, 1);
	cout << " completed." << endl;
	cout << "Reading " << name_wl_gifts;
	read_csv(name_wl_gifts, wlist_gifts, 0, 1);
	cout << " completed." << endl;
	
	// READ INIT CONF AND INIT STATES
	cout << "Reading " << name_initCF; 
	vector<vector<int>> answ_ini;
	read_csv(name_initCF, answ_ini, 1, 0);

	vector<int> states_fill_i(Ntypes, 0); // Filling position of each gift
	for (int i = 0; i < Nchild; i++) { // Build states and answ
		int c = answ_ini[i][0]; // Child id
		int g = answ_ini[i][1]; // Gift id
		answ[c] = g;
		states[g][states_fill_i[g]] = c;
		states_fill_i[g]++;
	}
	cout << " completed. " << "[" << answ_ini.size() << "x" << answ_ini[0].size() << "]" << endl;
	
    // CHECK SIZES
	cout << "Check sizes..." << endl;
    if(wlist_child.size() != Nchild) error("(init) child wish list size: (size Nchild)", 2, wlist_child.size(), Nchild);
	if(wlist_gifts.size() != Ntypes) error("(init) gifts wish list size: (size Nchild)", 2, wlist_gifts.size(), Ntypes);
	if (states.size() != Ntypes) error("(init) states number error: (size Nchild)", 2, states.size(), Ntypes);

    for(int i=0; i<Nchild; i++){ // wlist_child
        if(wlist_child[i].size() != Nchprf) error("(init) : child wish list prf size (size Nchprf)", 2, wlist_child[i].size(), Nchprf);
    }
	for(int i=0; i<Ntypes; i++){ // wlist_gifts
        if(wlist_gifts[i].size() != Ngfprf) error("(init) : gifts wish list prf size (size Nchprf)", 2, wlist_gifts[i].size(), Ngfprf);
    }
	
	// Check triplet and twin initial assignment consistent
	cout << "Total triplet inconsistent N:      " << triplet_consistent(answ) << endl;
	cout << "Total twin inconsistent N:         " << twin_consistent(answ) << endl;
	check_gift_count(answ);
	check_states();
	cout << "Order inconsistent N:              " << triplet_twin_order() << endl;	
}

Eff_R update_eff_r(long long int total_child_happiness_points, long long int total_gift_happiness_points) {
	Eff_R eff_r_cal = llround(child_gifts_happiness_ratio_3 * pow((double)total_child_happiness_points
		/ total_gift_happiness_points, 2));
#if defined(FIX_EFF_R)
	printf("Recommend effective ratio:      %12lld\n", eff_r_cal);
#else
	eff_r = eff_r_cal;  // Update effective ratio
	printf("New effective ratio:            %12lld\n", eff_r);
#endif // !ZERO_GIFT_R
	return eff_r_cal;
}

void reset_happiness_map(vector<vector<HappyInt>>& happiness, HappyInt non) {
	for (int i = 0; i < happiness.size(); i++) {
		for (int j = 0; j < happiness[0].size(); j++) {
			happiness[i][j] = non;
		}
	}
}

// Build Happiness Maps
void build_happines_all(vector<vector<HappyInt>>& happiness_c, vector<vector<HappyInt>>& happiness_g) {
	cout << "Build approx. happiness maps...";

	cout << " triplets";
	build_app_triplet_happiness_map(happiness_c);
	cout << " twins";
	build_app_twin_happiness_map(happiness_c);
	cout << " singles";
	build_app_single_happiness_map(happiness_c);
	cout << " gifts";
	build_app_gift_happiness_map(happiness_g);
	cout << " Done!" << endl;
}

void build_app_triplet_happiness_map(vector<vector<HappyInt>>& happiness_c) {
	for (int c = 0; c < Ntripl; c += 3) {
		for (int t = 0; t < 3; t++) {
			for (int j = 0; j < Nchprf; j++) {
				int g = wlist_child[c + t][j];
				HappyInt happy_points = triplet_r * ((Nchprf - j) * points_childIN - points_childOUT);
				happiness_c[c][g] += happy_points;
				happiness_c[c + 1][g] += happy_points;
				happiness_c[c + 2][g] += happy_points;
			}
		}
	}
}

void build_app_twin_happiness_map(vector<vector<HappyInt>>& happiness_c) {
	for (int c = Ntripl; c < Ntwins; c += 2) {
		for (int t = 0; t < 2; t++) {
			for (int j = 0; j < Nchprf; j++) {
				int g = wlist_child[c + t][j];
				HappyInt happy_points = twin_r * ((Nchprf - j) * points_childIN - points_childOUT);
				happiness_c[c    ][g] += happy_points;
				happiness_c[c + 1][g] += happy_points;
			}
		}
	}
}

void build_app_single_happiness_map(vector<vector<HappyInt>>& happiness_c) {
	for (int c = Ntwins; c < Nchild; c += 1) {
		for (int j = 0; j < Nchprf; j++) {
			int g = wlist_child[c][j];
			HappyInt happy_points = single_r * ((Nchprf - j) * points_childIN - points_childOUT);
			happiness_c[c][g] += happy_points;
		}
	}
}

void build_app_gift_happiness_map(vector<vector<HappyInt>>& happiness_g) {
	for (int g = 0; g < Ntypes; g += 1) {
		for (int j = 0; j < Ngfprf; j++) {
			int c = wlist_gifts[g][j];
			if (c >= Ntwins) {       // single
				HappyInt happy_points = gift_r * ((Ngfprf - j) * points_giftsIN - points_giftsOUT);
				happiness_g[c][g] += happy_points;
			}
			else if (c >= Ntripl) {  // twin
				int c2 = -1;
				get_twins(c, c2);
				HappyInt happy_points = gift_twin_r * ((Ngfprf - j) * points_giftsIN - points_giftsOUT);
				happiness_g[c][g] += happy_points;
				happiness_g[c2][g] += happy_points;
			}
			else {                   // triplet
				int c2 = -1, c3 = -1;
				get_triplets(c, c2, c3);
				HappyInt happy_points = gift_triplet_r * ((Ngfprf - j) * points_giftsIN - points_giftsOUT);
				happiness_g[c][g] += happy_points;
				happiness_g[c2][g] += happy_points;
				happiness_g[c3][g] += happy_points;
			}
		}
	}
}

// Happiness points averaged between triplets and twins scaled by single_r
inline HappyInt cal_single_child(int c, int g, const vector<vector<HappyInt>>& happiness_c) {
	return happiness_c[c][g];
}

inline HappyInt cal_single_gifts(int c, int g, const vector<vector<HappyInt>>& happiness_g) {
	return happiness_g[c][g];
}

// Happiness points averaged between triplets and twins
long long int cal_child_all(const vector<short>& answ, const vector<vector<HappyInt>>& happiness_c) {
	long long int total_child_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_child_happiness += happiness_c[c][answ[c]];
	}
	return llround(total_child_happiness / (double)single_r);
}

long long int cal_gift_all(const vector<short>& answ, const vector<vector<HappyInt>>& happiness_g) {
	long long int total_gifts_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_gifts_happiness += happiness_g[c][answ[c]];
	}
	return llround(total_gifts_happiness / (double)single_r);
}

// Calculate happiness given total child and gift happiness points
// C: total child happiness points
// G: total gift happiness points
double cal_all_from_points(long long int C, long long int G) {
	return (pow(C * 1.0 * child_gifts_happiness_ratio, 3)
		+ pow(G, 3)) / pow(1.0 * max_gifts_happiness, 3);
}

//double cal_all(const vector<short>& answ, const vector<vector<HappyInt>>& happiness_c, const vector<vector<HappyInt>>& happiness_g) {
//	long long int total_child_happiness = 0;
//	long long int total_gifts_happiness = 0;
//
//	for (int c = 0; c < answ.size(); c++) { // i: ichild
//		total_child_happiness += happiness_c[c][answ[c]];
//		total_gifts_happiness += happiness_g[c][answ[c]];
//	}
//
//	return (pow(total_child_happiness * 1.0 * child_gifts_happiness_ratio, 3) 
//		+ pow(total_gifts_happiness, 3)) / pow(1.0 * max_gifts_happiness, 3);
//}


// Original happienss points without averaged
int cal_single_child_ori(int c, int g) {
	for (int j = 0; j < Nchprf; j++) {
		if (wlist_child[c][j] == g) return (Nchprf - j) * points_childIN;
	}
	return points_childOUT;
}

int cal_single_gifts_ori(int c, int g) {
	for (int j = 0; j < Ngfprf; j++) {
		if (wlist_gifts[g][j] == c) return (Ngfprf - j) * points_giftsIN;
	}
	return points_giftsOUT;
}

long long int cal_child_all_ori(const vector<short>& answ) {
	long long int total_child_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_child_happiness += cal_single_child_ori(c, answ[c]);
	}
	return total_child_happiness;
}

long long int cal_gift_all_ori(const vector<short>& answ) {
	long long int total_gifts_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_gifts_happiness += cal_single_gifts_ori(c, answ[c]);
	}
	return total_gifts_happiness;
}

double cal_all_ori(const vector<short>& answ) {
	long long int total_child_happiness = 0;
	long long int total_gifts_happiness = 0;

	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_child_happiness += cal_single_child_ori(c, answ[c]);
		total_gifts_happiness += cal_single_gifts_ori(c, answ[c]);
	}

	return (pow(total_child_happiness * 1.0 * child_gifts_happiness_ratio, 3)
		+ pow(total_gifts_happiness, 3)) / pow(1.0 * max_gifts_happiness, 3);
}


// Update answ from states (global)
void states_to_answ() {
	for (int g = 0; g < Ntypes; g++) {
		for (int j = 0; j < Nqttys; j++) {
			answ[states[g][j]] = g;
		}
	}
}

// Update states from answ
void answ_to_states(const vector<int>& answ) {
	vector<int> states_fill_i(Ntypes, 0); // Filling position of each gift
	for (int c = 0; c < Nchild; c++) { // Build states and answ
		int g = answ[c]; // Gift id
		states[g][states_fill_i[g]] = c;
		states_fill_i[g]++;
	}
}


// Check states with non-negative childId
int check_states() {
	int n_negative = 0;
	for (int g = 0; g < Ntypes; g++) {
		for (int ic = 0; ic < Nqttys; ic++) {
			if (states[g][ic] < 0) {
				n_negative++;
			}
		}
	}
	cout << "Negative childId in states:        " << n_negative << endl;
	return n_negative;
}

// Check the order of triplets and twins in the states
int triplet_twin_order() {
	int n_inconsistent = 0;
	for (int g = 0; g < Ntypes; g++) {
		int j = 0;
		while (j < Nqttys) {
			int c = states[g][j];
			if (c < Ntripl) {
				if (j + 2 >= Nqttys || states[g][j + 1] != c + 1 || states[g][j + 2] != c + 2) {
					//cout << g << ", " << j << " - " << c << " triplet order inconsistent" << endl;
					n_inconsistent++;
				}
				else j += 2;
			}
			else if (c < Ntwins) {
				if (j + 1 >= Nqttys || states[g][j + 1] != c + 1) {
					//cout << g << ", " << j << " - " << c << " twin order inconsistent" << endl;
					n_inconsistent++;
				}
				else j += 1;
			}
			j++;
		}
	}
	return n_inconsistent;
}

int triplet_consistent(const vector<short>& answ) {
	int n_inconsistent = 0;
	for (int c = 0; c < Ntripl; c += 3) {
		if (answ[c] != answ[c + 1] || answ[c] != answ[c + 2]) {
			n_inconsistent++;
			//cout << c << " is not consistent: " << answ[c] << " " << answ[c + 1] << " " << answ[c + 2] << endl;
		}
	}
	return n_inconsistent;
}

int twin_consistent(const vector<short>& answ) {
	int n_inconsistent = 0;
	for (int c = Ntripl; c < Ntwins; c += 2) {
		if (answ[c] != answ[c + 1]) {
			n_inconsistent++;
			//cout << c << " is not consistent: " << answ[c] << " " << answ[c + 1] << endl;
		}
	}
	return n_inconsistent;
}

int check_negative(const vector<short>& answ) {
	int n_negative = 0;
	for (int i = 0; i < answ.size(); i++) {
		if (answ[i] < 0) n_negative++;
	}
	return n_negative;
}

// Check gift quantity in answ
bool check_gift_count(const vector<short>& answ) {
	vector<int> gift_count(Ntypes, 0);
	for (int i = 0; i < answ.size(); i++) {
		gift_count[answ[i]]++;
	}
	int min = *min_element(gift_count.begin(), gift_count.end());
	int max = *max_element(gift_count.begin(), gift_count.end());
	cout << "Gift count min, max:               " << min << ", " << max << endl;
	if (min == Nqttys && max == Nqttys) return true;
	else return false;
}

// Check remaining unassigned gift
int check_unassigned_gifts(const vector<int>& unassigned_gifts) {
	int remaining_gifts = 0;
	cout << "Unassigned gifts: ";
	for (int i = 0; i < unassigned_gifts.size(); i++) {
		if (unassigned_gifts[i] != 0) {
			cout << i << ":" << unassigned_gifts[i] << " ";
			remaining_gifts += unassigned_gifts[i];
		}
	}
	cout << "[" << remaining_gifts << "]" << endl;
	return remaining_gifts;
}

//bool verify_app_happy_points(vector<short>& answ, vector<vector<HappyInt>>& happiness,
//	long long int C, long long int G) {
//	long long int total_points_from_CG = (C - Nchild * points_childOUT) * single_r + (G - Ntypes * Nqttys*points_giftsOUT) * gift_r;
//	long long int total_points_from_happiness = 0;
//	for (int c = 0; c < answ.size(); c++) {
//		int g = answ[c];
//		total_points_from_happiness += happiness[c][g];
//	}
//	cout << "App points from C and G (shift)  = " << total_points_from_CG << endl;
//	cout << "App points from happiness matrix = " << total_points_from_happiness << endl;
//	cout << "App points difference            = " << total_points_from_happiness - total_points_from_CG << endl;
//	return total_points_from_happiness == total_points_from_CG;
//}

// Compare difference between two answ
int compare_answ_diff(const vector<short>& answ1, const vector<short>& answ2) {
	int n_diff = 0;  // N of swap child
	for (int c = 0; c < answ1.size(); c++) {
		if (answ1[c] != answ2[c]) n_diff++;
	}
	return n_diff;
}

// Print vector
void print_vector(vector<int>& vec) {
	for (int i = 0; i < vec.size(); i++) {
		cout << vec[i] << " ";
	}
	cout << endl;
}

void error(string errinfo, int nnum, double num1, double num2) {
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events

	cout << "\nError: " << errinfo;
	switch (nnum) {
	case 0:  cout << endl; break;
	case 1:  cout << ": " << num1 << endl; break;
	case 2:  cout << ": " << num1 << " " << num2 << endl; break;
	default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
#ifdef WINDOWS
	cin.get();
#endif // WINDOWS
	exit(1);
}

void read_csv(const char* filename, vector<vector<int>>& matrix, int header, int index_col) {
	/* header determine the number of header rows to skip (1 = skip first row, 0 = no skip) */
	/* index_col determine the number of index columns to skip (1 = skip first column, 0 = no skip) */

	ifstream file;
	file.open(filename);
	if (!file.is_open()) {
		string error_s = "Can't read "; error_s.append(filename);
		error(error_s, 0);
	}
	string line;
	int i = 0;
	while (getline(file, line)) {  // read till '\n'
		if (i >= header) {
			int j = 0;
			vector<int> row;
			istringstream templine(line);  // convert string to stream
			string data;
			while (getline(templine, data, ',')) {  // read till ','. Note ',' at end of line will not be included.
				if (j >= index_col) {
					row.push_back(atoi(data.c_str())); // Missing data will be treated as zero
				}
				j++;
			}
			matrix.push_back(row);
		}
		i++;
	}
	file.close();
}

void read_csv(const char* filename, vector<vector<short>>& matrix, int header, int index_col) {
	/* header determine the number of header rows to skip (1 = skip first row, 0 = no skip) */
	/* index_col determine the number of index columns to skip (1 = skip first column, 0 = no skip) */

	ifstream file;
	file.open(filename);
	if (!file.is_open()) {
		string error_s = "Can't read "; error_s.append(filename);
		error(error_s, 0);
	}
	string line;
	int i = 0;
	while (getline(file, line)) {  // read till '\n'
		if (i >= header) {
			int j = 0;
			vector<short> row;
			istringstream templine(line);  // convert string to stream
			string data;
			while (getline(templine, data, ',')) {  // read till ','. Note ',' at end of line will not be included.
				if (j >= index_col) {
					row.push_back(atoi(data.c_str())); // Missing data will be treated as zero
				}
				j++;
			}
			matrix.push_back(row);
		}
		i++;
	}
	file.close();
}

// Output result from answ with tag
void write_conf(vector<short>& answ, const char* tag) {
    ofstream CONF;
	char filepath[100];

	sprintf_s(filepath, "%s_%s.csv", output_path, tag);
	CONF.open(filepath);

    // Print out
    CONF << "ChildId,GiftId" << endl;
    for (int i=0; i < Nchild; i++) {
        CONF << i << "," << answ[i] << endl;
    }
    CONF.close();
}

// Output result from states with score
void write_conf(vector<short>& answ, double best_happy, bool isrestart) {
	static int Nout = 0;
	Nout++;

	ofstream CONF;
	char file_name[100];

	if (isrestart) {
		sprintf_s(file_name, "%s.csv", output_path);
		CONF.open(file_name);
	}
	else {
		sprintf_s(file_name, "%s_optimum_%02d_%8d.csv", output_path, Nout, (int)(best_happy * 100000000));
		CONF.open(file_name);
	}

	// Make ichild-igift list & check if quantity correct
	//int count = 0;
	//vector <int> check_gifts(Ntypes, 0);
	//vector <int> pred(Nchild); // ichild:igift
	//for (int i = 0; i<Ntypes; i++) {
	//	for (int j = 0; j<Nqttys; j++) {
	//		pred[states[i][j]] = i;

	//		count++;
	//		check_gifts[i] ++;
	//	}
	//}

	//if (count != Nchild) error("(write_conf) Nchild inconsist", 2, count, Nchild);
	//for (auto n : check_gifts) { if (n != Nqttys) error("(write_conf) gift Nqttys inconsist", 2, n, Nqttys); }

	// Print out
	CONF << "ChildId,GiftId" << endl;
	for (int i = 0; i<Nchild; i++) {
		CONF << i << "," << answ[i] << endl;
	}

	CONF.close();
}
