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
vector<vector<int>> wlist_child;  // Child's wishlist of gifts
#ifndef NO_GIFT_MAP
vector<vector<int>> wlist_gifts;  // Gift's wishlist of children
#endif // !NO_GIFT_MAP
vector<vector<int>> states(Ntypes, vector<int>(Nqttys, -1));  // Child of each gift assigned to [Took 4MB for int]
vector<int> answ(Nchild, -1); // gift of each child now have. [Took 4MB for int]
//vector<vector<int>> gift_list(Nchild, vector<int>(2, -1)); // gift of each child now have. coordinate in states {(giftId, index in gift), ...} [Took 8MB for int]

vector<vector<HappyInt>> happiness(Nchild, vector<HappyInt>(Ntypes, 0));  // Map of happyness given gift for each child. happiness[child][gift] == approx_happiness 
																	   // [Took 8GB for long long int as HappyInt]

#ifdef ZERO_GIFT_R

#else // !ZERO_GIFT_R
HappyInt eff_r = eff_r_init;  // effective happiness ratio r^3 C^2 / G^2, 
HappyInt single_r = eff_r * gift_r;
#endif // ZERO_GIFT_R

#ifdef MERGED
HappyInt triplet_r = single_r;
HappyInt twin_r = single_r;
#else
HappyInt triplet_r = single_r / 3;
HappyInt twin_r = single_r / 2;
#endif // MERGED

////////// GLOBAL FUNCTIONS //////////

void error(string errinfo, int nnum, double num1, double num2){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: " << errinfo;
	switch(nnum){
		case 0:  cout << endl; break;
		case 1:  cout << ": " << num1 << endl; break;
		case 2:  cout << ": " << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
	exit(1);
}

void read_csv(const char* filename, vector<vector<int>>& matrix, int header, int index_col){
	/* header determine the number of header rows to skip (1 = skip first row, 0 = no skip) */
	/* index_col determine the number of index columns to skip (1 = skip first column, 0 = no skip) */

	fstream file;
	file.open(filename);
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

double init(){
	srand(time(NULL));

    // READ WISH LISTS
    read_csv(name_wl_child, wlist_child, 0, 1);
	cout << "Reading " << name_wl_child << " completed." << endl;
#ifndef NO_GIFT_MAP
	read_csv(name_wl_gifts, wlist_gifts, 0, 1);
	cout << "Reading " << name_wl_gifts << " completed." << endl;
#endif // !NO_GIFT_MAP
	
	// READ INIT CONF AND INIT STATES
	vector<vector<int>> answ_ini;
	read_csv(name_initCF, answ_ini, 1, 0);

	vector<int> states_fill_i(Ntypes, 0); // Filling position of each gift
	for (int i = 0; i < Nchild; i++) { // Build states and answ
		int c = answ_ini[i][0]; // Child id
		int g = answ_ini[i][1]; // Gift id
		answ[c] = g;
		//answ[c][1] = states_fill_i[g];
		states[g][states_fill_i[g]] = c;
		states_fill_i[g]++;
	}
	cout << "Reading " << name_initCF << " completed. " << "[" << answ_ini.size() << "x" << answ_ini[0].size() << "]" << endl;


    // CHECK SIZES
	cout << "Check sizes..." << endl;
    if(wlist_child.size() != Nchild) error("(init) child wish list size: (size Nchild)", 2, wlist_child.size(), Nchild);
#ifndef NO_GIFT_MAP
	if(wlist_gifts.size() != Ntypes) error("(init) gifts wish list size: (size Nchild)", 2, wlist_gifts.size(), Ntypes);
#endif // !NO_GIFT_MAP
	if (states.size() != Ntypes) error("(init) states number error: (size Nchild)", 2, states.size(), Ntypes);

    for(int i=0; i<Nchild; i++){ // wlist_child
        if(wlist_child[i].size() != Nchprf) error("(init) : child wish list prf size (size Nchprf)", 2, wlist_child[i].size(), Nchprf);
    }
#ifndef NO_GIFT_MAP
	for(int i=0; i<Ntypes; i++){ // wlist_gifts
        if(wlist_gifts[i].size() != Ngfprf) error("(init) : gifts wish list prf size (size Nchprf)", 2, wlist_gifts[i].size(), Ngfprf);
    }
#endif // !NO_GIFT_MAP
	
	// Check triplet and twin initial assignment consistent
	cout << "Total triplet inconsistent N: " << triplet_consistent(answ) << endl;
	cout << "Total twin inconsistent N: " << twin_consistent(answ) << endl;
	check_gift_count(answ);
	cout << "Order inconsistent N: " << triplet_twin_order() << endl;

	// Build Happiness Maps
	build_happines_all();

	return cal_all(answ);
}

// Build Happiness Maps
void build_happines_all() {
	cout << "Build approx. happiness maps...";

	reset_happiness_map(0);  // Fill happiness map with 0
	cout << " triplets";
	build_app_triplet_happiness_map();
	cout << " twins";
	build_app_twin_happiness_map();
#ifndef NO_SINGLE_MAP
	cout << " singles";
	build_app_single_happiness_map();
#endif // !NO_SINGLE_MAP
#ifndef NO_GIFT_MAP
	cout << " gifts";
	build_app_gift_happiness_map();
#endif // !NO_GIFT_MAP
	cout << " Done!" << endl;
}

void reset_happiness_map(HappyInt non) {
	for (int i = 0; i < happiness.size(); i++) {
		for (int j = 0; j < happiness[0].size(); j++) {
			happiness[i][j] = non;
		}
	}
}

int triplet_consistent(vector<int> &answ) {
	int n_inconsistent = 0;
	for (int c = 0; c < Ntripl; c += 3) {
		if (answ[c] != answ[c + 1] || answ[c] != answ[c + 2]) {
			n_inconsistent++;
			cout << c << " is not consistent: " << answ[c] << " " << answ[c + 1] << " " << answ[c + 2] << endl;
		}
	}
	return n_inconsistent;
}

int twin_consistent(vector<int> &answ) {
	int n_inconsistent = 0;
	for (int c = Ntripl; c < Ntwins; c += 2) {
		if (answ[c] != answ[c + 1]) {
			n_inconsistent++;
			cout << c << " is not consistent: " << answ[c] << " " << answ[c + 1] << endl;
		}
	}
	return n_inconsistent;
}

// Check the order of triplets and twins in the states
int triplet_twin_order() {
	int n_inconsistent = 0;
	for (int g = 0; g < Ntypes; g++) {
		int j = 0;
		while (j < Nqttys) {
			int c = states[g][j];
			if (c < Ntripl) {
				if (j + 2 >= Nqttys || states[g][j + 1] != c+1 || states[g][j + 2] != c+2) {
					cout << g << ", " << j << " - " << c << " triplet order inconsistent" << endl;
					n_inconsistent++;
				}
				else j += 2;
			}
			else if (c < Ntwins) {
				if (j + 1 >= Nqttys || states[g][j + 1] != c+1) {
					cout << g << ", " << j << " - " << c << " twin order inconsistent" << endl;
					n_inconsistent++;
				}
				else j += 1;
			}
			j++;
		}		
	}
	return n_inconsistent;
}

// Check gift quantity in answ
bool check_gift_count(vector<int> &answ) {
	vector<int> gift_count(Ntypes, 0);
	for (int i = 0; i < answ.size(); i++) {
		gift_count[answ[i]]++;
	}
	int min = *min_element(gift_count.begin(), gift_count.end());
	int max = *max_element(gift_count.begin(), gift_count.end());
	cout << "Gift count min, max: " << min << ", " << max << endl;
	if (min == Nqttys && max == Nqttys) return true;
	else return false;
}

void build_app_triplet_happiness_map() {
	for (int c = 0; c < Ntripl; c += 3) {
		for (int t = 0; t < 3; t++) {
			for (int j = 0; j < Nchprf; j++) {
				int g = wlist_child[c + t][j];
				HappyInt happy_points = triplet_r * (-1 * points_childOUT + (Nchprf - j) * points_childIN);
				happiness[c    ][g] += happy_points;
#ifndef MERGED
				happiness[c + 1][g] += happy_points;
				happiness[c + 2][g] += happy_points;
#endif // !MERGED
			}
		}
	}
}

void build_app_twin_happiness_map() {
	for (int c = Ntripl; c < Ntwins; c += 2) {
		for (int t = 0; t < 2; t++) {
			for (int j = 0; j < Nchprf; j++) {
				int g = wlist_child[c + t][j];
				HappyInt happy_points = twin_r * (-1 * points_childOUT + (Nchprf - j) * points_childIN);
				happiness[c    ][g] += happy_points;
#ifndef MERGED
				happiness[c + 1][g] += happy_points;
#endif // !MERGED
			}
		}
	}
}

void update_eff_r(HappyInt total_child_happiness_points, HappyInt total_gift_happiness_points) {
#ifndef ZERO_GIFT_R
	eff_r = child_gifts_happiness_ratio_3 * pow((double)total_child_happiness_points
		/ total_gift_happiness_points, 2);  // Update effective ratio
	single_r = eff_r * gift_r;
#ifdef MERGED
	triplet_r = single_r;
	twin_r = single_r;
#else
	triplet_r = single_r / 3;
	twin_r = single_r / 2;
#endif // MERGED
	printf("New effective ratio: %d\n", eff_r);
#else
	printf("Recommend effective ratio: %.7e\n", child_gifts_happiness_ratio_3 * pow((double)total_child_happiness_points
		/ total_gift_happiness_points, 2))
#endif // !ZERO_GIFT_R
}

void build_app_single_happiness_map() {
	for (int c = Ntwins; c < Nchild; c += 1) {
		for (int j = 0; j < Nchprf; j++) {
			int g = wlist_child[c][j];
			HappyInt happy_points = single_r * (-1 * points_childOUT + (Nchprf - j) * points_childIN);
			happiness[c][g] += happy_points;
		}
	}
}

#ifndef NO_GIFT_MAP
void build_app_gift_happiness_map() {
	for (int g = 0; g < Ntypes; g += 1) {
		for (int j = 0; j < Ngfprf; j++) {
			int c = wlist_gifts[g][j];
			if (c >= Ntwins) {
				HappyInt happy_points = gift_r * (-1 * points_giftsOUT + (Ngfprf - j) * points_giftsIN);
				happiness[c    ][g] += happy_points;
			}
			else if (c >= Ntripl) {
				HappyInt happy_points = gift_twin_r * (-1 * points_giftsOUT + (Ngfprf - j) * points_giftsIN);
				happiness[c    ][g] += happy_points;
#ifndef MERGED
				happiness[c + 1][g] += happy_points;
#endif // !MERGED
			}
			else {
				HappyInt happy_points = gift_triplet_r * (-1 * points_giftsOUT + (Ngfprf - j) * points_giftsIN);
				happiness[c    ][g] += happy_points;
#ifndef MERGED
				happiness[c + 1][g] += happy_points;
				happiness[c + 2][g] += happy_points;
#endif // !MERGED
			}
		}
	}
}
#endif // !NO_GIFT_MAP


inline int cal_single_child(int c, int g) {
	for (int j = 0; j < Nchprf; j++) {
		if (wlist_child[c][j] == g) return (Nchprf - j) * points_childIN;
	}
	return points_childOUT;
}

#ifndef NO_GIFT_MAP
inline int cal_single_gifts(int c, int g) {
	for (int j = 0; j < Ngfprf; j++) {
		if (wlist_gifts[g][j] == c) return (Ngfprf - j) * points_giftsIN;
	}
	return points_giftsOUT;
}

// Sum of single child-gift happiness without normalization
inline long long int raw_single_happiness(int ichild, int igift) {
	return cal_single_child(ichild, igift) * single_r + cal_single_gifts(ichild, igift);
}

long long int cal_child_all(vector<int> &answ) {
	long long int total_child_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_child_happiness += cal_single_child(c, answ[c]);
	}
	return total_child_happiness;
}

long long int cal_gift_all(vector<int> &answ) {
	long long int total_gifts_happiness = 0;
	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_gifts_happiness += cal_single_gifts(c, answ[c]);
	}
	return total_gifts_happiness;
}

// New score for stage 2
double cal_all(vector<int> &answ) {
	long long int total_child_happiness = 0;
	long long int total_gifts_happiness = 0;

	for (int c = 0; c < answ.size(); c++) { // i: ichild
		total_child_happiness += cal_single_child(c, answ[c]);
		total_gifts_happiness += cal_single_gifts(c, answ[c]);
	}

	return (pow(total_child_happiness * 1.0 * child_gifts_happiness_ratio, 3) 
		+ pow(total_gifts_happiness, 3)) / pow(1.0 * max_gifts_happiness, 3);
}
#endif // !NO_GIFT_MAP

// Update answ from states
void states_to_answ() {
	for (int g = 0; g < Ntypes; g++) {
		for (int j = 0; j < Nqttys; j++) {
			answ[states[g][j]] = g;
		}
	}
}

// Update states from answ
void answ_to_states() {
	vector<int> states_fill_i(Ntypes, 0); // Filling position of each gift
	for (int c = 0; c < Nchild; c++) { // Build states and answ
		int g = answ[c]; // Gift id
		states[g][states_fill_i[g]] = c;
		states_fill_i[g]++;
	}
}

// Output result from answ with tag
void write_conf(vector<int> &answ, const char* tag) {
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
void write_conf(double best_happy, bool isrestart) {
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
	int count = 0;
	vector <int> check_gifts(Ntypes, 0);
	vector <int> pred(Nchild); // ichild:igift
	for (int i = 0; i<Ntypes; i++) {
		for (int j = 0; j<Nqttys; j++) {
			pred[states[i][j]] = i;

			count++;
			check_gifts[i] ++;
		}
	}

	if (count != Nchild) error("(write_conf) Nchild inconsist", 2, count, Nchild);
	for (auto n : check_gifts) { if (n != Nqttys) error("(write_conf) gift Nqttys inconsist", 2, n, Nqttys); }

	// Print out
	CONF << "ChildId,GiftId" << endl;
	for (int i = 0; i<Nchild; i++) {
		CONF << i << "," << pred[i] << endl;
	}

	CONF.close();
}
