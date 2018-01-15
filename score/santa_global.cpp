#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <unordered_map>

#include "santa_global.h"
#include "santa_par.h"

using namespace std;

////////// GLOBAL VARIABLES //////////
vector<vector<int>> wlist_child;
vector<vector<int>> wlist_gifts;
vector<vector<int>> states(Ntypes);
vector<int> unhappy_child;  // Child with poor happiness {childId1, ...}

vector<vector<int>> gift_list(Nchild); // gift of each child now have. coordinate in states {(giftId, index in gift), ...}

vector<unordered_map<int, int>> child_happy_map(Nchild);  // Map of happyness given gift for each child. child_happy_map[child][gift] == child_happiness
vector<unordered_map<int, int>> gift_happy_map(Ntypes); // Map of happyness given child for each gift. gift_happy_map[gift][child] == gift_happiness

////////// GLOBAL FUNCTIONS //////////

//inline double ran_generator(){
//	return rand()/((double) RAND_MAX+1.0);
//}

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
	while (getline(file, line))  // read till '\n'
	{
		if (i >= header) 
		{
			int j = 0;
			vector<int> row;
			istringstream templine(line);  // convert string to stream
			string data;
			while (getline(templine, data, ','))  // read till ','. Note ',' at end of line will not be included.
			{
				if (j >= index_col)
				{
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
    read_csv(name_wl_gifts, wlist_gifts, 0, 1);
    cout << "Reading " << name_wl_child << " and " << name_wl_gifts << " completed." << endl;

    // READ INIT CONF AND INIT STATES
    vector<vector<int>> conf;
	read_csv(name_initCF, conf, 1, 0);

	for(int i=0; i<Nchild; i++){ // build states
	    states[conf[i][1]].push_back(conf[i][0]);
		vector<int> row = {conf[i][1], (int)states[conf[i][1]].size() - 1};
		gift_list[conf[i][0]] = row;
    }
    cout << "Reading " << name_initCF << " completed." << endl;

    // CHECK SIZES
	cout << "Check sizes..." << endl;
    if(wlist_child.size() != Nchild) error("(init) child wish list size: (size Nchild)", 2, wlist_child.size(), Nchild);
    if(wlist_gifts.size() != Ntypes) error("(init) gifts wish list size: (size Nchild)", 2, wlist_gifts.size(), Ntypes);
    if(states.size() != Ntypes) error("(init) states number error: (size Nchild)", 2, states.size(), Ntypes);
    
    for(int i=0; i<Nchild; i++){ // wlist_child
        if(wlist_child[i].size() != Nchprf) error("(init) : child wish list prf size (size Nchprf)", 2, wlist_child[i].size(), Nchprf);
    }
    for(int i=0; i<Ntypes; i++){ // wlist_gifts
        if(wlist_gifts[i].size() != Ngfprf) error("(init) : gifts wish list prf size (size Nchprf)", 2, wlist_gifts[i].size(), Ngfprf);
    }
    for(int i=0; i<Ntypes; i++){ // states
        if(states[i].size() != Nqttys) error("(init) initial conf not has right numbers: (size Ntypes)", 2, states[i].size(), Ntypes);
        /*for(int j=0; j<Ngfprf; j++){
            if(states[i][j]<Ntwins) ptr_twins[states[i][j]]= j;
        }*/
    }

	// Build Happiness Maps
	cout << "Build happiness maps..." << endl;
	build_child_happiness_map();
	build_gift_happiness_map();

    // CALCULATE HAPPINESS
	cout << "Calculate total happiness..." << endl;
    double happiness = cal_all();

	//// Add unhappy child into unhappy_child list
	//cout << "Build unhappy child list" << endl;
	//for (int ichild = 0; ichild < Nchild; ichild++)	{
	//	int igift = gift_list[ichild][0];
	//	//if (cal_single_child_gifts(ichild, igift) < 0) {
	//	if (raw_single_happiness(ichild, igift) < 0) {
	//		unhappy_child.push_back(ichild);
	//	}
	//}
    return happiness;
}

void build_child_happiness_map() {
	for (int ichild = 0; ichild < Nchild; ichild++)	{
		for (int i = 0; i < Nchprf; i++) {
			child_happy_map[ichild][wlist_child[ichild][i]] = points_childIN * (Nchprf - i);
		}
	}
}

void build_gift_happiness_map() {
	for (int igift = 0; igift < Ntypes; igift++) {
		for (int i = 0; i < Ngfprf; i++) {
			gift_happy_map[igift][wlist_gifts[igift][i]] = points_giftsIN * (Ngfprf - i);
		}
	}
}

//int cal_single_child(int ichild, int igift){
//    for(int i=0; i<Nchprf; i++){ 
//		if(wlist_child[ichild][i]==igift) return points_childIN * (Nchprf-i);
//	}
//	return points_childOUT;
//}
//
//int cal_single_gifts(int ichild, int igift){
//    for(int i=0; i<Ngfprf; i++){
//		if(wlist_gifts[igift][i]==ichild) return points_giftsIN * (Ngfprf-i);
//	}
//	return points_giftsOUT;
//}

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

double cal_single_child_gifts(int ichild, int igift) {
	return (cal_single_child(ichild, igift) * child_gifts_happiness_ratio + cal_single_gifts(ichild, igift)) * 1.0 / max_gifts_happiness;
}

// Sum of single child-gift happiness without normalization
inline int raw_single_happiness(int ichild, int igift) {
	return cal_single_child(ichild, igift) * child_gifts_happiness_ratio 
		   + cal_single_gifts(ichild, igift);
}

long long int cal_child_all() {
	long long int total_child_happiness = 0;
	for (int i = 0; i<Ntypes; i++) { // i: igift
		for (int j = 0; j<Nqttys; j++) { // states[i][j]: ichild
			total_child_happiness += cal_single_child(states[i][j], i);
		}
	}
	return total_child_happiness;
}

long long int cal_gift_all() {
	long long int total_gifts_happiness = 0;
	for (int i = 0; i<Ntypes; i++) { // i: igift
		for (int j = 0; j<Nqttys; j++) { // states[i][j]: ichild
			total_gifts_happiness += cal_single_gifts(states[i][j], i);
		}
	}
	return total_gifts_happiness;
}

double cal_all(){
	long long int total_child_happiness = 0;
    long long int total_gifts_happiness = 0;

    for(int i=0; i<Ntypes; i++){ // i: igift
        for(int j=0; j<Nqttys; j++){ // states[i][j]: ichild
            total_child_happiness += cal_single_child(states[i][j], i);
            total_gifts_happiness += cal_single_gifts(states[i][j], i);
        }
    }

	return (total_child_happiness * child_gifts_happiness_ratio + total_gifts_happiness) * 1.0 / max_gifts_happiness;
}


// New score for stage 2
double cal_new_all() {
	long long int total_child_happiness = 0;
	long long int total_gifts_happiness = 0;

	for (int i = 0; i<Ntypes; i++) { // i: igift
		for (int j = 0; j<Nqttys; j++) { // states[i][j]: ichild
			total_child_happiness += cal_single_child(states[i][j], i);
			total_gifts_happiness += cal_single_gifts(states[i][j], i);
		}
	}

	return (pow(total_child_happiness * 1.0 * child_gifts_happiness_ratio, 3) 
		+ pow(total_gifts_happiness, 3)) / pow(1.0 * max_gifts_happiness, 3);
}

void write_conf(double best_happy, bool isrestart){
    static int Nout= 0;
    Nout ++;

    ofstream CONF;
	char file_name[100];

	if (isrestart) {
		sprintf_s(file_name, "%s.csv", output_path);
		CONF.open(file_name);
	}
    else{
        sprintf_s(file_name, "%s_optimum_%02d_%8d.csv", output_path, Nout, (int) (best_happy*100000000));
        CONF.open(file_name);
    }

    // Make ichild-igift list & check if quantity correct
    int count= 0;
    vector <int> check_gifts(Ntypes, 0);
    vector <int> pred(Nchild); // ichild:igift
    for(int i=0; i<Ntypes; i++){
        for(int j=0; j<Nqttys; j++){
            pred[states[i][j]]= i;
            
            count ++;
            check_gifts[i] ++;
        }
    }
    
    if(count != Nchild) error("(write_conf) Nchild inconsist", 2, count, Nchild);
    for(auto n:check_gifts){ if(n != Nqttys) error("(write_conf) gift Nqttys inconsist", 2, n, Nqttys);}

    // Print out
    CONF << "ChildId,GiftId" << endl;
    for(int i=0; i<Nchild; i++){
        CONF << i << "," << pred[i] << endl;
    }

    CONF.close();
}
