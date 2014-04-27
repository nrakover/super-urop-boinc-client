//============================================================================
// Name        : HMM_Client.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "client.h"
#include "timeslice.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "boinc_api.h"
#include "filesys.h"

using namespace std;

vector<vector<alfax::Timeslice> > GetData(const string&);

int main(int argc, char **argv) {
	boinc_init();

	vector<vector<double> > A;
	vector<vector<vector<double> > > B;
	vector<int> Sig;

	int N = alfax::Client::LoadParams(&A, &B, &Sig);

	vector<vector<alfax::Timeslice> > data = GetData("HMM_EM_data");

	alfax::Client *client = new alfax::Client(N, Sig, A, B, data);

	string operation = argv[1];

	if (operation.compare("ESTEP") == 0) {

		client->PartialE();

	} else if (operation.compare("LOGLIKELIHOOD") == 0) {

		client->PartialLL();

	}

	boinc_finish(0);
}

vector<vector<alfax::Timeslice> > GetData(const string &data_filepath){

	vector<vector<alfax::Timeslice> > data_points;
	ifstream myfile(data_filepath.c_str());
	if (myfile.is_open()) {
			string line;
			int count = 0;
			int i=0;
			while(getline(myfile, line)){
				int observation_count = 0;
				stringstream ss(line);
				string observation;
				vector<int> observations;
				while(getline(ss, observation, ';')){
					istringstream iss(observation);
					int n;
					iss >> n;

					observations.push_back(n);
					observation_count++;
				}
				alfax::Timeslice t = alfax::Timeslice(observations);
				if (data_points.size()==i){
					vector<alfax::Timeslice> timeslices;
					timeslices.push_back(t);
					data_points.push_back(timeslices);
				}
				else{
					data_points.at(i).push_back(t);
				}
				count++;
				if (count == 15){
					i++;
					count=0;
				}
			}
		}
	return data_points;
}
