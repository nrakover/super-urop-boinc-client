/*
 * client.cpp
 *
 *  Created on: Feb 22, 2014
 *      Author: nrakover
 */

#include "client.h"
#include "utils.h"

#include <cmath>
#include <fstream>

namespace alfax {

Client::Client(int N, vector<int> Sig, vector<vector<double> > A,
		vector<vector<vector<double> > > B, vector<vector<Timeslice> > data) :
		A_(A), B_(B), N_(N), Sigma_(Sig), data_(data) {
}

void Client::PartialE() {
	vector<vector<long double> > A_counts;
	vector<vector<vector<long double> > > B_counts;

	InitializeA_counts(&A_counts, N_);  // Initialize counts matrix A_counts
	InitializeB_counts(&B_counts, N_, &Sigma_); // Initialize counts matrix B_counts

	for (int data_index = 0; data_index < data_.size(); data_index++) { // iterate over data points
		vector<vector<long double> > Alpha_memo_table;
		Initialize_memo_table(&Alpha_memo_table, N_,
				data_.at(data_index).size()); //Initialize memoization table
		vector<vector<long double> > Beta_memo_table;
		Initialize_memo_table(&Beta_memo_table, N_,
				data_.at(data_index).size()); //Initialize memoization table

		vector<Timeslice> observation_sequence = data_.at(data_index);

		for (int observation_index = 0;
				observation_index < observation_sequence.size() - 1;
				observation_index++) {				// iterate over timeslices
			for (int p = 1; p < N_; p++) {// iterate over current state			// N_-1
				for (int q = 1; q < N_; q++) {// iterate over next state				// N_-1

					long double posterior_prob = PosteriorGivenObs(
							observation_index, p, q, &observation_sequence,
							&Alpha_memo_table, &Beta_memo_table);

					A_counts.at(p).at(q) += posterior_prob;
					for (int o = 0; o < Sigma_.size(); o++) {
						B_counts.at(p).at(o).at(
								observation_sequence.at(observation_index).observations_.at(
										o)) += posterior_prob;
					}
					if (observation_index == 0)
						A_counts.at(0).at(p) += posterior_prob;
//					else if (observation_index
//							== observation_sequence.size() - 2)
//						A_counts.at(q).at(N_ - 1) += posterior_prob;
				}
			}
		}
	}

	ExportECounts(&A_counts, &B_counts);
}

void Client::PartialLL() {
	long double log_likelihood = 0.0L;
		for (int data_index = 0; data_index < data_.size(); data_index++) {
			vector<Timeslice> observation_sequence = data_.at(data_index);
			vector<vector<long double> > Alpha_memo_table;
			Initialize_memo_table(&Alpha_memo_table, N_,
					observation_sequence.size());
			vector<vector<long double> > Beta_memo_table;
			Initialize_memo_table(&Beta_memo_table, N_,
					observation_sequence.size());

			log_likelihood += LogTotalProbOfObs(&observation_sequence,
					&Alpha_memo_table, &Beta_memo_table);
		}

		ExportLL("<file name>", log_likelihood);	//TODO: need to specify file
}

void Client::ExportECounts(vector<vector<long double> > *A_counts,
		vector<vector<vector<long double> > > *B_counts) {

	ofstream output_file;
	output_file.open("transitions.txt");

	for (int i = 0; i < N_; i++) {
		for (int j = 0; j < N_; j++) {
			output_file << A_counts->at(i).at(j) << ";";
		}
		output_file << "\n";
	}

	output_file.close();
	output_file.clear();

	output_file.open("emissions.txt");

	for (int i = 1; i < N_; i++) {		// N_-1
		for (int j = 0; j < Sigma_.size(); j++) {
			for (int k = 0; k < Sigma_.at(j); k++) {
				output_file << B_counts->at(i).at(j).at(k) << ";";
			}
			output_file << "|";
		}
		output_file << "\n";
	}

	output_file.close();
}

void Client::ExportLL(const string &dest_file, long double log_likelihood) {

	ofstream f;
	f.open(dest_file.c_str());
	f << log_likelihood;
	f.close();
}

int Client::LoadParams(const string &params_path, vector<vector<double> > *A,
		vector<vector<vector<double> > > *B, vector<int> *Sig) {

	string bar = "";
	if (params_path.substr(params_path.length() - 1,1).compare("/") != 0 )
		bar = "/";
	string transitions_path = params_path + bar + "transitions.txt";
	string emissions_path = params_path + bar + "emissions.txt";
	string line;

	ifstream params_file;
	params_file.open(transitions_path.c_str());	// First read transition matrix

	int N = 0;
	int start = 0;
	A->resize(1);
	std::getline(params_file, line);
	for (int i = 0; i < line.length(); i++) {
		if (line.substr(i, 1).compare(";") == 0) {
			N++;
			A->at(0).push_back(atof(line.substr(start, i - start).c_str()));
			start = i + 1;
		}
	}
	A->resize(N);

	int row = 1;
	while (std::getline(params_file, line)) {
		if (row > N - 1)
			break;
		int start = 0;
		for (int i = 0; i < line.length(); i++) {
			if (line.substr(i, 1).compare(";") == 0) {
				A->at(row).push_back(
						atof(line.substr(start, i - start).c_str()));
				start = i + 1;
			}
		}
		row++;
	}
	params_file.close();

	params_file.open(emissions_path.c_str());// Next read the emissions matrix

	B->resize(N);
	row = 1;
	while (std::getline(params_file, line)) {
		if (row > N - 1)
			break;
		B->at(row).resize(1);
		start = 0;
		int obs_support = 0;
		int obs_index = 0;
		for (int i = 0; i < line.length(); i++) {
			if (line.substr(i, 1).compare(";") == 0) {
				B->at(row).at(obs_index).push_back(
						atof(line.substr(start, i - start).c_str()));
				obs_support++;
				start = i + 1;
			} else if (line.substr(i, 1).compare("|") == 0) {
				if (row == 1)
					Sig->push_back(obs_support);
				obs_index++;
				B->at(row).resize(obs_index + 1);
				obs_support = 0;
				start = i + 1;
			}
		}
		B->at(row).resize(obs_index);
		row++;
	}

	return N;
}

long double Client::PosteriorGivenObs(int j, int y_j, int y_j_plus_1,
		vector<Timeslice> *observation_sequence,
		vector<vector<long double> > *Alpha_memo_table,
		vector<vector<long double> > *Beta_memo_table) {
	long double Z = LogTotalProbOfObs(observation_sequence, Alpha_memo_table,
			Beta_memo_table);

	double a = LogAlpha(y_j, j, Alpha_memo_table, observation_sequence);
	double b = log(A_[y_j][y_j_plus_1]);
	double c = LogEmissionProb(y_j, &observation_sequence->at(j));
	double d = LogBeta(y_j_plus_1, j + 1, Beta_memo_table,
			observation_sequence);

	double log_likelihood = a + b + c + d - Z;

	return exp(log_likelihood);
}

long double Client::LogTotalProbOfObs(vector<Timeslice> *observation_sequence,
		vector<vector<long double> > *Alpha_memo_table,
		vector<vector<long double> > *Beta_memo_table) {

	vector<double> sum_terms;
	sum_terms.resize(N_ - 1);				// N_-2
	for (int p = 1; p < N_; p++) {			// N_-1
		sum_terms.at(p - 1) = LogAlpha(p, 0, Alpha_memo_table,
				observation_sequence)
				+ LogBeta(p, 0, Beta_memo_table, observation_sequence);
	}
	return LogSumExp(&sum_terms);
}

long double Client::LogEmissionProb(int y, const Timeslice *x) {
	double b = 0.0;
	for (int i = 0; i < x->observations_.size(); i++) {
		b += log(B_.at(y).at(i).at(x->observations_.at(i)));
	}
	return b;
}

long double Client::LogAlpha(int p, int j,
		vector<vector<long double> > *memo_table, vector<Timeslice> *sequence) {
	if (j == 0)
		return log(A_[0][p]);

	if (memo_table->at(p).at(j) != -1.0)
		return memo_table->at(p).at(j);

	vector<double> sum_terms;
	sum_terms.resize(N_ - 1);					// N_-2

	for (int q = 1; q < N_; q++) {			// N_-1
		long double a = log(A_[q][p]);
		long double b = LogEmissionProb(q, &sequence->at(j - 1));
		long double c = LogAlpha(q, j - 1, memo_table, sequence);

		sum_terms.at(q - 1) = a + b + c;
	}
	double log_sum_exp = LogSumExp(&sum_terms);
	memo_table->at(p).at(j) = log_sum_exp;

	return log_sum_exp;
}

long double Client::LogBeta(int p, int j,
		vector<vector<long double> > *memo_table, vector<Timeslice> *sequence) {

	if (memo_table->at(p).at(j) != -1.0) {
		return memo_table->at(p).at(j);
	}

	long double b = LogEmissionProb(p, &sequence->at(j));

	if (j == (sequence->size() - 1)) {
		//return log(A_[p][N_ - 1]) + b;
		return b;
	}

	vector<double> sum_terms;
	sum_terms.resize(N_ - 1);			// N_-2

	for (int q = 1; q < N_; q++) {			// N_-1
		long double a = log(A_[p][q]);

		long double c = LogBeta(q, j + 1, memo_table, sequence);

		sum_terms.at(q - 1) = a + b + c;
	}
	double log_sum_exp = LogSumExp(&sum_terms);
	memo_table->at(p).at(j) = log_sum_exp;

	return log_sum_exp;
}

} // namespace alfax
