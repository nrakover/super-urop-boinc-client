/*
 * utils.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: nrakover
 */

#include "utils.h"

#include <vector>
#include <cmath>

using namespace std;

namespace alfax {

void InitializeA_counts(vector<vector<long double> > *A_counts, int N) {
	A_counts->resize(N);
	for (int i = 0; i < N; i++) {
		A_counts->at(i).resize(N);
		for (int j = 0; j < N; j++)
			A_counts->at(i).at(j) = 0.0;
	}
}

void InitializeB_counts(vector<vector<vector<long double> > > *B_counts,
		int N, vector<int> *Sigma) {
	B_counts->resize(N);
	for (int i = 0; i < N; i++) {
		B_counts->at(i).resize(Sigma->size());
		for (int j = 0; j < Sigma->size(); j++) {
			B_counts->at(i).at(j).resize(Sigma->at(j));
			for (int k = 0; k < Sigma->at(j); k++) {
				B_counts->at(i).at(j).at(k) = 0.0;
			}
		}
	}
}

void Initialize_memo_table(vector<vector<long double> > *memo_table,
		int N, int timeslices) {
	memo_table->resize(N);
	for (int i = 0; i < N; i++) {
		memo_table->at(i).resize(timeslices);
		for (int j = 0; j < timeslices; j++)
			memo_table->at(i).at(j) = -1.0;
	}
}

double LogSumExp(vector<double> *sum_terms) {
	double max = *max_element(sum_terms->begin(), sum_terms->end());

	double sum = 0.0;
	for (int i = 0; i < sum_terms->size(); i++) {
		sum += exp(sum_terms->at(i) - max);
	}

	return max + log(sum);
}

} // namespace alfax
