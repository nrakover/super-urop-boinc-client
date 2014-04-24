/*
 * utils.h
 *
 *  Created on: Feb 26, 2014
 *      Author: nrakover
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <vector>

using namespace std;

namespace alfax {

void InitializeA_counts(vector<vector<long double> > *A_counts, int N);

void InitializeB_counts(vector<vector<vector<long double> > > *B_counts,
		int N, vector<int> *Sigma);

void Initialize_memo_table(vector<vector<long double> > *memo_table,
		int N, int timeslices);

double LogSumExp(vector<double> *sum_terms);

} // namespace alfax


#endif /* UTILS_H_ */
