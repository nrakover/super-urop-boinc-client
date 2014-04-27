/*
 * client.h
 *
 *  Created on: Feb 22, 2014
 *      Author: nrakover
 */

#include <vector>
#include "timeslice.h"

#ifndef CLIENT_H_
#define CLIENT_H_

namespace alfax {

class Client {
public:
	Client(int N, vector<int> Sig, vector<vector<double> > A,
			vector<vector<vector<double> > > B,
			vector<vector<Timeslice> > data);

	void PartialE();

	void PartialLL();

	void ExportECounts(vector<vector<long double> > *A_counts,
			vector<vector<vector<long double> > > *B_counts);

	void ExportLL(const string &directory_path, long double log_likelihood);

	static int LoadParams(vector<vector<double> > *A,
			vector<vector<vector<double> > > *B, vector<int> *Sig);

	long double PosteriorGivenObs(int j, int y_j, int y_j_plus_1,
			vector<Timeslice> *observation_sequence,
			vector<vector<long double> > *Alpha_memo_table,
			vector<vector<long double> > *Beta_memo_table);

	long double LogTotalProbOfObs(vector<Timeslice> *observation_sequence,
			vector<vector<long double> > *Alpha_memo_table,
			vector<vector<long double> > *Beta_memo_table);

	long double LogEmissionProb(int y, const Timeslice *observations);

	long double LogAlpha(int p, int j, vector<vector<long double> > *memo_table,
			vector<Timeslice> *sequence);

	long double LogBeta(int p, int j, vector<vector<long double> > *memo_table,
			vector<Timeslice> *sequence);

	vector<vector<double> > A_;
	vector<vector<vector<double> > > B_;
	int N_;
	vector<int> Sigma_;
	vector<vector<Timeslice> > data_;
};

} // namespace alfax

#endif /* CLIENT_H_ */
