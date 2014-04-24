/*
 * timeslice.h
 *
 *  Created on: Feb 22, 2014
 *      Author: nrakover
 */

#include <vector>
using namespace std;

#ifndef TIMESLICE_H_
#define TIMESLICE_H_

namespace alfax {

class Timeslice {
public:
	Timeslice(vector<int> observations);

	Timeslice(int hidden_state, vector<int> observations);

	vector<int> observations_;
	int hidden_state_;
};

}  // namespace alfax


#endif /* TIMESLICE_H_ */
