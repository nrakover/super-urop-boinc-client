/*
 * timeslice.cpp
 *
 *  Created on: Feb 22, 2014
 *      Author: nrakover
 */

#include "timeslice.h"

namespace alfax {

Timeslice::Timeslice(vector<int> observations) :
		observations_(observations), hidden_state_(-1) {
}

Timeslice::Timeslice(int hidden_state, vector<int> observations) :
		observations_(observations), hidden_state_(hidden_state) {
}

}  // namespace alfax
