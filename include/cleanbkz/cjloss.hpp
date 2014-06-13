/*
    Copyright (C) 2014 Janos Follath.
 
    This file is part of cleanbkz.

    cleanbkz is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cleanbkz is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cleanbkz.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CJLOSS_H_
#define CJLOSS_H_

#include<NTL/mat_ZZ.h>
#include<NTL/vec_RR.h>

using namespace NTL;

/** 	@class
 	@brief class representing a knapsack problem with the purpose to solve it with enumeration and be able to check if the result is correct */

class cjloss {
	private:
		vec_ZZ solution;
		vec_ZZ values;
		ZZ sum;
		void randomize(long, ZZ);
		bool check();
	public:
		/** Generates (with the given seed) a random knapsack problem with maxlength bit long numbers and whose corresponding lattice has dimension dimension and density= dimension/maxlength */
		cjloss(long dimension, double density, long seed); 

		/** Computes the density of the generated lattice. */ 
		double get_density() const;

		/** Given the shortest vector in the scaled cjloss lattice, it computes and prints the solution to the corresponding subset sum problem 
		void print_solution(const vec_RR& shortest);*/

		/** Puts the knapsack problem in question to the given output stream. */
		friend std::ostream& operator<<(std::ostream&, const cjloss&);

		/** Returns the reduced and truncated basis of the cjloss lattice */
		mat_ZZ get_basis(int k); 
};


#endif
