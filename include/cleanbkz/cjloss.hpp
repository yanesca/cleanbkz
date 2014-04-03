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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CJLOSS_H_
#define CJLOSS_H_

#include<NTL/mat_ZZ.h>
#include<NTL/vec_RR.h>

using namespace NTL;

/** Class representing a knapsack problem with the purpose to solve it with enumeration and be able to check if the result is correct */
class cjloss {
	private:
		mat_ZZ basis;
		vec_ZZ solution;
		vec_ZZ values;
		ZZ sum;
		void randomize(long, ZZ);
		bool check();
	public:
		cjloss(long,double);
		mat_ZZ get_basis();
		double get_density() const;
		void print_solution(const vec_RR& shortest);
		friend std::ostream& operator<<(std::ostream&, const cjloss&);
};


#endif
