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

#include <sstream>
#include <iostream>
#include <NTL/LLL.h>
#include <NTL/mat_RR.h>
#include "cjloss.hpp"
#include "enumeration.hpp"

using namespace std;
using namespace NTL;

int main(int argc, char** argv) {
  	int dimension= 0;
	double density= 0.94;

  	if(argc<2)
		return 1;  

	stringstream ss(argv[1]);
	ss >> dimension;

  	cjloss m(dimension,density);
	mat_ZZ basis= m.get_basis();

	cout << "Original knapsack problem : "  << endl << m << endl;

	vec_RR result;
	enumerate(basis,NULL,result);

	m.print_solution(result);

	return 0;
}
