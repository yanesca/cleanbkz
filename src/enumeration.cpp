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

#include <iostream>
#include "enumeration.hpp"

using namespace std;

void enumerate(mat_ZZ& basis, double* pruning, vec_RR& result) {
	
	BKZ_QP1(basis, 0.99, 30);		

	//cout << "BKZ30 reduced basis: "  << endl  << basis << endl;

	mat_RR mu;
	vec_RR c;

	ComputeGS(basis,mu,c);

	//cout << "GS coeffitients: "  << endl  << mu << endl;
	//cout << "GS basis lengths: "  << c << endl;

	result.SetLength(basis.NumRows());
	for(int i= 0; i < result.length(); i++)
		result[i]= 0;

	cout << "short vector: "  << result << endl;
} 
