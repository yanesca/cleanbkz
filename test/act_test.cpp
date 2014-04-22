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

#include <iostream>
#include <sstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <NTL/LLL.h>

using namespace std;

int main(int argc, char** argv) {
	cjloss l(10, 0.94, 0);
	BKZ_QP1(l.basis, 0.99, 2);		
	cout << "Lattice basis: " << endl << l.basis << endl << endl;

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);

	for(int i= 0; i < mu1.NumRows(); i++)
		mu1[i][i]= 1;

	cout << "GS coeffs: " << endl << mu1 << endl;

	cout << "Squared lengths of GS basis vectors: " << endl << c1 << endl;
}
