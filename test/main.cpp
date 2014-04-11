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
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/enumeration.hpp>
#include <cleanbkz/boundtools.hpp>

using namespace std;
using namespace NTL;

int main(int argc, char** argv) {
  	int dimension= 0;
	//double density= 0.94;
	int method= 0; 
	int bkz_blocksize= 2;

  	if(argc<2)
		return 1;  

	stringstream ss(argv[1]);
	ss >> dimension;

	if (argc==3) {
		stringstream ss2(argv[2]);
		ss2 >> method;
		}

	if (argc==4) {
		stringstream ss3(argv[3]);
		ss3 >> bkz_blocksize;
		}

  	/*cjloss m(dimension,density);

	//cout << "Original knapsack problem : "  << endl << m << endl;

	double* Rvec= new double[dimension];
	for(int i= 0; i < dimension; i++)
		Rvec[i]= dimension;

	vec_RR result;
	if(method==0)	
		enumerate_epr(m.basis,bkz_blocksize,Rvec,result);
	else
		enumerate_ntl(m.basis,bkz_blocksize,Rvec,result);

	m.print_solution(result);*/

	double t_node, t_reduc;
	measure_epr(dimension, 2, 10000, 1000000, t_node, t_reduc); 
	cout << "t_node= " << t_node << endl;
	cout << "t_reduc= " << t_reduc << endl;

	return 0;
}
