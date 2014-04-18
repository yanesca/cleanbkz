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

double polytope_volume(double vols[], double bounds[], int dim);

int main(int argc, char** argv) {

	int dim= 80;
	int beta= 2;
	double t_node= 3.47193e-08; 
	double t_reduc= 0.101471;
	cjloss m(dim, 0.94, 0);
	double delta= 1e-1;
	unsigned long iterations= 1000;
	
	stringstream s;
	s << argv[1];	
	s >> iterations;
	s << argv[2];	
	s >> delta;

	BKZ_QP1(m.basis, 0.99, beta);		

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(m.basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
			conv(c[i], c1[i]);

	double* boundary= new double[dim];	

	generate_boundary(c, t_node, t_reduc, dim, boundary, dim-1, delta, iterations); 

	cout << endl << endl << "# " << dim << " dimensional boundary function: " << endl;
	for(int i= 0; i < dim; i++)
		cout << i << " " << boundary[i] << endl;
	cout << endl;

	return 0;
}
