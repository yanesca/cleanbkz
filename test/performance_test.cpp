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
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/boundary.hpp>
#include <NTL/LLL.h>
#include <sstream>

using namespace std;

extern RR t_extreme_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n);

extern RR fact_RR(int x);

extern RR ball_vol_RR(int k, RR r);

extern RR RR_PI;

extern RR* factorials;

extern void init_factorials(int up_to);

void rand_boundary(int dim, RR* boundary, double R);

int main(int argc, char** argv) {
	stringstream ss;
	int count= 10;
	long dim= 80;

	if (argc>=2) {
		ss << argv[1];
		ss >> dim;	
		ss.clear();		
		}

	if (argc>=3) {
		ss << argv[2];
		ss >> count;	
		ss.clear();		
		}

	init_factorials(2*dim+2);

	cjloss lattice(dim, 0.94, 2330);

	mat_RR mu1;
        vec_RR c1;
        ComputeGS(lattice.get_basis(2),mu1,c1);

        //lengths of the GS basis vectors
        RR* c= new RR[mu1.NumRows()];
        for(int i= 0; i < mu1.NumRows(); i++) {
        	c[i]= SqrRoot(c1[i]);
        	c[i].SetPrecision(RR_PRECISION);
        }
        
	RR* rand= new RR[dim];

	using namespace std;
  	clock_t begin = clock();
	for(int i= 0; i< count; i++) {
		rand_boundary(dim, rand, sqrt(dim));
		t_extreme_RR(rand, c, 7.35e-08, 47.603, dim);
	}

	clock_t end = clock();
  	cout << ((end - begin) / CLOCKS_PER_SEC) << endl;

	return 0;
}

void rand_boundary(int dim, RR* boundary, double R) {
	RR tmp;

	for(int i= 0; i< dim; i++) 
		random(boundary[i]);

	for (int c = 0 ; c < ( dim - 1 ); c++) 
    		for (int d = 0 ; d < dim - c - 1; d++) 
      			if (boundary[d] > boundary[d+1]) {
        			tmp = boundary[d];
        			boundary[d] = boundary[d+1];
        			boundary[d+1] = tmp;
      			}

	for(int i= 0; i< dim-1; i++) { 
		boundary[i]= boundary[i+1]*=R;		
		i++;
	}
	boundary[dim-1]= R;
	if(dim%2==0)
		boundary[dim-2]= R;
		
/*	for(int i= 0; i< dim; i++) 
		cout << boundary[i] << " ";
	cout << endl;*/
}

