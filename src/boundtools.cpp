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

#include <NTL/LLL.h>
#include <iostream>
#include <cleanbkz/boundtools.hpp>
#include <cleanbkz/cjloss.hpp>

using namespace std;
using namespace NTL;

//extern 
void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, const unsigned long termination, double &time); 

void measure_epr(int dimension, int beta, unsigned long samples, unsigned long termination, double& t_node, double& t_reduc) {
	double enu_time= 0;
	cjloss* lattice;
	mat_RR mu1;
	vec_RR c1;
	clock_t begin, end;
	vec_RR result;

	t_node= t_reduc= 0;

	double** mu= new double*[dimension];
	for(int i= 0; i < dimension; i++)
		mu[i]= new double[dimension]; 

	double* c= new double[dimension];

	double* Rvec= new double[dimension];
	for(int i= 0; i < dimension; i++)
		Rvec[i]= dimension;

	for(unsigned long i= 0; i < samples; i++) {	
		lattice= new cjloss(dimension,0.94,i);	
		//cout << "lattice: " << endl << lattice->basis << endl; 

		BKZ_QP1(lattice->basis, 0.99, beta);		

		begin= clock();
		ComputeGS(lattice->basis,mu1,c1);
		end= clock();
		t_reduc+= (double)(end - begin) / CLOCKS_PER_SEC;

		for(int i= 0; i < dimension; i++)
			for(int j= 0; j < dimension; j++)
				conv(mu[i][j], mu1[i][j]);

		for(int i= 0; i < dimension; i++)
			conv(c[i], c1[i]);

		enumerate_epr(mu, c, Rvec, dimension, result, termination, enu_time);

		t_node+= enu_time;

		//lattice->print_solution(result);
		//cout << "enu_time= "  << enu_time << endl; 

		delete lattice;
	}

	t_node/= samples;
	t_reduc/= samples;
}

