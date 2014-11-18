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

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#include <iostream>
#include <sstream>
#include <fstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/enumeration.hpp>
#include <cleanbkz/tools.hpp>
#include <NTL/LLL.h>

using namespace std;

extern void init_factorials(int);

extern void predict_nodes_RR(RR  Rvec[], double b_star_norm[], int n); 

extern void predict_nodes(double Rvec[], double b_star_norm[], int n); 
extern double ball_vol(int k, double r);
extern void profile_enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result); 

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim);

int main(int argc, char** argv) {
	int dim, bsize= 2;
	stringstream ss;
	bool l_cjloss= true;
	mat_ZZ basis;

	ss << argv[1];
	ss >> dim;
	ss.clear();

	init_factorials(2*dim+2);

	if (l_cjloss) {
		//cjloss l(dim, 0.94, time(NULL));
		cjloss l(dim, 0.94, 678);
		basis= l.get_basis(bsize);
	} else {
		ZZ determinant;
		GenPrime(determinant,dim);
		gen_randlat(basis,determinant,dim); 
		BKZ_QP1(basis, 0.99, bsize); 
	}

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], c1[i]);

	
	double** mu= new double*[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		mu[i]= new double[mu1.NumCols()]; 
	for(int i= 0; i < mu1.NumRows(); i++)
		for(int j= 0; j < mu1.NumCols(); j++)
			conv(mu[i][j], mu1[i][j]);

// Gaussian heuristic
	double gh= 1;
	double* gsghs= new double[dim]; 
	for(int i= 0; i < mu1.NumRows(); i++) {
		gh*= c[i];
		gsghs[i]= pow(gh/ball_vol(i+1, 1),1.0/(i+1));
		}
	gh= pow(gh/ball_vol(dim, 1),1.0/dim);
	if(l_cjloss)
		cout << "# Shortest vector length: " << sqrt(dim) << endl;
	cout << "# Gaussian heuristic: " << gh << endl;

// Bounding functions
	double R;
	if(l_cjloss)	
		R= dim;
	else
		R= gh*gh;

	double* full= new double[dim];
	for(int i= 0; i < dim; i++) 
		full[i]= R;

	double* sh= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		sh[2*i]= sh[2*i+1]= MIN(1,R*sqrt(1.05*2*i/dim));
	if(dim%2==0)
		sh[dim-1]= R;	
	else
		sh[dim-2]= sh[dim-1]= R;	

	double* linear= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		linear[2*i]= linear[2*i+1]= 2*(i+1)*R/dim;
	if(dim%2==0)
		linear[dim-1]= R;	
	else
		linear[dim-2]= linear[dim-1]= R;	

	double* step= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		if(i<dim/4)
			step[2*i]= step[2*i+1]= R/2;
		else
			step[2*i]= step[2*i+1]= R;	
	if(dim%2==0)
		step[dim-1]= R;	
	else
		step[dim-2]= step[dim-1]= R;	

	double* min= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		min[2*i]= min[2*i+1]= pow(gsghs[2*i+1],2);

	double* hmin= new double[dim];
	hmin[0]= c[dim-1]*c[dim-1]*0.5;
	for(int i= 1; i < dim; i++) {
		hmin[i]= hmin[i-1]*(1+1.0/(i+1));	
		hmin[i]= MIN(hmin[i], R);	
		if(i%2==1)
			hmin[i-1]= hmin[i];
	}

	double* manual= new double[dim];
	for(int i= 0; i < dim; i++) 
		manual[i]= R;

	for(int i= 2; i< argc; i++) {
		ss << argv[i];
		ss >> manual[i-2];
		ss.clear();
		}
	
	double* act= manual;

// Enumeration
//	vec_RR solution;	
//	profile_enumerate_epr(mu, c, act, dim, solution); 

// Prediction	
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

	for(int i= 0; i < dim; i++) { 
		full[i]= sqrt(full[i]);
		step[i]= sqrt(step[i]);
		linear[i]= sqrt(linear[i]);
		min[i]= sqrt(min[i]);
		hmin[i]= sqrt(hmin[i]);
		manual[i]= sqrt(manual[i]);
		}

	//predict_nodes(act, c, dim);
	RR* act_RR= new RR[dim];
	for(int i= 0; i < dim; i++)
		act_RR[i]= act[i];
	predict_nodes_RR(act_RR, c, dim);


/* GSA	
	cout << endl << endl;
	double alpha=1;
	for(int i= 0; i < mu1.NumRows(); i++) 
		alpha*= c[i];
	alpha= log(pow(alpha,2.0/dim));
	n_full_gsa(gh, alpha, bsize, dim);*/ 

}

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim) {
	basis.SetDims(dim, dim);
	
	ZZ s;
	s= time(NULL);
	SetSeed(s);
	for(int i= 0; i < dim; i++) {
		basis[i][0]= RandomBnd(determinant);
		basis[i][i]= 1;
	}

	basis[0][0]= determinant;
} 
