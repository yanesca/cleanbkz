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

extern RR p_succ(RR Rvec[], int n);
extern void init_factorials(int);
extern RR RR_PI; 
extern void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, unsigned long &termination, double &time);

extern void predict_nodes_RR(RR  Rvec[], double b_star_norm[], int n); 

extern void predict_nodes(double Rvec[], double b_star_norm[], int n); 
extern double ball_vol(int k, double r);
extern void profile_enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result); 

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim);

// The maximum of the coefficients of the random unimodular matrix
#define UNIMOD_SIZE 100 
mat_ZZ unimod(int dim){
	mat_ZZ L, U;
	L.SetDims(dim, dim);
	U.SetDims(dim, dim);

	for(int i= 0; i < dim; i++)
		L[i][i]= U[i][i]= 1;

	for(int i= 0; i < dim; i++)
		for(int j= i+1; j < dim; j++) {
			U[i][j]= RandomBnd(UNIMOD_SIZE);
			L[j][i]= RandomBnd(UNIMOD_SIZE);
			}

	return L*U;
}


int main(int argc, char** argv) {
	int dim, bsize= 2;
	stringstream ss;
	mat_ZZ basis;
	int lattices, bases;
	int* seeds;
	int* lsucc;
	int* psucc;

	// dimension
	ss << argv[1];
	ss >> dim;
	ss.clear();

	// number of random lattices  
	ss << argv[2];
	ss >> lattices;
	ss.clear();
	seeds= new int[lattices];
	lsucc= new int[lattices];
	psucc= new int[lattices];

	// number of bases per lattice 
	ss << argv[3];
	ss >> bases;
	ss.clear();

	// reduction blocksize
	if(argc>4) {
		ss << argv[4];
		ss >> bsize;
		ss.clear();
	}

	double* polybound= new double[dim];
	double* linbound= new double[dim];

	// Schneider's polynomial approximation boundary function 
	double div = (double) dim;
	double y,scale;
	scale= dim;
	for (double x = 1; x<=dim; x+=1){
		y = x/div;
		polybound[int(x)] = min(  -276.919520419850*pow(y,8) \
						 +1208.33321732860 *pow(y,7) \
						 -2128.79627565000 *pow(y,6) \
						 +1952.22801180000 *pow(y,5) \
						 -1012.11476080000 *pow(y,4) \
						 + 305.211610000000*pow(y,3) \
						 -  51.347076000000*pow(y,2) \
						 +   4.408932000000*pow(y,1) \
						 +   0.000914465000, 1.0);
		polybound[int(x)]*= scale;
	}

	for(int i= 0; i< dim; i+=2)
		polybound[i]= polybound[i+1]= floor(polybound[i+1]);
	
	for(int i= 1; i< dim-1; i+=2)
		if(polybound[i]>polybound[i+1])
			for(int j= i; j >= 0; j--) 
				if(polybound[j] > polybound[i+1])
					polybound[j]= polybound[i+1];
			
	// linear boudary function
	ZZ* act= new ZZ[dim];
	ZZ R2;
	conv(R2, dim);
	act[0]= act[1]= (2*R2)/dim;
	for(int i= 2; i < dim; i+=2) {
		act[i]= act[i+1]= act[i-1]+act[0];
	}
	act[dim-1]= R2;
	if(dim%2==0)
		act[dim-2]= R2;

	for(int i= 0; i < dim; i++) 
		conv(linbound[i],act[i]);


	// Computing success probabilities
	RR_PI.SetPrecision(RR_PRECISION);
	RR_PI= ComputePi_RR();	
	init_factorials(2*dim+2);
	RR* mod_f= new RR[dim];
	cout << "Linear boundary: " << endl;
	for(int j= 0; j < dim; j++) {
		cout << linbound[j] << " ";
		conv(mod_f[j], linbound[j]); 
		mod_f[j]= sqrt(mod_f[j]);
	}
	cout << endl << "Linear success rate: " << p_succ(mod_f, dim) << endl;

	cout << "Polynomial boundary: " << endl;
	for(int j= 0; j < dim; j++) {
		cout << polybound[j] << " ";
		conv(mod_f[j], polybound[j]); 
		mod_f[j]= sqrt(mod_f[j]);
	}
	cout << endl << "Polynomial success rate: " << p_succ(mod_f, dim) << endl;


	// Performing randomization followed by enumeration
	unsigned long termination= 0;
	mat_RR mu1;
	vec_RR c1;
	vec_RR solution;	
	double* c= new double[dim];
	double** mu= new double*[dim];
	double t_enum;

	for(int i= 0; i<lattices; i++) {

		seeds[i]= rand();
		lsucc[i]= psucc[i]= 0;
		cjloss l(dim, 0.94, seeds[i]);
		basis= l.get_basis(bsize);

		cout << "\nLattice " << seeds[i] << "\t" << flush;

		ComputeGS(basis,mu1,c1);
		for(int j= 0; j < mu1.NumRows(); j++) 
			conv(c[j], c1[j]);
	
		for(int j= 0; j < mu1.NumRows(); j++)
			mu[j]= new double[mu1.NumCols()]; 
		for(int j= 0; j < mu1.NumRows(); j++)
			for(int k= 0; k < mu1.NumCols(); k++)
				conv(mu[j][k], mu1[j][k]);

		enumerate_epr(mu, c, linbound, basis.NumRows(), solution, termination, t_enum);
		if(solution.length()!=0)
			lsucc[i]++;

		enumerate_epr(mu, c, polybound, basis.NumRows(), solution, termination, t_enum);
		if(solution.length()!=0)
			psucc[i]++;

		cout << "*" << flush;

		for(int l= 1; l<bases; l++) {
			basis= basis*unimod(basis.NumRows());
			LLL_XD(basis, 0.99);
			BKZ_FP(basis, 0.99, bsize); 
			
			ComputeGS(basis,mu1,c1);
			for(int j= 0; j < mu1.NumRows(); j++) 
				conv(c[j], c1[j]);
		
			for(int j= 0; j < mu1.NumRows(); j++)
				mu[j]= new double[mu1.NumCols()]; 
			for(int j= 0; j < mu1.NumRows(); j++)
				for(int k= 0; k < mu1.NumCols(); k++)
					conv(mu[j][k], mu1[j][k]);

			enumerate_epr(mu, c, linbound, basis.NumRows(), solution, termination, t_enum);
			if(solution.length()!=0)
				lsucc[i]++;

			enumerate_epr(mu, c, polybound, basis.NumRows(), solution, termination, t_enum);
			if(solution.length()!=0)
				psucc[i]++;

			cout << "*" << flush;

		}

	}

	double l, p;
	l= p= 0;
	cout << endl << endl << "Results: "<< endl << "\t  seed\tlinear\tpolynomial" << endl;
	for(int i= 0; i<lattices; i++){
		cout << "    " << seeds[i] << "\t  " << lsucc[i] << "\t  " << psucc[i] << endl;	
		l+=  lsucc[i];
		p+= psucc[i];
		}

	cout << endl << endl << "Overall: " << endl;
	cout << "Linear: " << l/(lattices*bases) << endl;
	cout << "Polynomial: " << p/(lattices*bases) << endl;
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
