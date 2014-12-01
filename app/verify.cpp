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

extern double ball_vol(int k, double r);

extern void profile_enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result); 

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim);

int main(int argc, char** argv) {
	int dim;
	mat_ZZ basis;
	vec_ZZ boundary;
	char* act_arg;

	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program verifies the node prediction." << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-l n\t\tThe name of the file containing the preprocessed lattice basis. (required)" << endl
 			<< "\t-b n\t\tThe name of the file containing the boundary function. (required)" << endl;
		return 0;
		}
 
	act_arg= get_cmd_option(argv, argv + argc, "-b");	
	if (act_arg) {
		ifstream basis_file(act_arg);
		if (basis_file.is_open()) {
			basis_file >> boundary;	
			cout << "# Bounding function: " << endl << "# " << boundary << endl;
		} else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-l");	
	if (act_arg) {
		ifstream basis_file(act_arg);
		if (basis_file.is_open()){
			basis_file >> basis;	
		} else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}

	dim= basis.NumRows();	

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
	for(int i= 0; i < mu1.NumRows(); i++) 
		gh*= sqrt(c[i]);
	gh= 1.05* pow(gh/ball_vol(dim, 1),1.0/dim);
	cout << "# Darmstadt bound: " << gh << endl;

// Enumeration
 	double* act= new double[dim];
	for(int i= 0; i<dim;i++)
		conv(act[i], boundary[i]);
	
	vec_RR solution;	
	profile_enumerate_epr(mu, c, act, dim, solution); 

// Prediction	
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

	RR* act_RR= new RR[dim];
	for(int i= 0; i < dim; i++) {
		conv(act_RR[i], boundary[i]);
		act_RR[i]= sqrt(act_RR[i]);
		}

	init_factorials(2*dim+2);
	predict_nodes_RR(act_RR, c, dim);

}
 
