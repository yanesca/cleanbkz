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
#include <fstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/enumeration.hpp>
#include <cleanbkz/tools.hpp>
#include <cleanbkz/cjloss.hpp>
#include <NTL/LLL.h>

using namespace std;

extern void predict_nodes_RR(RR  Rvec[], double b_star_norm[], int n); 
extern void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, unsigned long &termination, double &time);

extern void profile_enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result); 

mat_ZZ unimod(int dim){
	mat_ZZ L, U;
	L.SetDims(dim, dim);
	U.SetDims(dim, dim);

	for(int i= 0; i < dim; i++)
		L[i][i]= U[i][i]= 1;

	for(int i= 0; i < dim; i++)
		for(int j= i+1; j < dim; j++) {
			U[i][j]= RandomBnd(2);
			L[j][i]= RandomBnd(2);
			}

	return L*U;
}

int main(int argc, char** argv) {
	mat_ZZ basis;
	vec_RR bnd;
	char* act_arg;
	int beta= 0;

	
	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program performs lattice enumeration with extreme pruning:" << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-v \t\tPerforms a single enumeration round and compares the nodes to the prediction." << endl
 			<< "\t-l filename\tReads the lattice from the file filename. (This option is mandatory)" << endl
 			<< "\t-b filename\tReads the boundary function from the file filename. (This option is mandatory)" << endl;
		return 0;
	}

	act_arg= get_cmd_option(argv, argv + argc, "-l");	
	if (act_arg) {
		ifstream basis_file(act_arg);
		if (basis_file.is_open()) {
			basis_file >> basis;		
			basis_file >> beta;	
			if(beta<2) {
				cerr << "ERROR: the reduction block size given in the lattice file can't be less than 2. Aborting." << endl;
				return 1;
			}
		} else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 2;
			}
		}
	else {
		cerr << "ERROR: The option '-l' is mandatory. Aborting." << endl;
		return 1;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-b");	
	if (act_arg) {
		ifstream bnd_file(act_arg);
		if (bnd_file.is_open()) 
			bnd_file >> bnd;	
		else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 2;
			}
		}
	else {
		cerr << "ERROR: The option '-b' is mandatory. Aborting." << endl;
		return 1;
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

	double* boundary= new double[basis.NumRows()];
	for(int i= 0; i < basis.NumRows(); i++){
		conv(boundary[i], bnd[i]);	
		boundary[i]*= boundary[i];
		}

	vec_RR solution;	
	
	if (cmd_option_exists(argv, argv+argc, "-v")) {
		profile_enumerate_epr(mu, c, boundary, basis.NumRows(), solution); 

		for(int i= 0; i < mu1.NumRows(); i++) 
			conv(c[i], sqrt(c1[i]));

		predict_nodes_RR(bnd.elts(), c, basis.NumRows());

		return 0;
		}

	clock_t begin;	
	unsigned long termination= 0;
	double t_enum, t_reduc, t_all;
	t_all= t_enum= t_reduc= 0;
	enumerate_epr(mu, c, boundary, basis.NumRows(), solution, termination, t_enum);
	t_all+= t_enum;

	while(solution.length()==0) {
		cout << "Enumeration failed. Randomizing lattice basis..." << endl;

		begin= clock();	
		basis= basis*unimod(basis.NumRows());
		BKZ_QP1(basis, 0.99, beta); 

		ComputeGS(basis,mu1,c1);

		for(int i= 0; i < mu1.NumRows(); i++) 
			conv(c[i], c1[i]);
		
		for(int i= 0; i < mu1.NumRows(); i++)
			for(int j= 0; j < mu1.NumCols(); j++)
				conv(mu[i][j], mu1[i][j]);

		t_reduc+= clock()-begin;	

		cout << "Done." << endl;
		cout << "Performing enumeration..." << endl;

		begin= clock();	
		enumerate_epr(mu, c, boundary, basis.NumRows(), solution, termination, t_enum);
		t_all+= clock()-begin;
		cout << "Done." << endl;
		}

	cout << "Enumeration successful. Shortest vector coordinates in the basis: " << endl << solution << endl;
	
	double* sol= new double[basis.NumRows()];
	for(int i= 0; i < basis.NumRows(); i++)
		sol[i]= 0;

	double dbase, dsol;
	for(int i= 0; i< basis.NumRows(); i++)
		for(int j= 0; j< basis.NumCols(); j++) {	
			conv(dbase, basis[i][j]);
			conv(dsol, solution[i]);
			sol[j]+= dbase*dsol;
			}

	cout << "Coordinates in the lattice: " << endl << "[";
	for(int i= 0; i< basis.NumRows()-1; i++) 
		cout << sol[i] << " "; 	
	cout << sol[basis.NumRows()-1] << "]" << endl; 
	
	double sqrdLength= 0;
	for(int i= 0; i< basis.NumRows(); i++) 
		sqrdLength+= sol[i]*sol[i];

	cout << "Squared length of the shortest vector: " << sqrdLength << endl << endl; 

	cout << "Time spent with enumeration: " << t_all/CLOCKS_PER_SEC << endl;
	cout << "Time spent with reduction: " << t_reduc/CLOCKS_PER_SEC << endl;
	cout << "Total running time: " << (t_all + t_reduc)/CLOCKS_PER_SEC << endl;
	cout << "(Note that the basis was already reduced and thus the time of the first reduction is not included. Also the enumeration subrutine exits after finding the first acceptable vector and does not search through the whole pruned enumeration tree.)" << endl; 
}
 
