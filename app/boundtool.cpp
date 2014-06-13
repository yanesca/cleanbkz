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
#include <sstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/version.hpp>
#include <cleanbkz/tools.hpp>
#include <NTL/LLL.h>

using namespace std;

// Defined in boundary.cpp
extern double ball_vol(int k, double r);

// Defined in tools.cpp
char* get_cmd_option(char** begin, char** end, const string& option); 
bool cmd_option_exists(char** begin, char** end, const string& option);


int main(int argc, char** argv) {
	double v= 0;
	double t_node= 0; 
	double t_reduc= 0;
	double delta= 1e-1;
	unsigned long iterations= 1000;
	mat_ZZ basis;

	stringstream ss;
	char* act_arg;
	
	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program computes boundary functions. The parameters t_node, t_reduc (these are the running time of the enumeration and the prerocessing reduction algorithms on the current platform) and a preprocessed basis of the lattice are required  Program options:" << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-f n\t\tThe name of the file containing the preprocessed lattice basis. (required)" << endl
 			<< "\t-n n\t\tThe (avarage) time the enumeration algorithm takes to process a single node. (required)" << endl
			<< "\t-r n\t\tThe (avarege) running time of the preprocessing reduction algorithm. (required)" << endl
			<< "\t-c n\t\tNumber of random changes to make during the numerical optimization. (default: 1000)" << endl
			<< "\t-d n\t\tThe size of the single random changes to meke during the optimization. (default: 0.1)" << endl
 			<< "\t-l n\t\tOptimize for shortest vector with length square root of n. (the Gaussian heuristic is default)" << endl;
		return 0;
	}

	if (!cmd_option_exists(argv, argv+argc, "-f")) {
		cout << "Input basis is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-n")) {
		cout << "Parameter t_node is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-r")) {
		cout << "Parameter t_reduc is missing. Run again with -h for help." << endl;
		return 0;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-n");	
	if (act_arg) {
		ss << act_arg;
		ss >> t_node;
		ss.clear();
		if(t_node < 0) {
			cerr << "ERROR: invalid t_node. The running time should be greater than or equal to zero. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-r");	
	if (act_arg) {
		ss << act_arg;
		ss >> t_reduc;
		ss.clear();
		if(t_reduc < 0) {
			cerr << "ERROR: invalid t_reduc. The running time should be greater than or equal to zero. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-d");	
	if (act_arg) {
		ss << act_arg;
		ss >> delta;
		ss.clear();
		if(delta <= 0) {
			cerr << "ERROR: invalid delta. The optimization step should be greater than zero. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-c");	
	if (act_arg) {
		ss << act_arg;
		ss >> iterations;
		ss.clear();
		if(iterations <= 0) {
			cerr << "ERROR: invalid number of changes. The number of optimization steps should be greater than zero. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-l");	
	if (act_arg) {
		ss << act_arg;
		ss >> v;
		ss.clear();
		if(v <= 0) {
			cerr << "ERROR: invalid shortest vector length. The shortest vector length should be greater than zero. Aborting." << endl;
			return 1;
			}
		}


	act_arg= get_cmd_option(argv, argv + argc, "-f");	
	if (act_arg) {
		ifstream basis_file(act_arg);
		if (basis_file.is_open())
			basis_file >> basis;	
		else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(basis,mu1,c1);

	//lengths of the GS basis vectors
	RR* c= new RR[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) {
		c[i]= SqrRoot(c1[i]);
		c[i].SetPrecision(RR_PRECISION);
	}

	int dim= mu1.NumRows();
	double* boundary= new double[dim];	
	
	cout << "# Generated with cleanbkz " << CBKZ_VERSION << endl 
	<< "# Copyright (C) 2014 Janos Follath" << endl 
	<< "# This is free software with ABSOLUTELY NO WARRANTY." << endl << "#" << endl; 

	//length of the shortest vector in the cjloss lattice
	if(v>0) {
		v= sqrt(v);
		cout << "# Using supplied lambda : " << v << endl;
		}
	else {
		double tmp,gh= 1;
		double* gsghs= new double[dim]; 
		for(int i= 0; i < mu1.NumRows(); i++) {
			conv(tmp, c1[i]);
			gh*= tmp;
			gsghs[i]= pow(gh/ball_vol(i+1, 1),1.0/(i+1));
		}
		delete gsghs;
		gh= pow(gh/ball_vol(dim, 1),1.0/dim);
		cout << "# No lambda supplied, using Gaussian heuristic: " << gh << endl;
	}
		


	// TODO: csinálni egy változatot, ahol nem iterationt hanem thressholdot adunk meg
	double p_succ;
	double t_enum;	
	generate_boundary(c, t_node, t_reduc, dim, boundary, v, delta, iterations, p_succ, t_enum, false); 

	cout << "# basis: '" << act_arg << "' " << endl
	<< "# estimated enumeration time: " << t_enum << endl  
	<< "# success probability: " << p_succ << endl  
	<< "# boudary function: " << endl;

	vec_RR out;
	out.SetLength(dim);
	for(int i= 0; i < dim; i++)
		out[i]= boundary[i];		
	cout << "# " << out << endl << endl;

	for(int i= 0; i < dim; i++)
		cout << i << " " << boundary[i] << endl;
	cout << endl;

	return 0;
}
