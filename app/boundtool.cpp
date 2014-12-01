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
extern RR ball_vol_RR(int k, RR r);
extern void init_factorials(int up_to);
extern RR RR_PI;

// Defined in tools.cpp
char* get_cmd_option(char** begin, char** end, const string& option); 
bool cmd_option_exists(char** begin, char** end, const string& option);


int main(int argc, char** argv) {
	double t_node= 0; 
	double t_reduc= 0;
	unsigned long iterations= 1000;
	mat_ZZ basis;
	int dim;
	double* boundary;	
	ZZ v;
	v= 0;
	ZZ delta;
	delta= 1;

	stringstream ss;
	char* act_arg;
	
	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program computes boundary functions. The parameters t_node, t_reduc (these are the running time of the enumeration and the prerocessing reduction algorithms on the current platform) and a preprocessed basis of the lattice are required  Program options:" << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-f n\t\tThe name of the file containing the preprocessed lattice basis. (required)" << endl
 			<< "\t-n n\t\tThe (avarage) time the enumeration algorithm takes to process a single node. (required)" << endl
			<< "\t-r n\t\tThe (avarege) running time of the preprocessing reduction algorithm. (required)" << endl
			<< "\t-c n\t\tNumber of random changes to make during the numerical optimization. (default: 1000)" << endl
			<< "\t-d n\t\tThe size of the single random changes to meke during the optimization. (default: 1)" << endl
 			<< "\t-l n\t\tOptimize for shortest vector with length square root of n. (the Gaussian heuristic is default)" << endl
 			<< "\t-o n\t\tThe name of the file containing the information about previous computation that has to be continued. Just the -c switch can be used in this case." << endl;
		return 0;
	}

	if (!cmd_option_exists(argv, argv+argc, "-f") && !cmd_option_exists(argv, argv+argc, "-o")) {
		cout << "Input basis information is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-n") && !cmd_option_exists(argv, argv+argc, "-o")) {
		cout << "Parameter t_node is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-r") && !cmd_option_exists(argv, argv+argc, "-o")) {
		cout << "Parameter t_reduc is missing. Run again with -h for help." << endl;
		return 0;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-n");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-o")){
			cerr << "ERROR: invalid switch. In the presence of the -o switch the only other switch allowed is the -c" << endl;
			return 1;
			}
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
		if(cmd_option_exists(argv, argv+argc, "-o")){
			cerr << "ERROR: invalid switch. In the presence of the -o switch the only other switch allowed is the -c" << endl;
			return 1;
			}
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
		if(cmd_option_exists(argv, argv+argc, "-o")){
			cerr << "ERROR: invalid switch. In the presence of the -o switch the only other switch allowed is the -c" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> delta;
		ss.clear();
		if(delta < 1) {
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
		if(cmd_option_exists(argv, argv+argc, "-o")){
			cerr << "ERROR: invalid switch. In the presence of the -o switch the only other switch allowed is the -c" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> v;
		ss.clear();
		if(v <= 0) {
			cerr << "ERROR: invalid shortest vector length. The shortest vector length should be greater than zero. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-o");	
	if (act_arg) {
		ifstream infile(act_arg);
		if (infile.is_open()) { 
			continue_boundary_gen(infile, iterations, &boundary, dim); 
		} else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-f");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-o")){
			cerr << "ERROR: invalid switch. In the presence of the -o switch the only other switch allowed is the -c" << endl;
			return 1;
			}
		ifstream basis_file(act_arg);
		if (basis_file.is_open())
			basis_file >> basis;	
		else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}


	cout << "# Generated with cleanbkz " << CBKZ_VERSION << endl 
	<< "# Copyright (C) 2014 Janos Follath" << endl 
	<< "# This is free software with ABSOLUTELY NO WARRANTY." << endl << "#" << endl; 

	if(!cmd_option_exists(argv, argv+argc, "-o")){
		mat_RR mu1;
		vec_RR c1;
		ComputeGS(basis,mu1,c1);

		//lengths of the GS basis vectors
		RR* c= new RR[mu1.NumRows()];
		for(int i= 0; i < mu1.NumRows(); i++) {
			c[i].SetPrecision(RR_PRECISION);
			c[i]= SqrRoot(c1[i]);
		}

		dim= mu1.NumRows();
		boundary= new double[dim];	
	
		//length of the shortest vector in the cjloss lattice
		if(v>0) {
			//v= sqrt(v);
			cout << "# Using supplied lambda square: " << v << endl;
			}
		else {
			RR_PI.SetPrecision(RR_PRECISION);
			RR_PI= ComputePi_RR();	
			init_factorials(2*dim+1);
			RR one,gh,exp;
			one.SetPrecision(RR_PRECISION);
			gh.SetPrecision(RR_PRECISION);
			exp.SetPrecision(RR_PRECISION);
			gh= one= 1;	
			for(int i= 0; i < mu1.NumRows(); i++) {
				gh*= sqrt(c1[i]);
			}
			exp= 1.0/dim; 	
			pow(gh, gh/ball_vol_RR(dim, one), exp);
			gh*= 1.05;
			conv(v, gh*gh);
			cout << "# No lambda square supplied, using Gaussian heuristic: " << gh << endl;
		}
		

		double p_succ;
		double t_enum;	
		//NOTE: gs vector lengths and wanted vector length given
	
		generate_boundary(c, t_node, t_reduc, dim, boundary, v, delta, iterations, p_succ, t_enum, false); 

		cout << "# basis: '" << act_arg << "' " << endl
		<< "# estimated enumeration time: " << t_enum << endl  
		<< "# success probability: " << p_succ << endl; 
	}

	cout << "# boudary function: " << endl;
	

	vec_RR out;
	out.SetLength(dim);
	for(int i= 0; i < dim; i++)
		out[i]= boundary[i];		
	cout << "# " << out << endl << endl;

	for(int i= 0; i < dim; i++) {
		out[i].SetOutputPrecision(RR_PRECISION);
		cout << i << " " << out[i] << endl;
	}
	cout << endl;

	return 0;
}
