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
extern RR t_extreme_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n);
extern RR p_succ(RR Rvec[], int n);

// Defined in tools.cpp
char* get_cmd_option(char** begin, char** end, const string& option); 
bool cmd_option_exists(char** begin, char** end, const string& option);


int main(int argc, char** argv) {
	mat_ZZ basis;
	int dim;
	double* boundary;	
	ZZ R2;
	R2= 0;
	double t_node, t_reduc;
	t_node= t_reduc= 0;

	stringstream ss;
	char* act_arg;
	
	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program generates boundary functions for full enumeration, linear pruning and Schneider's polynomial pruning. A preprocessed basis of the lattice is required. Options -f -d and -s are mutually exclusive.  Program options:" << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-l n\t\tThe name of the file containing the preprocessed lattice basis. (required)" << endl
 			<< "\t-f n\t\tGenerate boundary without pruning" << endl
 			<< "\t-d n\t\tGenerate double step linear boundary" << endl
 			<< "\t-s n\t\tGenerate Schneiders boundary function (default)" << endl
 			<< "\t-n n\t\tThe (avarage) time the enumeration algorithm takes to process a single node. (optional)" << endl
			<< "\t-r n\t\tThe (avarege) running time of the preprocessing reduction algorithm. (optional)" << endl;
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

	if (!cmd_option_exists(argv, argv+argc, "-l")) {
		cout << "Input basis information is missing. Run again with -h for help." << endl;
		return 0;
		}

	if ((cmd_option_exists(argv, argv+argc, "-s") && cmd_option_exists(argv, argv+argc, "-d")) || (cmd_option_exists(argv, argv+argc, "-s") && cmd_option_exists(argv, argv+argc, "-f")) || (cmd_option_exists(argv, argv+argc, "-f") && cmd_option_exists(argv, argv+argc, "-d")) ) {
		cout << "Options -f -d and -s are mutually exclusive. Run again with -h for help." << endl;
		return 0;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-l");	
	if (act_arg) {
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
	//dim=110;
	boundary= new double[dim];	
	
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
	conv(R2, gh*gh);
	cout << "# Using Darmstadt bound (Gaussian heuristic times 1.05): " << gh << endl;
	cout << "# basis: '" << act_arg << "' " << endl;

	if(cmd_option_exists(argv, argv+argc, "-f")) {
		ZZ* act= new ZZ[dim];
		act[0]= act[1]= (2*R2)/dim;
		for(int i= 2; i < dim; i+=2) {
			act[i]= act[i+1]= act[i-1]+act[0];
			//cout << act[i] << endl;
		}
		act[dim-1]= R2;
		if(dim%2==0)
			act[dim-2]= R2;

		for(int i= 0; i < dim; i++) 
			conv(boundary[i],act[i]);
	} else if (cmd_option_exists(argv, argv+argc, "-d")) {
		for(int i= 0; i < dim; i++) 
			conv(boundary[i],R2);
	} else {
		/*double coeffs[]= {9.1e-4, 4e-2, -4e-3, 2.3e-4, -6.9e-6, 1.21e-7, -1.2e-9, 6.2e-12, -1.29e-14};
		for(int i= 1; i<=dim; i++) {
			boundary[i-1]= 0;
			for(int j= 0; j<=8; j++)	
				//boundary[i-1]+= pow(i*110/dim,j)*coeffs[j];
				boundary[i-1]+= pow(1.0*i/dim,j)*coeffs[j];
			//boundary[i-1]*=110/dim;*/
		double div = (double) dim;
		double y,scale;
		conv(scale,R2);
		for (double x = 1; x<=dim; x+=1){
			y = x/div;
			boundary[int(x)] = min(  -276.919520419850*pow(y,8) \
							 +1208.33321732860 *pow(y,7) \
							 -2128.79627565000 *pow(y,6) \
							 +1952.22801180000 *pow(y,5) \
							 -1012.11476080000 *pow(y,4) \
							 + 305.211610000000*pow(y,3) \
							 -  51.347076000000*pow(y,2) \
							 +   4.408932000000*pow(y,1) \
							 +   0.000914465000, 1.0);
			boundary[int(x)]*= scale;
		}

		for(int i= 0; i< dim; i+=2)
			boundary[i]= boundary[i+1]= floor(boundary[i+1]);
		
		for(int i= 1; i< dim-1; i+=2)
			if(boundary[i]>boundary[i+1])
				for(int j= i; j >= 0; j--) 
					if(boundary[j] > boundary[i+1])
						boundary[j]= boundary[i+1];
				
		
			
	}

	RR* act= new RR[dim];
	for(int i= 0; i< dim; i++)
		conv(act[i],boundary[i]);

	RR p_succ_r;
	RR t_enum;	
	if((t_node!=0)&&(t_reduc!=0))
		t_enum= t_extreme_RR(act, c, t_node, t_reduc, dim);


	RR* mod_f= new RR[dim];
	for(int j= 0; j < dim; j++) {
		conv(mod_f[j], act[j]); 
		mod_f[j]= sqrt(mod_f[j]);
	}
	p_succ_r= p_succ(mod_f, dim);

	cout << "# estimated enumeration time: " << t_enum << endl;
	cout << "# success probability: " << p_succ_r << endl; 
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
