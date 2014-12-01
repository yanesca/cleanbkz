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

#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <NTL/LLL.h>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/tools.hpp>
#include <cleanbkz/version.hpp>

using namespace std;
using namespace NTL;

//TODO: random lattices

void measure_epr(int dimension, int beta, unsigned long samples, unsigned long termination, 
		double& t_node, double& t_reduc, unsigned long& zeroes, unsigned long& bias);

void measure_epr(mat_ZZ basis, int beta, unsigned long termination, double& t_node, double& t_reduc);

// Defined in boundary.cpp
extern RR ball_vol_RR(int k, RR r);
extern void init_factorials(int up_to);
extern RR RR_PI;

int main(int argc, char** argv) {
  	int dimension= 0;
	int bkz_blocksize= 2;
	int nodes= 1000000;
	int samples= 100; 
	int start= 0; 
	int end= 0;
	stringstream ss;
	char* act_arg;
	clock_t begin, finish;

	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program uses cjloss lattices to measure the parameters t_node and t_reduc required for the computation of boundary functions. These are the running time of the enumeration and reduction algorithms on the current platform (in the case the cjloss lattices the dimension of the preprocessing is one higher than the enumeration). Program options:" << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-l filename\tPerform tests on the given lattice. (no other switch is allowed just the -b and the -k)" << endl
 			<< "\t-d n\t\tMeasure running times in dimension n." << endl
 			<< "\t-s n\t\tStart the measurements in dimension n. It is ignored when -d is given." << endl
 			<< "\t-e n\t\tContinue measuring in all dimensions until dimension n with a step five. It is ignored when -d is given." << endl
			<< "\t-b n\t\tAbort enumeration after processing n nodes. (default: 10,000,000)" << endl
			<< "\t-n n\t\tNumber of experiments to make. (default: 100)" << endl
			<< "\t-k n\t\tThe blocksize of BKZ used for preprocessing. (default: 2)" << endl
			<< "\t-t\t\tPredict the running time of the experiments." << endl;
		return 0;
	}

	//TODO: exception handling
	act_arg= get_cmd_option(argv, argv + argc, "-k");	
	if (act_arg) {
		ss << act_arg;
		ss >> bkz_blocksize;
		ss.clear();
		if(bkz_blocksize < 2) {
			cerr << "ERROR: Invalid blocksize. The BKZ blocksize should be greater than or equal to two. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-d");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-l")){
			cerr << "ERROR: invalid switch. In the presence of the -l switch the only other switch allowed is the -b and the -k" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> dimension;
		ss.clear();
		if(dimension < 30) {
			cerr << "ERROR: Dimension too low for mesurements. Please measure higher dimensions and extrapolate the results." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-n");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-l")){
			cerr << "ERROR: invalid switch. In the presence of the -l switch the only other switch allowed is the -b and the -k" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> samples;
		ss.clear();
		if(samples < 1) {
			cerr << "ERROR: Number of samples less than one. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-b");	
	if (act_arg) {
		ss << act_arg;
		ss >> nodes;
		ss.clear();
		if(nodes < 1) {
			cerr << "ERROR: Bound on processed nodes is less than one. Aborting." << endl;
			return 1;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-s");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-l")){
			cerr << "ERROR: invalid switch. In the presence of the -l switch the only other switch allowed is the -b and the -k" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> start;
		ss.clear();
		if(start < 45) {
			cerr << "ERROR: Start dimension too low for mesurements, setting start dimension to 45." << endl;
			start= 45;
			}
		}

	act_arg= get_cmd_option(argv, argv + argc, "-e");	
	if (act_arg) {
		if(cmd_option_exists(argv, argv+argc, "-l")){
			cerr << "ERROR: invalid switch. In the presence of the -l switch the only other switch allowed is the -b and the -k" << endl;
			return 1;
			}
		ss << act_arg;
		ss >> end;
		ss.clear();
		if(end < start) {
			cerr << "ERROR: End dimension lower than start dimension. Aborting" << endl;
			return 1;
			}
		}

	mat_ZZ basis;
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

	cout << "cleanbkz " << CBKZ_VERSION << endl 
	<< "Copyright (C) 2014 Janos Follath" << endl 
	<< "This is free software with ABSOLUTELY NO WARRANTY." << endl << endl; 

	unsigned long bias, zeroes;
	double t_node, t_reduc;
	if (cmd_option_exists(argv, argv+argc, "-d")) {
		if (cmd_option_exists(argv, argv+argc, "-t")) {
			cout << "Predicting experiment time for dimension " << dimension << ": "; 
			begin= clock();
			measure_epr(dimension, bkz_blocksize, 10, nodes, t_node, t_reduc, zeroes, bias); 
			finish= clock();
			cout << ((double)(finish - begin) / CLOCKS_PER_SEC)/10*samples << "s" << endl;
		}

		cout << "Measuring dimension " << dimension << ":" << endl;
		measure_epr(dimension, bkz_blocksize, samples, nodes, t_node, t_reduc, zeroes, bias); 
		cout << "t_node= " << t_node << endl;
		cout << "t_reduc= " << t_reduc << endl;
		cout << "biased/zeroes: " << bias << "/" << zeroes << endl << endl;
	} else if(cmd_option_exists(argv, argv+argc, "-l")){
		cout << "Measuring timing for: " << get_cmd_option(argv, argv + argc, "-l") << endl;

		measure_epr(basis, bkz_blocksize, nodes, t_node, t_reduc);

		cout << "t_node= " << t_node << endl;
		cout << "t_reduc= " << t_reduc << endl;
	} else {
		if (cmd_option_exists(argv, argv+argc, "-t")) {
			unsigned long int total= 0;
			for(int i= start; i <= end; i+=5) {
				cout << "Predicting experiment time for dimension " << i << ": "; 
				begin= clock();
				measure_epr(i, bkz_blocksize, 10, nodes, t_node, t_reduc, zeroes, bias); 
				finish= clock();
				cout << ((double)(finish - begin) / CLOCKS_PER_SEC)/10*samples << "s" << endl;
				total+= (unsigned long int) (((double)(finish - begin) / CLOCKS_PER_SEC)/10*samples);
				}
			cout << "Expected total running time of the experiments: " << total/3600 << "h" 
					<< (total%3600)/60 << "m" << (total%36000)%60 << "s" << endl << endl;
			}

		for(int i= start; i <= end; i+=5) {
			cout << "Measuring dimension " << i << ":" << endl;
			measure_epr(i, bkz_blocksize, samples, nodes, t_node, t_reduc, zeroes, bias); 
			cout << "t_node= " << t_node << endl;
			cout << "t_reduc= " << t_reduc << endl;
			cout << "biased/zeroes: " << bias << "/" << zeroes << endl << endl;
			}

		cout << "(Sometimes the enumeration doesn't have to inspect the prescribed number of nodes or even nodes enough " 
			<< "to make the runing time measureable, the number of these cases are in the \"biased/zeroes\" line." 
			<< "The first number measures only the biased cases with nonzero running time.)" << endl << endl;
	}

	return 0;
}

extern void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, 
		unsigned long &termination, double &time); 

void measure_epr(int dimension, int beta, unsigned long samples, unsigned long termination, 
		double& t_node, double& t_reduc, unsigned long& zeroes, unsigned long& bias) {
	double enu_time= 0;
	cjloss* lattice;
	mat_ZZ tmp,basis;
	mat_RR mu1;
	vec_RR c1;
	clock_t begin, end;
	vec_RR result;
	unsigned long nodes;

	t_node= t_reduc= 0;
	zeroes= bias= 0;

	double** mu= new double*[dimension];
	for(int i= 0; i < dimension; i++)
		mu[i]= new double[dimension]; 

	double* c= new double[dimension];

	double* Rvec= new double[dimension];
	for(int i= 0; i < dimension; i++)
		Rvec[i]= dimension;

	for(unsigned long i= 0; i < samples; i++) {	
		lattice= new cjloss(dimension,0.94,i);	

		tmp= lattice->get_basis(0);		
		begin= clock();
		LLL_XD(tmp, 0.99);		
		BKZ_FP(tmp, 0.99, beta);		
		end= clock();
		t_reduc+= (double)(end - begin) / CLOCKS_PER_SEC;

		basis.SetDims(tmp.NumRows()-1, tmp.NumRows()-1);
		for(int i= 0; i<basis.NumRows(); i++)
			for(int j= 0; j<basis.NumCols(); j++)
				basis[i][j]= tmp[i][j];

		ComputeGS(basis,mu1,c1);

		for(int i= 0; i < dimension; i++)
			for(int j= 0; j < dimension; j++)
				conv(mu[i][j], mu1[i][j]);

		for(int i= 0; i < dimension; i++)
			conv(c[i], c1[i]);

		nodes= termination;
		enumerate_epr(mu, c, Rvec, dimension, result, nodes, enu_time);

		if(enu_time < 1.0/CLOCKS_PER_SEC/termination) 
			zeroes++;
		else if (nodes!=0)
			bias++;
		t_node+= enu_time;

		delete lattice;
	}

	t_node/= samples;
	t_reduc/= samples;
}

void measure_epr(mat_ZZ basis, int beta, unsigned long termination, double& t_node, double& t_reduc) {
	double enu_time= 0;
	mat_RR mu1;
	vec_RR c1;
	clock_t begin, end;
	vec_RR result;
	unsigned long nodes;

	t_node= t_reduc= 0;

	int dimension= basis.NumRows();

	double** mu= new double*[dimension];
	for(int i= 0; i < dimension; i++)
		mu[i]= new double[dimension]; 

	double* c= new double[dimension];

	begin= clock();
	LLL_XD(basis, 0.99);		
	BKZ_FP(basis, 0.99, beta);		
	end= clock();

	ComputeGS(basis,mu1,c1);

	for(int i= 0; i < dimension; i++)
		for(int j= 0; j < dimension; j++)
			conv(mu[i][j], mu1[i][j]);

	for(int i= 0; i < dimension; i++)
		conv(c[i], c1[i]);


	RR_PI.SetPrecision(RR_PRECISION);
	RR_PI= ComputePi_RR();	
	init_factorials(2*dimension+1);
	RR one,gh,exp;
	one.SetPrecision(RR_PRECISION);
	gh.SetPrecision(RR_PRECISION);
	exp.SetPrecision(RR_PRECISION);
	gh= one= 1;	
	for(int i= 0; i < mu1.NumRows(); i++) {
		gh*= sqrt(c1[i]);
		}
	exp= 1.0/dimension; 	
	pow(gh, gh/ball_vol_RR(dimension, one), exp);
	gh*= 1.05;

	cout << "Using Darmstadt bound: " << gh << endl;

	double* Rvec= new double[dimension];
	for(int i= 0; i < dimension; i++)
		conv(Rvec[i], gh*gh);

	t_reduc= (double)(end - begin) / CLOCKS_PER_SEC;
	nodes= termination;
	enumerate_epr(mu, c, Rvec, dimension, result, nodes, enu_time);

	if(enu_time < 1.0/CLOCKS_PER_SEC/termination) 
		cerr << "Warning! Zero enumeration time!" << endl;
	else if (nodes!=0)
		cerr << "Warning! Nodes traversed: " << (termination-nodes) << endl;
	t_node= enu_time;

}

