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
extern void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, 
		unsigned long &termination, double &time); 


void measure_epr(mat_ZZ& basis, int beta, unsigned long termination, double& t_node, double& t_reduc); 

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

int main(int argc, char** argv) {

	int maxsize= 0;
	int dim= 0;
	int beta= 0;
	int seed= 0;
	double density= 0;

	stringstream ss;
	char* act_arg;
	
	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program generates random lattices. The presence of -c or -m switch is required." << endl
 			<< "\t-d n\t\tThe dimension of the lattice. (required)" << endl
			<< "\t-c n\t\tGenerates a CJLOSS lattice of density of n." << endl
			<< "\t-m n\t\tGenerates a random lattice with entries of length at most n. It is ignored when -c is given." << endl
 			<< "\t-h \t\tPrint this help." << endl
 			<< "\t-s n\t\tUse n as the seed of the pseudorandom generator." << endl
 			<< "\t-p filename\tPreprocess an existing lattice basis" << endl
 			<< "\t-t \t\tMeasure and output enumeration timing information" << endl
			<< "\t-k n\t\tReduces the generated basis with blocksize of n." << endl;
		return 0;
	}

	if (!cmd_option_exists(argv, argv+argc, "-d") && !cmd_option_exists(argv, argv+argc, "-p")) {
		cout << "Lattice dimension is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-m") && !cmd_option_exists(argv, argv+argc, "-c") && !cmd_option_exists(argv, argv+argc, "-p")) {
		cout << "ERROR: neither maximum size of entries nor density is provided. Aborting." << endl;
		return 0;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-d");	
	if (act_arg) {
		ss << act_arg;
		ss >> dim;
		ss.clear();
		if(dim < 4) {
			cerr << "ERROR: invalid dimesion. The dimension should be at least four. Aborting." << endl;
			return 1;
			}
		}
	else if (cmd_option_exists(argv, argv+argc, "-d")) {	
		cerr << "ERROR: can't parse option '-d'. Aborting." << endl;
		return 1;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-m");	
	if (act_arg) {
		ss << act_arg;
		ss >> maxsize;
		ss.clear();
		if(maxsize <= 0) {
			cerr << "ERROR: invalid maximum size. The maximum size of entries should be greater than zero. Aborting." << endl;
			return 1;
			}
		}
	else if (cmd_option_exists(argv, argv+argc, "-m")) {	
		cerr << "ERROR: can't parse option '-m'" << endl;
		return 1;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-k");	
	if (act_arg) {
		ss << act_arg;
		ss >> beta;
		ss.clear();
		if(beta < 2) {
			cerr << "ERROR: invalid blocksize. The block size of the reduction should be greater than one. Aborting." << endl;
			return 1;
			}
		}
	else if (cmd_option_exists(argv, argv+argc, "-k")) {	
		cerr << "ERROR: can't parse option '-k'" << endl;
		return 1;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-c");	
	if (act_arg) {
		ss << act_arg;
		ss >> density;
		ss.clear();
		if(density <= 0) {
			cerr << "ERROR: invalid density. The cjloss lattice density should be greater than zero. Aborting." << endl;
			return 1;
			}
		}
	else if (cmd_option_exists(argv, argv+argc, "-c")) {	
		cerr << "ERROR: can't parse option '-c'" << endl;
		return 1;
		}

	act_arg= get_cmd_option(argv, argv + argc, "-s");	
	if (act_arg) {
		ss << act_arg;
		ss >> seed;
		ss.clear();
		}
	else if (cmd_option_exists(argv, argv+argc, "-s")) {	
		cerr << "ERROR: can't parse option '-s'" << endl;
		return 1;
		}

	cerr << "cleanbkz " << CBKZ_VERSION << endl 
	<< "Copyright (C) 2014 Janos Follath" << endl 
	<< "This is free software with ABSOLUTELY NO WARRANTY." << endl << endl; 

	mat_ZZ m;
	if (cmd_option_exists(argv, argv+argc, "-c")) {	
		cjloss l(dim, density, seed);
		m= l.get_basis(beta);
		} 
	else if (cmd_option_exists(argv, argv+argc, "-m")) {
		ZZ determinant;
		GenPrime(determinant,maxsize);
		gen_randlat(m,determinant,dim); 
		}
	else if (cmd_option_exists(argv, argv+argc, "-p")) {
		act_arg= get_cmd_option(argv, argv + argc, "-p");	
		ifstream basis_file(act_arg);
		if (basis_file.is_open())
			basis_file >> m;	
		else {
			cerr << "ERROR: can't open input file: '" << act_arg << "'. Aborting." << endl;
			return 1;
			}
		}
	
	double t_node, t_reduc;
	if(beta > 1 && !(cmd_option_exists(argv, argv+argc, "-c"))) {	

		if(cmd_option_exists(argv, argv+argc, "-t")) {
			measure_epr(m, beta, 10000000, t_node, t_reduc);
		}else{
			LLL_XD(m, 0.99);		
			BKZ_FP(m, 0.99, beta);		
		}
	}

	cout << m << endl << beta << endl;
	if(cmd_option_exists(argv, argv+argc, "-t")) {
		cout << t_node << endl;
		cout << t_reduc << endl;
	}

	return 0;
}

void measure_epr(mat_ZZ& basis, int beta, unsigned long termination, double& t_node, double& t_reduc) {
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

	//cout << "Using Darmstadt bound: " << gh << endl;

	double* Rvec= new double[dimension];
	for(int i= 0; i < dimension; i++)
		conv(Rvec[i], gh*gh);

	t_reduc= (double)(end - begin) / CLOCKS_PER_SEC;
	nodes= termination;
	enumerate_epr(mu, c, Rvec, dimension, result, nodes, enu_time);

	/*if(enu_time < 1.0/CLOCKS_PER_SEC/termination) 
		cerr << "Warning! Zero enumeration time!" << endl;
	else if (nodes!=0)
		cerr << "Warning! Nodes traversed: " << (termination-nodes) << endl;*/
	t_node= enu_time;

}

