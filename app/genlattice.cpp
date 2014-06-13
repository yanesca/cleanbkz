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
#include <sstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/version.hpp>
#include <cleanbkz/tools.hpp>
#include <NTL/LLL.h>

using namespace std;

//TODO: refaktorálni a paraméterfeldolgozást: csinálni hozzá int és double változatot, hibaüzenetekkel és kivételkezeléssel

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
			<< "\t-k n\t\tReduces the generated basis with blocksize of n." << endl;
		return 0;
	}

	if (!cmd_option_exists(argv, argv+argc, "-d")) {
		cout << "Lattice dimension is missing. Run again with -h for help." << endl;
		return 0;
		}

	if (!cmd_option_exists(argv, argv+argc, "-m") && !cmd_option_exists(argv, argv+argc, "-c")) {
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
	
	if(beta > 1 && !(cmd_option_exists(argv, argv+argc, "-c")))	
		BKZ_QP1(m, 0.99, beta);		

	cout << m << endl << beta;

	return 0;
}
