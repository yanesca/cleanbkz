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
#include <iostream>
#include <algorithm>
#include <cleanbkz/boundtools.hpp>
#include <cleanbkz/version.hpp>

using namespace std;

char* get_cmd_option(char** begin, char** end, const string& option) {
    char** itr= find(begin, end, option);

    if (itr != end && ++itr != end) {
        return *itr;
    }

    return 0;
}

bool cmd_option_exists(char** begin, char** end, const string& option) {
    return find(begin, end, option) != end;
}

int main(int argc, char** argv) {
  	int dimension= 0;
	int bkz_blocksize= 2;
	int nodes= 1000000;
	int samples= 1000; 
	int start= 0; 
	int end= 0;
	stringstream ss;
	char* act_arg;
	clock_t begin, finish;

	if (argc==1 || cmd_option_exists(argv, argv+argc, "-h")) {
		cout << "This program measures the parameters t_node and t_reduc required for the computation of boundary functions. These are the running time of the enumeration and reduction algorithms on the current platform. Program options:" << endl
 			<< "\t-d n\t\tMeasure running times in dimension n." << endl
 			<< "\t-s n\t\tStart the measurements in dimension n. It is ignored when -d is given" << endl
 			<< "\t-e n\t\tContinue measuring in all dimensions until dimension n with a step five. It is ignored when -d is given" << endl
			<< "\t-b n\t\tAbort enumeration after processing n nodes. (default: 10,000,000)" << endl
			<< "\t-n n\t\tNumber of experiments to make. (default: 1000)" << endl
			<< "\t-k n\t\tThe blocksize of BKZ used for preprocessing. (default: 2)" << endl;
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
		ss << act_arg;
		ss >> end;
		ss.clear();
		if(end < start) {
			cerr << "ERROR: End dimension lower than start dimension. Aborting" << endl;
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
		}
	else {
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
			cout << "Expected total running time of the experiments: " << total/3600 << "h" << (total%3600)/60 << "m" << (total%36000)%60 << "s" << endl << endl;
			}

		for(int i= start; i <= end; i+=5) {
			cout << "Measuring dimension " << i << ":" << endl;
			measure_epr(i, bkz_blocksize, samples, nodes, t_node, t_reduc, zeroes, bias); 
			cout << "t_node= " << t_node << endl;
			cout << "t_reduc= " << t_reduc << endl;
			cout << "biased/zeroes: " << bias << "/" << zeroes << endl << endl;
			}
		}

	cout << "(Sometimes the enumeration doesn't have to inspect the prescribed number of nodes or even nodes enough to make the runing time measureable, the number of these cases are in the \"biased/zeroes\" line. The first number measures only the biased cases with nonzero running time.)" << endl << endl;

	return 0;
}
