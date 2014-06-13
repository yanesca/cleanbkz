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

#include <cmath>
#include <ctime>
#include <iostream>
#include <cleanbkz/cjloss.hpp>
#include <NTL/LLL.h>

using namespace std;

cjloss::cjloss(long dimension, double density, long seed) {
	values.SetLength(dimension);
	solution.SetLength(dimension);
	// computing the maximal value to reach the prescribed density
	ZZ max;
	conv(max,floor(pow(2, dimension/density)));

	ZZ s;
	conv(s,seed); 
	SetSeed(s);
	randomize(dimension,max);
	while(!check())
		randomize(dimension,max);

}

mat_ZZ cjloss::get_basis(int k) {
	mat_ZZ basis, ret;
	
	basis.SetDims(values.length()+1, values.length()+1);
	for(int i= 0; i<basis.NumRows(); i++)
		for(int j= 0; j<basis.NumCols(); j++)
			basis[i][j]= 0;

	ZZ N;
	conv(N,ceil(sqrt(basis.NumRows()-1)));
	for(int i= 0; i<basis.NumRows()-1; i++){
		basis[i][i]= 2;
		basis[i][basis.NumRows()-1]= 2*values[i]*N;
		basis[basis.NumCols()-1][i]= 1;
		}
		
	basis[basis.NumRows()-1][basis.NumCols()-1]= 2*sum*N;
	
	if(k<2)
		ret= basis;
	else {
		BKZ_QP1(basis, 0.99, k); 

		ret.SetDims(basis.NumRows()-1, basis.NumRows()-1);

		for(int i= 0; i<ret.NumRows(); i++)
			for(int j= 0; j<ret.NumCols(); j++)
				ret[i][j]= basis[i][j];
	}

	return ret;
}


/* generates a random knapsack problem with equally many zeroes and ones in the solution (condition applied in the extreme pruning article) */ 
void cjloss::randomize(long dimension, ZZ max){
	values[0]= max;
	solution[0]= 0;
	// Generating random values 
	for(int i= 1; i< dimension; i++){
		do {
			values[i]= RandomBnd(max);
		} while(values[i]==0);
		solution[i]= 0;
	}	

	// Generating a solution vector 
	int ones= dimension/2; 
	long index;
	while(ones>0){
		index= RandomBnd(dimension);
		if(solution[index]==0) {
			solution[index]= 1;
			ones--;	
		}	
	}

	// Computing the sum 
	sum= 0;
	for(int i= 0; i< dimension; i++)
		sum+= solution[i]*values[i];
	}

/* Checks if t/n < s < (n-1)t/n holds. This is condition sugested by the cjloss article. */ 
bool cjloss::check() {
	ZZ right;	
	ZZ first;
	ZZ second; 
	
	first= 0;
	for(int i= 0;i<values.length();i++)
		first+=values[i];	
	
	second= (values.length()-1)*first;
	right= values.length()*sum;
	if(first>right || second<right)
		return false;

	return true;
} 

double cjloss::get_density() const {
	ZZ max;	

	max= 0;	
	for(int i= 0; i<values.length(); i++)
		if(max<values[i])
			max= values[i];
	
	double dmax;
	conv(dmax,max);
	return values.length()/log2(dmax);
}

ostream& operator<<(ostream& os, const cjloss& obj) {
	//os << obj.basis << endl << endl;

	cout << "Density: " << obj.get_density() << endl;
	cout << "Values: " << obj.values << endl;
	cout << "Sum: " << obj.sum << endl;

	return os;
}

/*void cjloss::print_solution(const vec_RR& shortest){
	double* sol= new double[basis.NumCols()];
	double dbase, dsol;

	if(basis.NumRows()!=shortest.length()) {
		cout << "Invalid solution length: " <<  shortest.length()  << endl;
		cout << "Number of basis vectors: " << basis.NumRows()  << endl;
		return;
	}

	for(int j= 0; j< basis.NumCols(); j++)	
		sol[j]= 0;

	for(int i= 0; i< basis.NumRows(); i++)
		for(int j= 0; j< basis.NumCols(); j++) {	
			conv(dbase, basis[i][j]);
			conv(dsol, shortest[i]);
			sol[j]+= dbase*dsol;
			}

	cout << "Solution coordinates in the cjloss basis: " << shortest << endl;
	cout << "Scaled coordinates in the cjloss lattice: [";
	for(int i= 0; i< basis.NumCols()-1 ; i++) 
		cout << sol[i] << " "; 	
	cout << sol[basis.NumCols()-1] << "]" << endl; 
	
	double sqrdLength= 0;
	for(int i= 0; i< basis.NumCols(); i++) {
		sqrdLength+= sol[i]*sol[i];
		sol[i]= -(sol[i]-1)/2;
	}
 
	cout << "Squared length (found/solution): " << sqrdLength << " / " << basis.NumRows()-1 << endl << endl; 

	cout << "Original solution: " << solution << endl;
	cout << "Solution: ["; 
	for(int i= 0; i< basis.NumCols()-2 ; i++) 
		cout << sol[i] << " "; 	
	cout << sol[basis.NumCols()-2] << "]" << endl; 
}*/

