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
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/enumeration.hpp>
#include <NTL/LLL.h>
#include <NTL/new.h>

using namespace std;

unsigned long count_vectors_in_ball(double** mu, double *b, double R, int n); 
extern double ball_vol(int k, double r);

int main(int argc, char** argv) {
	int dim, bsize= 2;
	stringstream ss;
	double roverl;
	int seed= 0;

	ss << argv[1];
	ss >> dim;
	ss.clear();

	ss << argv[2];
	ss >> roverl;
	ss.clear();

	ss << argv[3];
	ss >> seed;
	ss.clear();


	cjloss l(dim, 0.94, seed);
//	cjloss l(dim, 0.94, time(NULL));
	double lambda= sqrt(dim);

	//cout << " " << l.get_basis(1) << endl;

	mat_ZZ basis= l.get_basis(1);

	mat_RR mu1;
	vec_RR c1;
	BKZ_QP1(basis, 0.99, bsize); 
	ComputeGS(basis,mu1,c1);

	double* cs= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(cs[i], sqrt(c1[i]));


// Gaussian heuristic
	double determinant= 1;
	for(int i= 0; i < mu1.NumRows(); i++) 
		determinant*= cs[i];
	//double gh= pow(determinant/ball_vol(dim, 1),1.0/dim);
	//cout << "# Shortest vector length: " << sqrt(dim-1) << endl;
	//cout << "# Gaussian heuristic: " << gh << endl;
	//cout << "# GH/lambda Ratio: " << gh/sqrt(dim-1) << endl;

// Enumeration

	//double R= lambda*roverl;
	double R= lambda;
	double** mu= new double*[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		mu[i]= new double[mu1.NumCols()]; 
	for(int i= 0; i < mu1.NumRows(); i++)
		for(int j= 0; j < mu1.NumCols(); j++)
			if(i==j)
				mu[i][j]= 1;
			else
				conv(mu[i][j], mu1[i][j]);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		conv(c[i], c1[i]);

	double* prune= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		prune[i]= R*R;

	vec_RR result;
	basis= l.get_basis(1);
	enumerate_epr(basis, bsize, prune, result);
	//cout << "Shortest vector: " << result << endl;
	//cout << "Reduced basis: " << endl << basis << endl;

	cout << "\"Enumeration\"" << endl;
	for(int i= 1; i<=roverl; i++) {
		R= lambda*i;
		//cout << "Number of lattice vectors inside of an " << R << " radius ball: " << count_vectors_in_ball(mu, c, R, dim) << endl;
		cout << i << " " <<  count_vectors_in_ball(mu, c, R, dim) << endl;
		}
		
		cout << endl << endl;

	cout << "\"Gaussian Heuristic\"" << endl;
	for(int i= 1; i<=roverl; i++) {
		R= lambda*i;
		//cout << "Gaussian heuristic: " << ball_vol(dim,R)/determinant << endl;
		cout << i << " " << ball_vol(dim,R)/determinant << endl;
		}

}

unsigned long count_vectors_in_ball(double** mu, double *b, double R, int n) {
	unsigned long* prof= new unsigned long[n];
	for(int i= 0; i < n; i++)
		prof[i]= 1;

	unsigned long count= 1;

	double* Rvec= new double[n];
	for(int i= 0; i < n; i++)
		Rvec[i]= R*R;

	double** sigmamat= new double*[n+1];
	for(int i= 0; i < n+1; i++)
		sigmamat[i]= new double[n]; 
	for(int i= 0; i < n+1; i++)
		for(int j= 0; j < n; j++)
			sigmamat[i][j]= 0;

	int* rvec= new int[n+1];
	for(int i= 0; i < n+1; i++)
		rvec[i]= i-1;

	double* rhovec= new double[n+1];
	for(int i= 0; i < n+1; i++)
		rhovec[i]= 0;
		
	int* vvec= new int[n];	
	vvec[0]= 1;
	for(int i= 1; i < n; i++)
		vvec[i]= 0;

	double* cvec= new double[n];	
	for(int i= 0; i < n; i++)
		cvec[i]= 0;

	int* wvec= new int[n];	
	for(int i= 0; i < n; i++)
		wvec[i]= 0;

	int last_nonzero= 0;

	int k= 0;

	while(true) {
		rhovec[k]= rhovec[k+1]+(vvec[k]-cvec[k])*(vvec[k]-cvec[k])*b[k];	

		if((k!=0) && (rhovec[k] <= Rvec[n-k-1])) {

			prof[n-k-1]++;

			k--;
			rvec[k]= rvec[k]>rvec[k+1]?rvec[k]:rvec[k+1];		
			for(int i= rvec[k+1]; i>k; i--)
				sigmamat[i][k]= sigmamat[i+1][k] + vvec[i]*mu[i][k];
			cvec[k]= -sigmamat[k+1][k];

			vvec[k]= lround(cvec[k]);
			wvec[k]= 1;

		} else {
			if((k==0) && (rhovec[0] <= Rvec[n-1])){ 
				prof[n-k-1]++;
				count++;
				}

			k++;
			if(k==n)
				break;
			rvec[k]= k;

			if(k>=last_nonzero){
				last_nonzero= k;
				vvec[k]++;
			} else {
				if(vvec[k]>cvec[k])
					vvec[k]-= wvec[k];
				else
					vvec[k]+= wvec[k];
				wvec[k]++;
			}
		}
	}

	delete [] Rvec;
	for(int i= 0; i < n+1; i++)
		delete [] sigmamat[i]; 
	delete [] sigmamat;
	delete [] rvec; 
	delete [] vvec;
	delete [] wvec;
	delete [] cvec;
	delete [] rhovec;

	return count;
}
