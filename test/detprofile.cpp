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
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/enumeration.hpp>
#include <cleanbkz/tools.hpp>
#include <NTL/LLL.h>

using namespace std;

double t_extreme_reference(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n); 
double t_extreme(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n); 
double n_full(double R, double b_star_norm[], int n);
double n_full_gsa(double R, double scale, int bsize, int n);
double ball_vol(int k, double r);

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim);
ZZ autoscale(mat_ZZ& basis, double Rvec[]);

int main(int argc, char** argv) {
	int dim, bsize= 20;
	stringstream ss;

	ss << argv[1];
	ss >> dim;
	ss.clear();


//	cjloss l(dim, 0.94, 10);
//	cjloss l(dim, 0.94, time(NULL));

//	cout << "# " << l.basis << endl;

	ZZ determinant;
//	GenPrime(determinant,dim);
	GenPrime(determinant,floor(dim));
	mat_ZZ basis;
	gen_randlat(basis,determinant,dim); 

	//basis= l.basis;

	mat_RR mu1;
	vec_RR c1;
	BKZ_QP1(basis, 0.99, bsize); 
	ComputeGS(basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

// Gaussian heuristic
	double gh= 1;
	double* gsghs= new double[dim]; 
	for(int i= 0; i < mu1.NumRows(); i++) {
		gh*= c[i];
		gsghs[i]= pow(gh/ball_vol(i+1, 1),1.0/(i+1));
		}
	gh= pow(gh/ball_vol(dim, 1),1.0/dim);
	//cout << "# Shortest vector length: " << sqrt(dim-1) << endl;
	//cout << "# Gaussian heuristic: " << gh << endl;
	//cout << "# GH/lambda Ratio: " << gh/sqrt(dim-1) << endl;

// Bounding functions
	//double R= dim-1;
	double R= gh*gh;

	double* full= new double[dim];
	for(int i= 0; i < dim; i++) 
		full[i]= R;

	double* sh= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		sh[2*i]= sh[2*i+1]= MIN(1,R*sqrt(1.05*2*i/dim));
	if(dim%2==0)
		sh[dim-1]= R;	
	else
		sh[dim-2]= sh[dim-1]= R;	

	double* linear= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		linear[2*i]= linear[2*i+1]= 2*(i+1)*R/dim;
	if(dim%2==0)
		linear[dim-1]= R;	
	else
		linear[dim-2]= linear[dim-1]= R;	

	double* step= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		if(i<dim/3)
			step[2*i]= step[2*i+1]= R/2;
		else
			step[2*i]= step[2*i+1]= R;	
	if(dim%2==0)
		step[dim-1]= R;	
	else
		step[dim-2]= step[dim-1]= R;	

	double* min= new double[dim];
	for(int i= 0; i < dim/2; i++) 
		min[2*i]= min[2*i+1]= pow(gsghs[2*i+1],2);

	double* act= linear;

	//cout << "# Scale: " << autoscale(basis, act) << endl;

// Enumeration
	vec_RR solution;	
	//enumerate_ntl(l.basis, bsize, full, solution);
	enumerate_epr(basis, bsize, act, solution);
	cout << "# Solution length: " << solution.length() << endl << endl << endl;


// Prediction	

	ComputeGS(basis,mu1,c1);

	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

	for(int i= 0; i < dim; i++) { 
		full[i]= sqrt(full[i]);
		step[i]= sqrt(step[i]);
		linear[i]= sqrt(linear[i]);
		min[i]= sqrt(min[i]);
		}

	t_extreme_reference(act, c, 1, 1, dim);


// GSA	
	cout << endl << endl;
	double alpha=1;
	for(int i= 0; i < mu1.NumRows(); i++) 
		alpha*= c[i];
	alpha= log(pow(alpha,2.0/dim));
//	n_full_gsa(gh, alpha, bsize, dim); 


/*	
	cout << "# Rank: " << BKZ_QP1(l.basis, 0.99, 2) << endl; 

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);

	//cout << l.basis << endl;

	
	conv(R, sqrt(c1[0]));
	for(int i= 0; i < dim; i++) 
		full[i]= R ;

	//t_extreme_reference(full, c, t_node, t_reduc, dim); 
	//t_extreme(full, c, 1, 1, dim); 
	t_extreme_reference(full, c, 1, 1, dim); 


	//Full test
//	cout << "# Performing full test" << endl;

	cout << "# nfull: " << n_full(full[dim-1], c, dim) << endl << endl << endl; 
//	cout << "# nfull (gsa): " << n_full_gsa(full[dim-1], c[0], 1/(1.0219*1.0219), dim) << endl; 

//	t_extreme_reference(full, c, 1, 1, dim);
	
	cout << "# BKZ 30 rank: " << BKZ_QP1(l.basis, 0.99, 30) << endl; 
	
	ComputeGS(l.basis,mu1,c1);

	//cout << l.basis << endl;

	for(int i= 0; i < dim; i++) 
		full[i]= sqrt(full[i]);

	for(int i= 0; i < mu1.NumRows(); i++) {
		conv(c[i], sqrt(c1[i]));
		//cout << i << " " << c[i] << endl;
		}

	t_extreme_reference(full, c, 1, 1, dim);*/

/*	conv(R, c1[0]);
	for(int i= 0; i < dim; i++) 
		full[i]= R;

	enumerate_epr(l.basis, 2, full, solution);
	cout << "# Solution length: " << solution.length() << endl << endl << endl;*/
}

ZZ autoscale(mat_ZZ& basis, double Rvec[]){
	ZZ scale;

	scale= 1;
	for(int i= 0; i < basis.NumRows(); i++)
		if(sqrt(i+1)/Rvec[i]>scale)
			conv(scale, sqrt((i+1)/Rvec[i]));

	// conversion computes floor and ceil is needed
	scale+= 1;
	if(scale == 1)
		return scale;

	for(int i= 0; i <  basis.NumRows(); i++)
		for(int j= 0; j < basis.NumCols(); j++)
			basis[i][j]*= scale;

	double rscale;
	conv(rscale, scale);
	rscale*= rscale;
	for(int i= 0; i <  basis.NumRows(); i++)
		Rvec[i]*= rscale;
	
	return scale;
}

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim) {
	basis.SetDims(dim, dim);
	
	for(int i= 0; i < dim; i++) {
		basis[i][0]= RandomBnd(determinant);
		basis[i][i]= 1;
	}

	basis[0][0]= determinant;
} 
