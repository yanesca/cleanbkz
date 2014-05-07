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
double n_full_gsa(double R, double b1_norm, double alpha, int n); 
double ball_vol(int k, double r);

int main(int argc, char** argv) {
	/*double t_node= 6.19e-08;
	double t_reduc= 43.248;
	cjloss l(80, 0.94, 0);

	cout << "Computing reduced basis." << endl;
	BKZ_QP1(l.basis, 0.99, 30);		

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);
	
	RR* c= new RR[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) {
		c[i]= SqrRoot(c1[i]);
		c[i].SetPrecision(RR_PRECISION);
	}

	int dim= mu1.NumRows();
	double* boundary= new double[dim];	

	double R;
	RR uBound;
	uBound= dim-1;
	ceilPrec(uBound, SqrRoot(uBound), 53);
	conv(R, uBound);
	cout << R << ", " << R*R << endl;

	cout << "Computing boundary function.";
	double p_succ;
	double t_enum;	
	generate_boundary(c, t_node, t_reduc, dim, boundary, R, 0.1, 5000, p_succ, t_enum, true); 
	cout << " Done." << endl;
	cout << "\tSuccess probability: " << p_succ << endl;
	cout << "\tPredicted enumeration time: " << t_enum << endl << endl;*/

/*	ifstream basis_file("s20d80k30.lat");	
	ifstream bound_file("s20d80k30.bnd");	

	mat_ZZ basis;
	basis_file >> basis;
	//cout << "Succesfully read basis: " << endl << basis << endl;

	vec_RR bounds;
	bound_file >> bounds;	
	//cout << "Succesfully read boundary: " << endl << bounds << endl;

	double* bounds_d= new double[bounds.length()];
	for(int i= 0; i < bounds.length(); i++)
		conv(bounds_d[i], bounds[i]*bounds[i]);
	
	cout << "Performing pruned enumeration." << endl;
	
	vec_RR solution;	
	enumerate_epr(basis, 30, bounds_d, solution);

	cout << "Solution length: " << solution.length() << endl;*/

	int dim;
	stringstream ss;

	ss << argv[1];
	ss >> dim;
	ss.clear();

	cjloss l(dim, 0.94, 10);

// Enumeration
	double R= dim-1;
	cout << "# Shortest vector length: " << R << endl;

	double* act= new double[dim];
	for(int i= 0; i < dim; i++) 
		act[i]= R;

	vec_RR solution;	
	enumerate_ntl(l.basis, 2, act, solution);
	enumerate_epr(l.basis, 2, act, solution);
	cout << "# Solution length: " << solution.length() << endl << endl << endl;


// Prediction	
	for(int i= 0; i < dim; i++) 
		act[i]= sqrt(act[i]);

	mat_RR mu1;
	vec_RR c1;
	//BKZ_QP1(l.basis, 0.99, 2); 
	ComputeGS(l.basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

	t_extreme_reference(act, c, 1, 1, dim);

/*	
	cout << "# Rank: " << BKZ_QP1(l.basis, 0.99, 2) << endl; 

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);

	//cout << l.basis << endl;

	double gh= 1;
	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) {
		conv(c[i], sqrt(c1[i]));
		gh*= c[i];
		//cout << i << " " << c[i] << endl;
		}
	gh= pow(gh/ball_vol(dim, 1),1.0/dim);
	cout << "# Gaussian heuristic: " << gh << endl;

	
	conv(R, sqrt(c1[0]));
	for(int i= 0; i < dim; i++) 
		act[i]= R ;

	//t_extreme_reference(act, c, t_node, t_reduc, dim); 
	//t_extreme(act, c, 1, 1, dim); 
	t_extreme_reference(act, c, 1, 1, dim); 


	//Full test
//	cout << "# Performing full test" << endl;

	cout << "# nfull: " << n_full(act[dim-1], c, dim) << endl << endl << endl; 
//	cout << "# nfull (gsa): " << n_full_gsa(act[dim-1], c[0], 1/(1.0219*1.0219), dim) << endl; 

//	t_extreme_reference(act, c, 1, 1, dim);
	
	cout << "# BKZ 30 rank: " << BKZ_QP1(l.basis, 0.99, 30) << endl; 
	
	ComputeGS(l.basis,mu1,c1);

	//cout << l.basis << endl;

	for(int i= 0; i < dim; i++) 
		act[i]= sqrt(act[i]);

	for(int i= 0; i < mu1.NumRows(); i++) {
		conv(c[i], sqrt(c1[i]));
		//cout << i << " " << c[i] << endl;
		}

	t_extreme_reference(act, c, 1, 1, dim);*/

/*	conv(R, c1[0]);
	for(int i= 0; i < dim; i++) 
		act[i]= R;

	enumerate_epr(l.basis, 2, act, solution);
	cout << "# Solution length: " << solution.length() << endl << endl << endl;*/
}
