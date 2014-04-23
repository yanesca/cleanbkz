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
#include <cleanbkz/cjloss.hpp>
#include <NTL/LLL.h>

using namespace std;

extern double t_extreme_reference(double boundary[], double b_star_norm[], double t_node, double t_reduc, int n); 
extern double t_extreme(double boundary[], double b_star_norm[], double t_node, double t_reduc, int n); 

extern double polytope_volume(double vols[], double bounds[], int dim);

int main(int argc, char** argv) {
	double* vols= new double[15];

	// Unit tests for the polytope volume computation. Reference values were computed with vinci.
	//TODO: dimension 15 multiple

	cout << "Testing dimension 5:\t";
	double bounds5[]= { 0.1, 0.2, 0.3, 0.4, 0.5 };
	double vinci_vol5= 1.08e-04;
	if (polytope_volume(vols, bounds5, 5)-vinci_vol5 < 0.01e-04)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;
	
	cout << "Testing dimension 6:\t";
	double bounds6[]= { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
	double vinci_vol6= 2.334305555556e-05;
	if (polytope_volume(vols, bounds6, 6)-vinci_vol6 < 1e-17)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;
	
	cout << "Testing dimension 9:\t";
	double bounds9[]= { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
	double vinci_vol9= 2.755731922399e-07;
	if (polytope_volume(vols, bounds9, 9)-vinci_vol9 < 1e-19)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;

	cout << "Testing dimension 10:\t";
	double bounds10[]= { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95 };
	double vinci_vol10= 5.120005762235e-08;
	if (polytope_volume(vols, bounds10, 10)-vinci_vol10 < 1e-20)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;


	// Unit tests for the enumeration running time estimator function
	//TODO: cjloss destruktort megirni és itt meghivni
	//TODO: setupot és teardownt irni a teszthez

/*	
	cjloss l(80, 0.94, 0);
	BKZ_QP1(l.basis, 0.99, 2);		

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
			conv(c[i], SqrRoot(c1[i]));

	double* boundary= new double[mu1.NumRows()];	
	boundary[0]= boundary[1]= 2;
	for(int i= 2; i < mu1.NumRows()-1; i+=2)
		boundary[i]= boundary[i+1]= boundary[i-1]+boundary[0];
	boundary[mu1.NumRows()-1]= mu1.NumRows();

	double t_node= 3.47193e-08;
	double t_reduc= 0.101471;
	if(t_extreme(boundary, c, t_node, t_reduc, mu1.NumRows())-t_extreme_reference(boundary, c, t_node, t_reduc, mu1.NumRows()) < 1e+62)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;

	cout << "Production: " << t_extreme(boundary, c, t_node, t_reduc, mu1.NumRows()) << endl;
	cout << "Reference: " << t_extreme_reference(boundary, c, t_node, t_reduc, mu1.NumRows()) << endl;
*/
	return 0;
}
