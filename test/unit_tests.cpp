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
#include <cleanbkz/boundary.hpp>
#include <NTL/LLL.h>

using namespace std;

extern RR t_extreme_reference_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n, int prec);
extern RR t_extreme_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n, int prec);

extern double t_extreme_reference(double boundary[], double b_star_norm[], double t_node, double t_reduc, int n); 
extern double t_extreme(double boundary[], double b_star_norm[], double t_node, double t_reduc, int n); 

extern double polytope_volume(double vols[], double bounds[], int dim);

inline double fact(double x) {
	return (x < 2 ? 1 : x * fact(x - 1));
}

extern RR polytope_volume_RR(RR vols[], RR bounds[], int dim, int prec);

extern RR fact_RR(int x, int prec);

extern RR ball_vol_RR(int k, RR r, int prec);

extern RR RR_PI;
extern double ball_vol(int k, double r); 

int main(int argc, char** argv) {
	RR* vols= new RR[71];

	int prec= RR_PRECISION;	

	/* Unit tests for the polytope volume computation. Reference values are simplices or computed with vinci. */

	cout << "Testing polytope volume computation:" << endl;

	bool vinci_test= true;
	RR vinci_vol;
	vinci_vol.SetPrecision(prec);
	RR* bounds= new RR[15]; 

	bounds[0]= 0.1; bounds[1]= 0.2; bounds[2]= 0.3;  bounds[3]= 0.4;  bounds[4]= 0.5;
	vinci_vol= 1.08e-04;
	if (abs(polytope_volume_RR(vols, bounds, 5, prec)-vinci_vol) > 1e-19) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 5" << endl;
	}	
	
	bounds[0]= 0.1; bounds[1]= 0.2; bounds[2]= 0.3;  bounds[3]= 0.4;  bounds[4]= 0.5; bounds[5]= 0.6;
	vinci_vol= 2.334305555556e-05;
	if (abs(polytope_volume_RR(vols, bounds, 6, prec)-vinci_vol) > 1e-17) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 6" << endl;
	}	
	
	bounds[0]= 0.1; bounds[1]= 0.2; bounds[2]= 0.3;  bounds[3]= 0.4;  bounds[4]= 0.5; bounds[5]= 0.6; 
	bounds[6]= 0.7; bounds[7]= 0.8; bounds[8]= 0.9;
	vinci_vol= 2.755731922399e-07;
	if (abs(polytope_volume_RR(vols, bounds, 9, prec)-vinci_vol) > 1e-19) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 9" << endl;
	}	

	bounds[0]= 0.1; bounds[1]= 0.2; bounds[2]= 0.3;  bounds[3]= 0.4;  bounds[4]= 0.5; bounds[5]= 0.6; 
	bounds[6]= 0.7; bounds[7]= 0.8; bounds[8]= 0.9; bounds[9]= 0.95;
	vinci_vol= 5.120005762235e-08;
	if (abs(polytope_volume_RR(vols, bounds, 10, prec)-vinci_vol) > 1e-20) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 10" << endl;
	}	

	//TODO: dimension 15 multiple

	if(vinci_test)
		cout << "\tVinci test PASSED." << endl;

	RR epsilon, x, y;
	epsilon.SetPrecision(prec);
	x= 2; 

	RR one;
	one.SetPrecision(prec);
	one= 1;
	
	RR* simplex= new RR[70];
	for(int i= 0; i < 70; i++)
		simplex[i]= one;

	RR pv, sv;	
	bool simplex_test= true;
	for(int i= 1; i < 71; i++) {
		pv= polytope_volume_RR(vols, simplex, i, prec);
		sv= 1/fact_RR(i, prec);
		y= pv.exponent()+95;
		pow(epsilon, x, y);
		// Checking if relative error is greater than prec-95 binary digits
		if( abs(pv - sv) > epsilon) {
			simplex_test= false;
			cout << "\tSimplex test failed in dimension " << i << endl;
			}
		/*cout << "Dimension " << i << endl;
		cout << "\t pol_vol: " << polytope_volume_RR(vols, simplex, i, prec) << endl;
		cout << "\t simplex formula: " << 1/fact_RR(i, prec) << endl << endl;
		cout << "\t diff: " << abs(polytope_volume_RR(vols, simplex, i, prec) - 1/fact_RR(i, prec)) << endl;
		cout << "\t epsilon: " << epsilon << endl;*/
		}

	if(simplex_test)
		cout << "\tSimplex test PASSED." << endl;

	// Unit tests for the n dimensional ball volume computation
	cout << "Testing ball volume computation:" << endl;

	bool ball_test= true;
	double r_d= 1;
	RR r;	
	r= r_d;
	RR_PI.SetPrecision(prec);
	RR_PI= ComputePi_RR();

	y= -50;
	pow(epsilon, x, y);
	for(int i= 1; i < 71; i++) {
		if( abs(ball_vol_RR(i, r, prec)-ball_vol(i, r_d)) > epsilon) {
			ball_test= false;
			cout << "\tBall test failed in dimension " << i << endl;
		}
	}

	if(ball_test)
		cout << "\tBall test PASSED." << endl;

	// Unit tests for node number estimation

	/* Unit tests for the enumeration running time estimator function
	//TODO: cjloss destruktort megirni és itt meghivni
	//TODO: setupot és teardownt irni a teszthez

		
	cjloss l(80, 0.94, 0);
	BKZ_QP1(l.basis, 0.99, 2);		

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(l.basis,mu1,c1);

	double* c_d= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		conv(c_d[i], SqrRoot(c1[i]));

	double* boundary_d= new double[mu1.NumRows()];	
	boundary_d[0]= boundary_d[1]= 2;
	for(int i= 2; i < mu1.NumRows()-1; i+=2)
		boundary_d[i]= boundary_d[i+1]= boundary_d[i-1]+boundary_d[0];
	boundary_d[mu1.NumRows()-1]= mu1.NumRows();

	RR* c= new RR[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) {
		c[i].SetPrecision(prec);
		c[i]= SqrRoot(c1[i]);
	}

	RR* boundary= new RR[mu1.NumRows()];	
	boundary[0].SetPrecision(prec);
	boundary[1].SetPrecision(prec);
	boundary[0]= boundary[1]= 2;
	for(int i= 2; i < mu1.NumRows()-1; i+=2) {
		boundary[i].SetPrecision(prec);
		boundary[i]= boundary[i+1]= boundary[i-1]+boundary[0];
	}
	boundary[mu1.NumRows()-1]= mu1.NumRows();

	double t_node= 3.47193e-08;
	double t_reduc= 0.101471;
	if(t_extreme(boundary, c, t_node, t_reduc, mu1.NumRows())-t_extreme_reference(boundary, c, t_node, t_reduc, mu1.NumRows()) < 1e-10)
		cout << "PASSED" << endl;
	else cout << "FAILED" << endl;

	//cout << "Production (double): " << t_extreme(boundary_d, c_d, t_node, t_reduc, mu1.NumRows()) << endl;
	cout << "Reference: " << t_extreme_reference_RR(boundary, c, t_node, t_reduc, mu1.NumRows(), prec) << endl; 
	cout << "Reference (double): " << t_extreme_reference(boundary_d, c_d, t_node, t_reduc, mu1.NumRows()) << endl; 
	cout << "Production: " << t_extreme(boundary, c, t_node, t_reduc, mu1.NumRows(), prec) << endl;
	cout << "Diff: " << t_extreme(boundary, c, t_node, t_reduc, mu1.NumRows(), prec) - t_extreme_reference_RR(boundary, c, t_node, t_reduc, mu1.NumRows(), prec) << endl; */

	return 0;
}
