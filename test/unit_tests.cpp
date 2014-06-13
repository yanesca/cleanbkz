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

extern RR fact_RR(int x);

extern RR ball_vol_RR(int k, RR r);

extern RR RR_PI;
extern double ball_vol(int k, double r); 

extern double integral_even(int ltilde, int l, double tvec[], double vvec[]);
extern RR integral_even_RR(int h, int l, RR tvec[], RR vvec[]); 

extern double integral_odd(int ltilde, int l, double tvec[], double vvec[]);
extern RR integral_odd_RR(int h, int l, RR tvec[], RR vvec[]); 

void genbounds(int l, double* tvec);

int main(int argc, char** argv) {

	/* Unit tests for the polytope volume computation. Reference values are simplices or computed with vinci. */

	double* vols_d= new double[70];

	bool vinci_test= true;
	double vinci_vol_d;
	double* bounds_d= new double[40]; 

	bounds_d[0]= 0.1; bounds_d[1]= 0.2; bounds_d[2]= 0.3;  bounds_d[3]= 0.4;  bounds_d[4]= 0.5;
	vinci_vol_d= 1.08e-04;
	if (abs(integral_even(5,5,bounds_d,vols_d)-vinci_vol_d) > 1e-19) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 5" << endl;
		cout << integral_even(5,5,bounds_d,vols_d) << endl; 
		cout << vinci_vol_d << endl; 
	}

	bounds_d[0]= 0.1; bounds_d[1]= 0.2; bounds_d[2]= 0.3;  bounds_d[3]= 0.4;  bounds_d[4]= 0.5; bounds_d[5]= 0.6;
	vinci_vol_d= 2.334305555556e-05;
	if (abs(integral_even(6,6,bounds_d,vols_d)-vinci_vol_d) > 1e-17) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 6" << endl;
		cout << integral_even(6,6,bounds_d,vols_d) << endl; 
		cout << vinci_vol_d << endl; 
	}

	bounds_d[0]= 0.1; bounds_d[1]= 0.2; bounds_d[2]= 0.3;  bounds_d[3]= 0.4;  bounds_d[4]= 0.5; bounds_d[5]= 0.6;
	bounds_d[6]= 0.7; bounds_d[7]= 0.8; bounds_d[8]= 0.9;
	vinci_vol_d= 2.755731922399e-07;
	if (abs(integral_even(9,9,bounds_d,vols_d)-vinci_vol_d) > 1e-19) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 9" << endl;
		cout << integral_even(9,9,bounds_d,vols_d) << endl; 
		cout << vinci_vol_d << endl; 
	}

	bounds_d[0]= 0.1; bounds_d[1]= 0.2; bounds_d[2]= 0.3;  bounds_d[3]= 0.4;  bounds_d[4]= 0.5; bounds_d[5]= 0.6;
	bounds_d[6]= 0.7; bounds_d[7]= 0.8; bounds_d[8]= 0.9; bounds_d[9]= 0.95;
	vinci_vol_d= 5.120005762235e-08;
	if (abs(integral_even(10,10,bounds_d,vols_d)-vinci_vol_d) > 1e-20) {
		vinci_test= false;
		cout << "\tVinci test FAILED in dimension 10" << endl;
		cout << integral_even(10,10,bounds_d,vols_d) << endl; 
		cout << vinci_vol_d << endl; 
	}

	if(vinci_test)
		cout << "\tVinci test PASSED." << endl;

	RR* vols= new RR[71];

	int prec= RR_PRECISION;	

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
		pv= integral_even_RR(i, i, simplex, vols);
		sv= 1/fact_RR(i);
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

	RR* linear= new RR[70];
	double* linear_d= new double[70];
	for(int i= 0; i < 70; i++)
		linear[i]= linear_d[i]= (i+1)/70.0;

	RR rv, dv;	
	bool double_test= true;
	for(int i= 1; i < 40; i++) {
		rv= integral_even_RR(i, i, linear, vols);
		dv= integral_even(i, i, linear_d, vols_d);
		y= rv.exponent()+175;
		pow(epsilon, x, y);
		// Checking if relative error is greater than prec-175 binary digits
		if( abs(rv - dv) > epsilon) {
			double_test= false;
			cout << "\tDouble test (even) failed in dimension " << i << endl;
			}
		/*cout << "Dimension " << i << endl;
		cout << "\t RR: " << rv << endl;
		cout << "\t double: " << dv << endl;
		cout << "\t diff: " << abs(rv-dv) << endl;
		cout << "\t epsilon: " << epsilon << endl << endl;*/
		}
	
	if(double_test)
		cout << "\tDouble test (even) PASSED." << endl;

	double_test= true;
	for(int i= 1; i < 40; i++) {
		rv= integral_odd_RR(i, i, linear, vols);
		dv= integral_odd(i, i, linear_d, vols_d);
		y= rv.exponent()+155;
		pow(epsilon, x, y);
		// Checking if relative error is greater than prec-155 binary digits
		if( abs(rv - dv) > epsilon) {
			double_test= false;
			cout << "\tDouble test (odd) failed in dimension " << i << endl;
			}
		/*cout << "Dimension " << i << endl;
		cout << "\t RR: " << rv << endl;
		cout << "\t double: " << dv << endl;
		cout << "\t diff: " << abs(rv-dv) << endl;
		cout << "\t epsilon: " << epsilon << endl << endl;*/
		}
	
	if(double_test)
		cout << "\tDouble test (odd) PASSED." << endl;


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
		if( abs(ball_vol_RR(i, r)-ball_vol(i, r_d)) > epsilon) {
			ball_test= false;
			cout << "\tBall test failed in dimension " << i << endl;
		}
	}

	if(ball_test)
		cout << "\tBall test PASSED." << endl;

	return 0;
}

