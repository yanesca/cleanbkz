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
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <NTL/RR.h>
#include <cleanbkz/boundary.hpp>

using namespace std;

inline double fact(double x) {
	return (x < 2 ? 1 : x * fact(x - 1));
}

RR fact_RR(int x, int prec) {
	RR ret;
	ret.SetPrecision(prec);
	ret= 1;

	for(int i= 2; i <= x; i++) 
		ret*= i;

	return ret;
	}

RR polytope_volume_RR(RR vols[], RR bounds[], int dim, int prec) {
	RR tmp, pwr, ret;

	tmp.SetPrecision(prec);
	pwr.SetPrecision(prec);
	ret.SetPrecision(prec);

	ret= 1;
	if(dim==0) 
		return ret;

	ret= 0;

	vols[dim-1]= polytope_volume_RR(vols, bounds, dim-1, prec);

	for(int i= 0; i < dim; i++) {
		tmp= dim-i;
		pow(pwr, bounds[i], tmp); 
		ret+= pwr*vols[i]/fact_RR(dim-i,prec)*((dim-i+1)%2==0?1:-1);
	}
		
	return ret;
}

RR RR_PI;

RR ball_vol_RR(int k, RR r, int prec) {
	RR pwr, exp, pipow;
	pwr.SetPrecision(prec);
	exp.SetPrecision(prec);
	pipow.SetPrecision(prec);
	
	exp= k;
	pow(pwr, r, exp);
 
	if (k%2==0) { 
		exp= k/2;
		pow(pipow, RR_PI, exp); 
		return pipow / fact_RR(k/2, prec) * pwr;
	} else {
		exp= k/2;
		pow(pipow, 4*RR_PI, exp);
		return 2 * pipow * fact_RR(k/2, prec) / fact_RR(k, prec) * pwr;
	}
}

double ball_vol(int k, double r) {

	/*cout << "#" << endl << "# ball_vol: R= " << r << " dim= " << k << endl;	
	cout << "# Gamma: " << pow(M_PI, k/2.0) / tgamma(k/2.0+1) * pow(r, k) << endl;
	if (k%2==0) 
		cout << "# Form: " << pow(M_PI, k/2) / fact(k/2) * pow(r, k) << endl;
	else
		cout << "# Form: " << 2 * pow(4*M_PI, k/2) * fact(k/2) / fact(k) * pow(r, k) << endl << "#" << endl;*/ 

	if (k%2==0) 
		return pow(M_PI, k/2) / fact(k/2) * pow(r, k);
	else
		return 2 * pow(4*M_PI, k/2) * fact(k/2) / fact(k) * pow(r, k);
}

double polytope_volume(double vols[], double bounds[], int dim) {
	double ret= 0;

	if(dim==0) 
		return 1;

	vols[dim-1]= polytope_volume(vols, bounds, dim-1);

	for(int i= 0; i < dim; i++)
		ret+= pow(bounds[i], dim-i)*vols[i]/fact(dim-i)*((dim-i+1)%2==0?1:-1);
		
	return ret;
}

RR t_full_RR(RR R, RR b_star_norm[], double t_node, double t_reduc, int n, int prec) {
	RR N, V_act, denom;
	N.SetPrecision(prec);
	V_act.SetPrecision(prec);
	denom.SetPrecision(prec);
	denom= 1;
	for(int i= 1; i<= n; i++) {
		V_act= ball_vol_RR(i, R, prec); 

		denom*= b_star_norm[n-i];
		N+= V_act/denom;
	}

	N/= 2;
	return t_reduc+t_node*N; 
	}

RR t_extreme_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n, long prec) {
	int pdim= (n%2==0)?n/2:n/2+1;
	RR p_succ;
	p_succ.SetPrecision(prec);
	
	RR* polytope_vols= new RR[pdim+1];
	RR* bounds= new RR[pdim];

	for(int i= 0; i< n; i+=2) {
		bounds[i/2].SetPrecision(prec);	
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);
	}

	RR N, V_act, denom;
	N.SetPrecision(prec);
	V_act.SetPrecision(prec);
	denom.SetPrecision(prec);
	denom= 1;
	polytope_vols[pdim]= polytope_volume_RR(polytope_vols, bounds, pdim, prec);
	for(int i= 1; i<= n; i++) {
		if(i%2==0)
			V_act= ball_vol_RR(i, Rvec[i-1], prec) * polytope_vols[i/2] * fact_RR(i/2, prec); 
		else
			V_act= ball_vol_RR(i, Rvec[i-1], prec) * (polytope_vols[i/2] * fact_RR(i/2,prec) + polytope_vols[i/2+1]*fact_RR(i/2+1,prec)) / 2; 

		//cout << "#  denom_" << i-1 << " = " << denom << endl;

		denom*= b_star_norm[n-i];

	/*	cout << "#  denom_" << i << " = " << denom << endl;
		cout << "#  bstarnorm_" << n-i << " = " << b_star_norm[n-i] << endl;
		cout << "#  V_" << i << " = " << V_act/denom << endl;
		cout << "#  Vball_" << i << "(" << Rvec[i-1] << ") = " << ball_vol_RR(i, Rvec[i-1], prec) << endl << endl;*/


		N+= V_act/denom;
	}

	p_succ= polytope_vols[pdim-1]*fact_RR(pdim-1, prec);

	delete [] bounds;
	delete [] polytope_vols;

	N/= 2;
	return (t_reduc+t_node*N)/p_succ; 
}

RR t_extreme(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n, long prec){ 
	return t_extreme_RR(Rvec, b_star_norm, t_node, t_reduc, n, prec);
	}

RR t_extreme_reference_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n, int prec) {
	int pdim= (n%2==0)?n/2:n/2+1;
	RR p_succ;
	p_succ.SetPrecision(prec);
	
	RR* polytope_vols= new RR[pdim+1];
	RR* bounds= new RR[pdim];

	for(int i= 0; i< n; i+=2) {
		bounds[i/2].SetPrecision(prec);	
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);
	}

	RR N, V_act, denom;
	N.SetPrecision(prec);
	V_act.SetPrecision(prec);
	denom.SetPrecision(prec);
	polytope_vols[pdim]= polytope_volume_RR(polytope_vols, bounds, pdim, prec);
	for(int i= 1; i<= n; i++) {
		if(i%2==0)
			V_act= ball_vol_RR(i, Rvec[i-1], prec) * polytope_vols[i/2] * fact_RR(i/2, prec); 
		else
			V_act= ball_vol_RR(i, Rvec[i-1], prec) * (polytope_vols[i/2] * fact_RR(i/2,prec) + polytope_vols[i/2+1]*fact_RR(i/2+1,prec)) / 2; 

		denom= 1;
		for(int j= n-i; j < n; j++) 
			denom*= b_star_norm[j];

		N+= V_act/denom;
	}

	p_succ= polytope_vols[pdim-1]*fact_RR(pdim-1, prec);

	delete [] bounds;
	delete [] polytope_vols;

	N/= 2;
	return (t_reduc+t_node*N)/p_succ; 
}


double n_full(double R, double b_star_norm[], int n) {
	double N, V_act, denom;

	N= 0;
	denom= 1;
	for(int i= 1; i<= n; i++) {
		V_act= ball_vol(i, R); 

		//cout << n-i << "bstar: " << b_star_norm[n-i] << endl;
		//denom*= b_star_norm[i-1];
		denom*= b_star_norm[n-i];
		N+= V_act/denom;
	}
	//cout << endl;

	N/= 2;
	return N; 
	}


double n_full_gsa(double R, double b1_norm, double alpha, int n) {
	double N, V_act, denom;

	N= 0;
	denom= 1;
	for(int i= 1; i<= n; i++) {
		V_act= ball_vol(i, R); 

		//cout << n-i << "gsa: " <<  b1_norm*pow(alpha,n-i) << endl;
		//cout << n-i << "alpha: " <<  alpha << " pow: " << pow(alpha,n-i) << endl;
		denom*= b1_norm*pow(alpha,n-i);
		N+= V_act/denom;
	}
	//cout << endl;

	N/= 2;
	return N; 
}

double pvol_and_scale(double Rvec[], int k) {
	double ret;
	int l= (k+1)/2-1;
	double* bounds= new double[l];
	double* polytope_vols= new double[l];


/*	cout << "khalf: " << khalf << endl;
	cout << "Rvec: [ ";
	for(int i= 0; i< 2*khalf; i++)
		cout << Rvec[i] << " ";	
	cout << "] "<< endl;*/

	for(int i= 0; i< l; i++)
		bounds[i]= Rvec[2*i]*Rvec[2*i]/(Rvec[k-1]*Rvec[k-1]);

	/*cout << "bounds: [ ";
	for(int i= 0; i< khalf; i++)
		cout << bounds[i] << " ";	
	cout << "] "<< endl << endl;*/

	ret= polytope_volume(polytope_vols, bounds, l);

	delete [] bounds;
	delete [] polytope_vols;

	return ret;
}

// The case where n is even
double t_extreme_reference(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n) {
	double N= 0;
	double V_act; 	
	double denom;
	//double N_next, N_prev;

	cout << "# T_extreme: " << endl << "# GS-lengths: [ ";
	for(int i= 0; i < n; i++)
		cout << b_star_norm[i] << " ";	
	cout << "]" << endl << "# Bounds: [ ";
	for(int i= 0; i < n; i++)
		cout << Rvec[i] << " ";	
	cout << "]" << endl << "#" << endl;

	cout << "\"Predicted\"" << endl << endl;

//	N_next= 0;
//	N_prev= 1/b_star_norm[n-1];
	for(int k= 1; k<= n; k++) {
		/*if(k%2==0) {
			//V_act= ball_vol(k, Rvec[k-1]) * pvol_and_scale(Rvec, k/2) * fact(k/2); 
			N_prev= N_next;
			N+= N_prev;
			cout << k-1 << " " << N_prev << endl; 
		} else {
			denom= 1;
			for(int i= n-k-1; i < n; i++) 
				denom*= b_star_norm[i];

			N_next= ball_vol(k+1, Rvec[k]) * pvol_and_scale(Rvec, k/2+1) * fact(k/2+1) / denom;  
			N+= (N_prev+N_next)/2;
			cout << k-1 << " " << (N_prev+N_next)/2 << endl; 

			//V_act= ball_vol(k-1, Rvec[k-2]) * pvol_and_scale(Rvec, k/2)*fact(k/2) * 2 * Rvec[k]; // upper bound: V_{k-1}*2*R_k 
			//V_act= ball_vol(k, Rvec[k-1]) * ( pvol_and_scale(Rvec, k/2)*fact(k/2) +  pvol_and_scale(Rvec, k/2+1)*fact(k/2+1)) / 2; 

			cout << "# Simvol_" << k/2+1 << " = " <<  fact(k/2+1) << endl;
			cout << "# Polvol_" << k/2+1 << " = " <<  pvol_and_scale(Rvec, k/2+1) << endl;
			cout << "# V_ball(" << k+1 << ", " << Rvec[k] << ") = " << ball_vol(k+1, Rvec[k]) << endl;
			cout << "# |b^*_" << n-k-1 << "| = " << b_star_norm[n-k-1] << endl;
			cout << "# V_" << k+1 << " = " << ball_vol(k+1, Rvec[k]) * pvol_and_scale(Rvec, k/2+1) * fact(k/2+1) << endl;
			cout << "# denom_" << k << " = " << denom << endl;
//		}*/
		denom= 1;
		//for(int i= 0; i < k; i++) 
		for(int i= n-k; i < n; i++) 
			denom*= b_star_norm[i]; 

		int psd= (k+1)/2-1;
		V_act= ball_vol(k, Rvec[k-1]) * pvol_and_scale(Rvec, k) * fact(psd); 

		N+= V_act/denom/2;

		/*cout << "# Simvol_" << psd << " = " <<  fact(psd) << endl;
		cout << "# Polvol_" << psd << " = " <<  pvol_and_scale(Rvec, k) << endl;
		cout << "# V_ball(" << k << ", " << Rvec[k-1] << ") = " << ball_vol(k, Rvec[k-1]) << endl;
		cout << "# |b^*_" << n-k << "| = " << b_star_norm[n-k] << endl;
		cout << "# V_" << k << " = " << V_act << endl;
		cout << "# denom_" << k << " = " << denom << endl;
		cout << "# N_" << k << " = " << V_act/denom/2 << endl; */
		cout << k-1 << " " << V_act/denom/2 << endl; 
	}

	//N/= 2;
	cout << "# Predicted nodes: " << N << endl;
	return (t_reduc+t_node*N)/(pvol_and_scale(Rvec, n/2)*fact(n/2)); 
}

double t_extreme(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n) {
	int pdim= (n%2==0)?n/2:n/2+1;
	
	double* polytope_vols= new double[pdim+1];
	double* bounds= new double[pdim];

	for(int i= 0; i< n; i+=2)
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);

	double N= 0;
	double V_act, denom= 1;
	polytope_vols[pdim]= polytope_volume(polytope_vols, bounds, pdim);
	for(int i= 1; i<= n; i++) {
		if(polytope_vols[i/2] * fact(i/2) > 1)
			cout << " # " << polytope_vols[i/2] * fact(i/2) << endl;

		if(i%2==0)
			V_act= ball_vol(i, Rvec[i-1]) * polytope_vols[i/2] * fact(i/2); 
		else
			V_act= ball_vol(i, Rvec[i-1]) * (polytope_vols[i/2] + polytope_vols[i/2+1])*fact(i/2+1) / 2; 
			//V_act= ball_vol(i, Rvec[i-1]) * (polytope_vols[i/2] + polytope_vols[i/2+1])*fact(i/2+1) / 2; 

		denom*= b_star_norm[n-i];
		N+= V_act/denom;
		/*cout << "Simvol_" << i/2 << " = " <<  fact(i/2) << endl;
		cout << "Polvol_" << i/2 << " = " <<  polytope_vols[i/2] << endl;
		cout << "V_ball(" << i << ", " << Rvec[i-1] << ") = " << ball_vol(i, Rvec[i-1])<< endl;
		cout << "|b^*_" << n-i << "| = " << b_star_norm[n-i] << endl;
		cout << "V_" << i << " = " << V_act << endl;
		cout << "denom_" << i << " = " << denom << endl;
		cout << "N_" << i << " = " << V_act/denom << endl << endl;*/
	}
	
	N/= 2;

	cout << "Predicted nodes: " << N << endl;
	return (t_reduc+t_node*N)/(polytope_vols[pdim-1]*fact(pdim-1)); 
}

RR compute_p_succ(RR Rvec[], int n, int prec){
	int pdim= (n%2==0)?n/2:n/2+1;
	RR p_succ;
	p_succ.SetPrecision(prec);
	
	RR* polytope_vols= new RR[pdim+1];
	RR* bounds= new RR[pdim];

	for(int i= 0; i< n; i+=2) {
		bounds[i/2].SetPrecision(prec);	
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);
	}

	polytope_vols[pdim]= polytope_volume_RR(polytope_vols, bounds, pdim, prec);

	p_succ= polytope_vols[pdim-1]*fact_RR(pdim-1, prec);

	delete [] bounds;
	delete [] polytope_vols;

	return p_succ;
}

RR p_succ(RR Rvec[], int n, int prec){
	return compute_p_succ(Rvec, n, prec);
	}

static 
void generate_boundary_step(RR b_star_norm[], double t_node, double t_reduc, int n, RR* act, double delta, unsigned long iterations, RR& t_enum, int& changes) {
	RR* mod= new RR[n];
	
	for(int j= 0; j < n; j++)
		mod[j]= act[j]; 
	
	int change;
	int sign;
	RR time;
	changes= 0;
	for(unsigned long i= 0; i < iterations; i++) {
		bool flag= true;	
		while(flag) {
			change= 2*(rand()%(n/2));

			sign= rand()%2;
			if(sign!=1) sign--;
			
			flag= false;
			// Ensuring nonnegativity of the bounding function
			// (negative bounds don't make sense
			if((change==0) && (mod[change]+sign*delta<=0))
				flag= true;

			// Making sure that R doesn't changes neither in the odd nor in the even dimensional case
			// (Changing R would result either in p_succ= 0 or the algorithm returning with a vector longer than R)
			if((change==n-1) || (change==n-2)) 
				flag= true;
		} 

		mod[change]= mod[change+1]+= sign*delta;

		// Ensuring monotonity
		if((change < n-2) && mod[change] > mod[change+2])
			mod[change]= mod[change+1]= mod[change+2];
		if((change > 0) && mod[change] < mod[change-1])
			for(int j= change-1; j >= 0; j--) 
				if(mod[j] > mod[change])
					mod[j]= mod[change];
				else
					break;
					 

		time= t_extreme_RR(mod, b_star_norm, t_node, t_reduc, n, RR_PRECISION);

		if(time < t_enum) {
			t_enum= time;

			for(int j= 0; j < n; j++)
				act[j]= mod[j]; 
	
			changes++;
		} else
			for(int j= 0; j < n; j++)
				mod[j]= act[j]; 
			
	}

	delete [] mod;
}

void generate_boundary(RR b_star_norm[], double t_node, double t_reduc, int n, double Rvec[], double R, double delta, unsigned long iterations, double& p_succ, double& t_enum_d, bool quiet) {
	int changes;
	RR* act= new RR[n];
	
	RR_PI.SetPrecision(RR_PRECISION);
	RR_PI= ComputePi_RR();	

	RR t_enum;
	
	act[0].SetPrecision(RR_PRECISION);
	act[1].SetPrecision(RR_PRECISION);
	act[0]= act[1]= 2*R/n;
	for(int i= 2; i < n; i+=2) {
		act[i].SetPrecision(RR_PRECISION);
		act[i+1].SetPrecision(RR_PRECISION);
		act[i]= act[i+1]= act[i-1]+act[0];
	}
	act[n-1].SetPrecision(RR_PRECISION);
	act[n-1]= R;
	t_enum= t_extreme_RR(act, b_star_norm, t_node, t_reduc, n, RR_PRECISION);

	if(!quiet) {
		RR tmp;
		tmp= R;
		cout << "# Full enumeration time: " << t_full_RR(tmp, b_star_norm, t_node, t_reduc, n, RR_PRECISION) << endl;
		cout << "# Double step linear pruning: " << t_enum << endl;
	}

	srand(time(NULL));

	generate_boundary_step(b_star_norm, t_node, t_reduc, n, act, delta, iterations, t_enum, changes);
	
	for(int j= 0; j < n; j++)
		conv(Rvec[j], act[j]); 
			
	if(!quiet) 
		cout << "# Changes: " << changes << endl;

	conv(p_succ, compute_p_succ(act, n, RR_PRECISION));
	conv(t_enum_d, t_enum);

	delete [] act;
}
