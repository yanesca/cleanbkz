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
#include <fstream>
#include <NTL/RR.h>
#include <cleanbkz/boundary.hpp>

using namespace std;

inline double fact(double x) {
	return (x < 2 ? 1 : x * fact(x - 1));
}

RR fact_RR(int x) {
	RR ret;
	ret.SetPrecision(RR_PRECISION);
	ret= 1;

	for(int i= 2; i <= x; i++) 
		ret*= i;

	return ret;
	}


double ball_vol(int k, double r) {
	if (k%2==0) 
		return pow(M_PI, k/2) / fact(k/2) * pow(r, k);
	else
		return 2 * pow(4*M_PI, k/2) * fact(k/2) / fact(k) * pow(r, k);
}

RR RR_PI;

RR ball_vol_RR(int k, RR r) {
	RR pwr, exp, pipow;
	pwr.SetPrecision(RR_PRECISION);
	exp.SetPrecision(RR_PRECISION);
	pipow.SetPrecision(RR_PRECISION);
	
	exp= k;
	pow(pwr, r, exp);
 
	if (k%2==0) { 
		exp= k/2;
		pow(pipow, RR_PI, exp); 
		return pipow / fact_RR(k/2) * pwr;
	} else {
		exp= k/2;
		pow(pipow, 4*RR_PI, exp);
		return 2 * pipow * fact_RR(k/2) / fact_RR(k) * pwr;
	}
}

RR integral_even_RR(int h, int l, RR tvec[], RR vvec[]) {
	RR tmp, pwr, ret;

	tmp.SetPrecision(RR_PRECISION);
	pwr.SetPrecision(RR_PRECISION);
	ret.SetPrecision(RR_PRECISION);

	ret= 1;
	if(h==0) 
		return ret;

	ret= 0;

	vvec[h-1]= integral_even_RR(h-1, l, tvec, vvec);

	for(int i= 0; i < h; i++) {
		tmp= h-i;
		pow(pwr, tvec[l-h], tmp); 
		ret+= pwr*vvec[i]/fact_RR(h-i)*((h-i-1)%2==0?1:-1);
	}
		
	return ret;
}


double integral_even(int h, int l, double tvec[], double vvec[]) {
	double ret= 0;

	if(h==0) 
		return 1;

	vvec[h-1]= integral_even(h-1, l, tvec, vvec);

	for(int i= 0; i < h; i++)
		ret+= pow(tvec[l-h],h-i)*vvec[i]/fact(h-i)*((h-i-1)%2==0?1:-1);
		
	return ret;
}

RR integral_odd_RR(int h, int l, RR tvec[], RR vvec[]) {
	RR tmp, pwr, ret, pow2, two;

	tmp.SetPrecision(RR_PRECISION);
	pwr.SetPrecision(RR_PRECISION);
	pow2.SetPrecision(RR_PRECISION);
	ret.SetPrecision(RR_PRECISION);
	two.SetPrecision(RR_PRECISION);

	two= 2;
	ret= 1;
	if(h==0) 
		return ret;

	ret= 0;

	vvec[h-1]= integral_odd_RR(h-1, l, tvec, vvec);

	tmp= (2*h+1)/2.0;
	pow(pwr, 1-tvec[l-h], tmp); 
	tmp= 2*h+1;
	pow(pow2, two, tmp); 

	ret-= (pwr*pow2/fact_RR(2*h+2))*fact_RR(h+1);
	for(int i= 1; i < h; i++) {
		tmp= h-i;
		pow(pwr, tvec[l-h], tmp); 
		ret+= pwr*vvec[i]/fact_RR(h-i)*((h-i-1)%2==0?1:-1);
	}
		
	if(h!=l)		
		return ret;

	return ret + pow2/fact_RR(2*h+2)*fact_RR(h+1);
}


double integral_odd(int h, int l, double tvec[], double vvec[]) {
	double ret= 0;

	if(h==0) 
		return 1;

	vvec[h-1]= integral_odd(h-1, l, tvec, vvec);

	ret-= pow(1-tvec[l-h], (2*h+1)/2.0)*pow(2,2*h+1)*fact(h+1)/fact(2*h+2);
	for(int i= 1; i < h; i++)
		ret+= pow(tvec[l-h],h-i)*vvec[i]/fact(h-i)*((h-i-1)%2==0?1:-1);

	if(h!=l)		
		return ret;

	return ret + pow(2,2*h+1)*fact(h+1)/fact(2*h+2);
}

RR t_full_RR(RR R, RR b_star_norm[], double t_node, double t_reduc, int n, int prec) {
	RR N, V_act, denom;
	N.SetPrecision(prec);
	V_act.SetPrecision(prec);
	denom.SetPrecision(prec);
	denom= 1;
	for(int i= 1; i<= n; i++) {
		V_act= ball_vol_RR(i, R); 

		denom*= b_star_norm[n-i];
		N+= V_act/denom;
	}

	N/= 2;
	return t_reduc+t_node*N; 
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

RR ci_prob_RR(RR Rvec[], int k) {
	RR ret;
	int l= k/2;
	RR* bounds= new RR[l];
	RR* vvec= new RR[l];

	for(int i= 0; i< l; i++)
		bounds[i]= Rvec[2*i]*Rvec[2*i]/(Rvec[k-1]*Rvec[k-1]);

	if(k%2==0)
		ret= integral_even_RR(l, l, bounds, vvec)*fact_RR(l);
	else
		ret= integral_odd_RR(l, l, bounds, vvec)*fact_RR(2*l+2)/(power2_RR(2*l+1)*fact_RR(l+1));

	delete [] bounds;
	delete [] vvec;

	return ret;
}

double ci_prob(double Rvec[], int k) {
	//TODO: ez a jó, a másikat is igy kijavitani
	double ret;
	int l= k/2;
	double* bounds= new double[l];
	double* vvec= new double[l];

	/*cout << "# khalf: " << l << endl;
	cout << "# Rvec: [ ";
	for(int i= 0; i< 2*l; i++)
		cout << Rvec[i] << " ";	
	cout << "] "<< endl;*/

	for(int i= 0; i< l; i++)
		bounds[i]= Rvec[2*i]*Rvec[2*i]/(Rvec[k-1]*Rvec[k-1]);

	/*cout << "# bounds: [ ";
	for(int i= 0; i< l; i++)
		cout << bounds[i] << " ";	
	cout << "] "<< endl << "# " << endl;*/

	if(k%2==0)
		ret= integral_even(l, l, bounds, vvec)*fact(l);
	else
		ret= integral_odd(l, l, bounds, vvec)*fact(2*l+2)/(pow(2,2*l+1)*fact(l+1));

	delete [] bounds;
	delete [] vvec;

	return ret;
}

RR p_succ(RR Rvec[], int n){
	RR ret;
	int l= (n%2==1)?n/2:n/2-1;
	RR* bounds= new RR[l];
	RR* vvec= new RR[l];

	for(int i= 0; i< l; i++)
		bounds[i]= Rvec[2*i]*Rvec[2*i]/(Rvec[n-1]*Rvec[n-1]);

	if(n%2==0)
		ret= integral_even_RR(l, l, bounds, vvec)*fact_RR(l);
	else
		ret= integral_odd_RR(l, l, bounds, vvec)*fact_RR(2*l+2)/(power2_RR(2*l+1)*fact_RR(l+1));

	delete [] bounds;
	delete [] vvec;

	return ret;
}

RR t_extreme_RR(RR Rvec[], RR b_star_norm[], double t_node, double t_reduc, int n) {
	RR N;
	RR V_act; 	
	RR denom;

	N= 0;
	for(int k= 1; k<= n; k++) {
		denom= 1;

		for(int i= n-k; i < n; i++) 
			denom*= b_star_norm[i]; 

		V_act= ball_vol_RR(k, Rvec[k-1]) * ci_prob_RR(Rvec, k); 

		N+= V_act/denom;
	}
	

	N/= 2;
	return (t_reduc+t_node*N)/p_succ(Rvec, n); 
}

/*double t_extreme(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n) {
	double N= 0;
	double V_act; 	
	double denom;

	for(int k= 1; k<= n; k++) {
		denom= 1;

		for(int i= n-k; i < n; i++) 
			denom*= b_star_norm[i]; 

		V_act= ball_vol(k, Rvec[k-1]) * ci_prob(Rvec, k); 

		N+= V_act/denom;
	}
	

	N/= 2;
	return (t_reduc+t_node*N)/ci_prob(Rvec, n/2); 
}*/

void predict_nodes(double Rvec[], double b_star_norm[], int n) {
	double N= 0;
	double V_act; 	
	double denom;

	cout << "# T_extreme: " << endl << "# GS-lengths: [ ";
	for(int i= 0; i < n; i++)
		cout << b_star_norm[i] << " ";	
	cout << "]" << endl << "# Bounds: [ ";
	for(int i= 0; i < n; i++)
		cout << Rvec[i] << " ";	
	cout << "]" << endl << "#" << endl;

	cout << "\"Predicted\"" << endl << endl;

	ifstream myfile ("ratio.tmp");
	unsigned long nodes;
	unsigned long sum= 0;

	for(int k= 1; k<= n; k++) {
		denom= 1;
		for(int i= n-k; i < n; i++) 
			denom*= b_star_norm[i]; 

		V_act= ball_vol(k, Rvec[k-1]) * ci_prob(Rvec, k); 

		N+= V_act/denom/2;

		cout << "# N_" << k << " = " << V_act/denom/2 << endl; 
		myfile >> nodes; sum+= nodes; 
		cout << "# Measured nodes: " << nodes << endl;
		cout << "# Difference: " << nodes - V_act/denom/2 << endl;		
		cout << "# Ratio: " << nodes/(V_act/denom/2) << endl;
		cout << "# GH_" << k << " = " << pow(denom/ball_vol(k, 1),1.0/k) << endl;
		cout << "# R/GH in dim " << k << " = " << Rvec[k-1]/ pow(denom/ball_vol(k, 1),1.0/k) << endl;

		cout << k-1 << " " << V_act/denom/2 << endl; 
	}
	
	myfile.close();

	cout << "# Predicted nodes: " << N << endl;
	cout << "# Measured nodes: " << sum << endl;
}

void predict_nodes_RR(RR Rvec[], double b_star_norm[], int n) {
	RR N;
	N= 0;
	RR V_act; 	
	RR denom;

	RR_PI.SetPrecision(RR_PRECISION);
	RR_PI= ComputePi_RR();	

	cout << "# T_extreme: " << endl << "# GS-lengths: [ ";
	for(int i= 0; i < n; i++)
		cout << b_star_norm[i] << " ";	
	cout << "]" << endl << "# Bounds: [ ";
	for(int i= 0; i < n; i++)
		cout << Rvec[i] << " ";	
	cout << "]" << endl << "#" << endl;

	cout << "\"Predicted\"" << endl << endl;

	ifstream myfile ("ratio.tmp");
	unsigned long nodes;
	unsigned long sum= 0;

	for(int k= 1; k<= n; k++) {
		denom= 1;
		for(int i= n-k; i < n; i++) 
			denom*= b_star_norm[i]; 

		V_act= ball_vol_RR(k, Rvec[k-1]) * ci_prob_RR(Rvec, k); 

		N+= V_act/denom/2;

		cout << "# Ball_" << k << " = " << ball_vol_RR(k, Rvec[k-1]) << endl;
		cout << "# P_" << k << " = " << ci_prob_RR(Rvec, k) << endl;
		cout << "# N_" << k << " = " << V_act/denom/2 << endl; 
		myfile >> nodes; sum+= nodes; 
		cout << "# Measured nodes: " << nodes << endl;
		cout << "# Difference: " << nodes - V_act/denom/2 << endl;		
		//cout << "# Ratio: " << nodes/(V_act/denom/2) << endl;
		//cout << "# GH_" << k << " = " << pow(denom/ball_vol(k, 1),1.0/k) << endl;
		//cout << "# R/GH in dim " << k << " = " << Rvec[k-1]/ pow(denom/ball_vol(k, 1),1.0/k) << endl;

		cout << k-1 << " " << V_act/denom/2 << endl; 
	}
	
	myfile.close();

	cout << "# Predicted nodes: " << N << endl;
	cout << "# Measured nodes: " << sum << endl; 
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
					 

		time= t_extreme_RR(mod, b_star_norm, t_node, t_reduc, n);

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

void generate_boundary(RR b_star_norm[], double t_node, double t_reduc, int n, double Rvec[], double R, double delta, unsigned long iterations, double& p_succ_v, double& t_enum_d, bool quiet) {
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
	t_enum= t_extreme_RR(act, b_star_norm, t_node, t_reduc, n);

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

	conv(p_succ_v, p_succ(act, n));
	conv(t_enum_d, t_enum);

	delete [] act;
}
