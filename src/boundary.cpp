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

using namespace std;

inline double fact(double x) {
	return (x < 2 ? 1 : x * fact(x - 1));
}

// TODO: 3D-ben tesztelni
// TODO: magasabb dimenziós teszteket irni hozzá
//inline 
double ball_vol(int k, double r) {
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
		
	if(ret<0) {
		cout << " # pol_vol in dimension " << dim << " is: " << ret << endl << " # with bounds: ";
		for(int i= 0; i < dim; i++)
			cout << bounds[i] << " ";
		cout << endl;
	}

	return ret;
}

double t_extreme_reference(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n) {
	int pdim= (n%2==0)?n/2:n/2+1;
	
	double* polytope_vols= new double[pdim+1];
	double* bounds= new double[pdim];

	for(int i= 0; i< n; i+=2)
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);

	double N= 0;
	double V_act, denom;
	polytope_vols[pdim]= polytope_volume(polytope_vols, bounds, pdim);
	for(int i= 1; i<= n; i++) {
		if(i%2==0)
			V_act= ball_vol(i, Rvec[i-1]) * polytope_vols[i/2+1]; 
		else
			V_act= ball_vol(i, Rvec[i-1]) * (polytope_vols[i/2]+polytope_vols[i/2+1]) / 2; 

		denom= 1;
		for(int j= n-i; j < n; j++) 
			denom*= b_star_norm[j];

		N+= V_act/denom;
	}

	return (t_reduc+t_node*N)/polytope_vols[pdim]/2; 
}

// TODO: megirni a naiv változatot és tesztelni
// TODO: megirni a t_full -t és összehasonlitani 
double t_extreme(double Rvec[], double b_star_norm[], double t_node, double t_reduc, int n) {
	int pdim= (n%2==0)?n/2:n/2+1;
	
	double* polytope_vols= new double[pdim+1];
	double* bounds= new double[pdim];

	for(int i= 0; i< n; i+=2)
		bounds[i/2]= Rvec[i]*Rvec[i]/(Rvec[n-1]*Rvec[n-1]);

	double N= 0;
	polytope_vols[pdim]= polytope_volume(polytope_vols, bounds, pdim);
	for(int i= 1; i<= n; i++) {
		if(i%2==0)
			N+= ball_vol(i, Rvec[i-1]) * polytope_vols[i/2+1]; 
		else
			N+= ball_vol(i, Rvec[i-1]) * (polytope_vols[i/2]+polytope_vols[i/2+1]) / 2; 

		N/= b_star_norm[i-1];
		//cout << "# norm of b_star_" << i << ": " << b_star_norm[i-1] << endl;
		//cout << "# vol_ball_" << i << ": " << ball_vol(i, R[i-1]) << endl;
		/*if(ball_vol(i, Rvec[i-1])<0)
			cout << "# ball_vol(" << i << ", " <<  Rvec[i-1] << ") = " << ball_vol(i, Rvec[i-1]) << endl;*/
	}

	//cout << "# N: " << N << endl;
	//cout << "# red/node: " << t_reduc << " / " << t_node << endl;
	//cout << "# p_succ: " << polytope_vols[pdim] << endl;
	/*cout << " # Bounds: ";
	for(int i= 0; i< pdim; i++) 
		cout << bounds[i] << " ";
	cout << endl;*/
	//cout << "# t_extreme: " << (t_reduc+t_node*N)/polytope_vols[pdim]/2 << endl;
	return (t_reduc+t_node*N)/polytope_vols[pdim]/2; 
}

void generate_boundary(double b_star_norm[], double t_node, double t_reduc, int n, double Rvec[], double R, double delta, unsigned long iterations, double& p_succ, double& t_enum) {
	//TODO: visszaadni a végső p_succ értéket
	/* Itt a deltától függően kellene kerekiteni a dolgokat, mert egyáltalán nem biztos, hogy elég kicsi lesz a delta ahhoz, hogy szép sima legyen a függvény */
	Rvec[0]= Rvec[1]= 2*R/n;
	for(int i= 2; i < n-1; i+=2)
		Rvec[i]= Rvec[i+1]= Rvec[i-1]+Rvec[0];
	Rvec[n-1]= R;
	t_enum= t_extreme(Rvec, b_star_norm, t_node, t_reduc, n);

	srand(time(NULL));

	int change;
	int sign;
	double time;
	for(unsigned long i= 0; i < iterations; i++) {
	
		bool flag= true;	
		while(flag) {
			change= 2*(rand()%(n/2));

			sign= rand()%2;
			if(sign!=1) sign--;
			
			flag= false;
			if((change==0) && (Rvec[change]+sign*delta<=0))
				flag= true;
			if((change!=0) && (Rvec[change]+sign*delta<Rvec[change-1])) 
				flag= true;
			if((change==n-1) && (Rvec[change]+sign*delta>R))
				flag= true;
			if((change!=n-1) && (Rvec[change]+sign*delta>Rvec[change+1]))
				flag= true;
		} 

		Rvec[change]= Rvec[change+1]+= sign*delta;

		time= t_extreme(Rvec, b_star_norm, t_node, t_reduc, n);
		//cout << "# t_" << i << ": " << time << endl;

		if(time < t_enum)
			t_enum= time;
		else
			Rvec[change]= Rvec[change+1]-= sign*delta;
	}

	p_succ= -1;
				
}
