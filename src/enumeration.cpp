/*
    Copyright (C) 2007 Victor Shoup. 
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

#include <NTL/fileio.h>
#include <NTL/vec_double.h>
#include <NTL/new.h>
#include <iostream>
#include <cleanbkz/enumeration.hpp>

using namespace std;

static void enumerate_ntl(double** mu, double *c, double* prune, int jj, int kk, int m, vec_RR& result) {
	int s, t;
   	//double cbar, eta;
	double t1;

   double *ctilda;
   ctilda = NTL_NEW_OP double[m+2];
   if (!ctilda) Error("ENUMERATE: out of memory");

   double *vvec;
   vvec = NTL_NEW_OP double[m+2];
   if (!vvec) Error("ENUMERATE: out of memory");

   double *yvec;
   yvec = NTL_NEW_OP double[m+2];
   if (!yvec) Error("ENUMERATE: out of memory");

   double *uvec;
   uvec = NTL_NEW_OP double[m+2];
   if (!uvec) Error("ENUMERATE: out of memory");

   double *utildavec;
   utildavec = NTL_NEW_OP double[m+2];
   if (!utildavec) Error("ENUMERATE: out of memory");

   long *Deltavec;
   Deltavec = NTL_NEW_OP long[m+2];
   if (!Deltavec) Error("ENUMERATE: out of memory");

   long *deltavec;
   deltavec = NTL_NEW_OP long[m+2];
   if (!deltavec) Error("ENUMERATE: out of memory");

         //cbar = c[jj];
         utildavec[jj] = uvec[jj] = 1;
   
         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;
   
         s = t = jj;
         deltavec[jj] = 1;
   
         for (int i = jj+1; i <= kk+1; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }

	int nodes= 0;
         while (t <= kk) {
       
		nodes++; 
            ctilda[t] = ctilda[t+1] + 
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

		//cout << nodes << "\tt: " << t << "\tctilda: " << ctilda[t] << "\tprune (kk-t-1): " << prune[kk-t-2] << endl; 

            ForceToMem(&ctilda[t]);  // prevents an infinite loop
   
            if ((t>jj) && (ctilda[t] <= prune[kk-t])) {
            // if (ctilda[t] <= prune[kk-t]) {
               //if (t > jj) {
                  t--;
                  t1 = 0;
                  for (int i = t+1; i <= s; i++)
                     t1 += utildavec[i]*mu[i][t];
                  yvec[t] = t1;
                  t1 = -t1;
                  if (t1 >= 0)
                     t1 = ceil(t1-0.5);
                  else
                     t1 = floor(t1+0.5);
                  utildavec[t] = vvec[t] = t1;
                  Deltavec[t] = 0;
                  if (utildavec[t] > -yvec[t]) 
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
               /*}
               else {
                  //cbar = ctilda[jj];
		  break;
                  for (int i = jj; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }*/
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec[t] = -Deltavec[t];
               if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         } 

	result.SetLength(m);
	for(int i= 0; i < result.length(); i++) 
		conv(result[i],uvec[i]);

	cout << "# NTL nodes: " << nodes << endl;

   delete [] ctilda;
   delete [] vvec;
   delete [] yvec;
   delete [] uvec;
   delete [] utildavec;
   delete [] Deltavec;
   delete [] deltavec;
} 

void enumerate_ntl(mat_ZZ& basis, int beta, double* prune, vec_RR& result) {
	
	BKZ_QP1(basis, 0.99, beta);		

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(basis,mu1,c1);

	double** mu= new double*[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		mu[i]= new double[mu1.NumCols()]; 
	for(int i= 0; i < mu1.NumRows(); i++)
		for(int j= 0; j < mu1.NumCols(); j++)
			conv(mu[i][j], mu1[i][j]);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
			conv(c[i], c1[i]);

	enumerate_ntl(mu, c, prune, 0, mu1.NumRows()-1, mu1.NumRows(), result);
}

void enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result, unsigned long &termination, double &time) {

	cout << "# Enumerate: " << endl << "# GS-squares: [ ";
	for(int i= 0; i < n; i++)
		cout << b[i] << " ";	
	cout << "]" << endl << "# Bound squares: [ ";
	for(int i= 0; i < n; i++)
		cout << Rvec[i] << " ";	
	cout << "]" << endl << "#" << endl;

	bool pruning= true; 
	if(Rvec==NULL) {
		pruning= false;
		Rvec= new double[n];
		for(int i= 0; i < n; i++)
			Rvec[i]= b[0];
		}

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

	unsigned long nodes= 0;
	unsigned long* prof= new unsigned long[n];
	for(int i= 0; i < n; i++)
		prof[i]= 0;

	int k= 0;
	clock_t end,begin= clock();	

	cout << "\"Measured\"" << endl << endl;

	while((termination==0)||(nodes < termination)) {
		rhovec[k]= rhovec[k+1]+(vvec[k]-cvec[k])*(vvec[k]-cvec[k])*b[k];	

		nodes++;
		prof[n-k-1]++;
		//cout << nodes << "\tk: " << k << "\trhovec: " << rhovec[k] << "\tRvec (n-k-1): " << Rvec[n-k-1] << endl; 
		/*cout << "vvec: " << endl;
		for(int i= 0; i < n; i++)
			cout << vvec[i] << " ";
		cout << endl;*/

		//if(rhovec[k] <= Rvec[n-k-1]) {
		if((k!=0) && (rhovec[k] <= Rvec[n-k-1])) {
			/*if(k==0) 
				break;
			else {*/
				k--;
				rvec[k]= rvec[k]>rvec[k+1]?rvec[k]:rvec[k+1];		
				for(int i= rvec[k+1]; i>k; i--)
					sigmamat[i][k]= sigmamat[i+1][k] + vvec[i]*mu[i][k];
				cvec[k]= -sigmamat[k+1][k];

				vvec[k]= lround(cvec[k]);
				wvec[k]= 1;

			//}
		} else {
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

	cout << "# Nodes: " << nodes << endl;

	end= clock();
	time= (double) (end-begin) / CLOCKS_PER_SEC / nodes;
	//if(termination!=0)
		termination-= nodes;

	for(int i= 0; i < n; i++)
		cout << i << " " << prof[i] << endl;

	if(k==0) {
		result.SetLength(n);
		for(int i= 0; i < result.length(); i++) 
			conv(result[i],vvec[i]);
	} else
		result.SetLength(0);

	if(!pruning)
		delete [] Rvec;
	for(int i= 0; i < n+1; i++)
		delete [] sigmamat[i]; 
	delete [] sigmamat;
	delete [] rvec; 
	delete [] vvec;
	delete [] wvec;
	delete [] cvec;
	delete [] rhovec;
}

void enumerate_epr(mat_ZZ& basis, int beta, double* prune, vec_RR& result) {
	double time;

	
//	cout << "Original: " << endl << basis << endl;

	BKZ_QP1(basis, 0.99, beta);		
	
//	cout << "Rank: " << rank << endl;
//	cout << "Reduced: " << endl << basis << endl;

	//cout << "Basis: " << endl << basis;

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(basis,mu1,c1);

//	cout << "Reduced: " << endl << mu1 << endl;

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

	unsigned long termination= 0;
	enumerate_epr(mu, c, prune, mu1.NumRows(), result, termination, time);
	
	cout << "# Nodes processed: " << -termination << endl;
}
