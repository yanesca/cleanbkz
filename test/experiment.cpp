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

//TODO: specifikáció: bemenet egy dimenzió és lambda szorzói az egyes sublattice-ekhez

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#include <iostream>
#include <sstream>
#include <fstream>
#include <cleanbkz/boundary.hpp>
#include <cleanbkz/cjloss.hpp>
#include <cleanbkz/enumeration.hpp>
#include <cleanbkz/tools.hpp>
#include <NTL/LLL.h>
#include <NTL/new.h>

using namespace std;

double enumerate_spec(double** mu, double *c, int jj, int kk, int m); 

extern void predict_nodes_RR(RR  Rvec[], double b_star_norm[], int n); 

extern void predict_nodes(double Rvec[], double b_star_norm[], int n); 
extern double ball_vol(int k, double r);
extern void profile_enumerate_epr(double** mu, double *b, double* Rvec, int n, vec_RR& result); 

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim);

int main(int argc, char** argv) {
	int dim, bsize;
	stringstream ss;
	mat_ZZ basis;
	ZZ determinant;

// Processing commandline parameters (first is the dimension, the others the lambda multipliers)

	ss << argv[1];
	ss >> dim;
	ss.clear();

	double* muls= new double[dim];
	for(int i= 0; i< dim; i++)
		muls[i]= 0;

	for(int i= 2; i< argc; i++) {
		ss << argv[i];
		ss >> muls[i-2];
		ss.clear();
		}

// Initializing fixed parameters

	SetSeed(to_ZZ(0xa4d5330d));
	GenPrime(determinant,dim);
	bsize= 2;

// Creating random lattice
	
	gen_randlat(basis,determinant,dim); 
	BKZ_QP1(basis, 0.99, bsize); 

		cjloss l(dim, 0.94, 678);
		//basis= l.get_basis(bsize);

	mat_RR mu1;
	vec_RR c1;
	ComputeGS(basis,mu1,c1);

	double* c= new double[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], c1[i]);

	
	double** mu= new double*[mu1.NumRows()];
	for(int i= 0; i < mu1.NumRows(); i++)
		mu[i]= new double[mu1.NumCols()]; 
	for(int i= 0; i < mu1.NumRows(); i++)
		for(int j= 0; j < mu1.NumCols(); j++)
			conv(mu[i][j], mu1[i][j]);

	double* lprof= new double[dim];
	for(int i= 0; i< dim; i++)
		 lprof[i]= enumerate_spec(mu, c, dim-1-i, dim-1, dim); 


// Bounding function construction

	double* bound= new double[dim];
	for(int i= 0; i < dim; i++) 
		if(muls[i]==0)
			bound[i]= lprof[dim-1];
		else
			bound[i]= lprof[i]*muls[i];

	for(int i= 0; i < dim; i++) 
		cout << bound[i] << endl;
	

/*/ Enumeration
	vec_RR solution;	
	profile_enumerate_epr(mu, c, act, dim, solution); 

// Prediction	

	for(int i= 0; i < mu1.NumRows(); i++) 
		conv(c[i], sqrt(c1[i]));

	for(int i= 0; i < dim; i++)  
		manual[i]= sqrt(manual[i]);
		

	//predict_nodes(act, c, dim);
	RR* act_RR= new RR[dim];
	for(int i= 0; i < dim; i++)
		act_RR[i]= act[i];
	predict_nodes_RR(act_RR, c, dim);
*/

}

void gen_randlat(mat_ZZ& basis, ZZ determinant, int dim) {
	basis.SetDims(dim, dim);
	
	ZZ s;
	s= time(NULL);
	SetSeed(s);
	for(int i= 0; i < dim; i++) {
		basis[i][0]= RandomBnd(determinant);
		basis[i][i]= 1;
	}

	basis[0][0]= determinant;
} 

double enumerate_spec(double** mu, double *c, int jj, int kk, int m) {
	int s, t;
   	double cbar;
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

         cbar = c[jj];
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

         while (t <= kk) {
       
            ctilda[t] = ctilda[t+1] + 
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];


            ForceToMem(&ctilda[t]);  // prevents an infinite loop
   
            if (ctilda[t] < cbar) {
               if (t > jj) {
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
               }
               else {
                  cbar = ctilda[jj];
                  for (int i = jj; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec[t] = -Deltavec[t];
               if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         } 

	/*vec_RR result;
	result.SetLength(m);
	for(int i= 0; i < result.length(); i++) 
		conv(result[i],uvec[i]);*/
	//cout << result << endl;

   delete [] ctilda;
   delete [] vvec;
   delete [] yvec;
   delete [] uvec;
   delete [] utildavec;
   delete [] Deltavec;
   delete [] deltavec;

	return cbar;
} 

