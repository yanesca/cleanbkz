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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <NTL/fileio.h>
#include <NTL/vec_double.h>
#include <NTL/new.h>
#include <iostream>
#include "enumeration.hpp"

using namespace std;


void enumerate(mat_ZZ& basis, double* prune, vec_RR& result) {
	
	BKZ_QP1(basis, 0.99, 30);		

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

	enumerate(mu, c, NULL, 0, mu1.NumRows()-1, mu1.NumRows(), result);
}


void enumerate(double** mu, double *c, double* prune, int jj, int kk, int m, vec_RR& result) {
	int s, t;
   	double cbar, t1, eta;

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

         long enum_cnt = 0;
   
         while (t <= kk) {
         
            ctilda[t] = ctilda[t+1] + 
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

            ForceToMem(&ctilda[t]);  // prevents an infinite loop
   
            if (prune != NULL) {
               eta = 0;
            }
            else
               eta = 0;
   
            if (ctilda[t] < cbar - eta) {
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

	result.SetLength(m);
	for(int i= 0; i < result.length(); i++) 
		conv(result[i],uvec[i]);

   delete [] ctilda;
   delete [] vvec;
   delete [] yvec;
   delete [] uvec;
   delete [] utildavec;
   delete [] Deltavec;
   delete [] deltavec;
} 

