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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef ENUMERATION_H_
#define ENUMERATION_H_

#include<NTL/LLL.h>

using namespace NTL;

/** @file Contains enumeration algorithms for finding the shortest vector in a lattice. */

/** Searches for the shortest vector with pruned enumeration. It performs an enumeration and returns with the SHORTEST vector found in the domain determined by the boundary function. It is based on the NTL's enumeration algorithm, extended with extreme pruning. (This i probably the algorithm published in this paper:  C. P. Schnorr and M. Euchner "Lattice basis reduction: Improved practical algorithms and solving subset sum problems", Mathematical Programming, 2 August 1994, Volume 66, Issue 1-3, pp 181-199. */
void enumerate_ntl(
	mat_ZZ& basis, 		//!< The lattice basis 
	int beta,		//!< The blocksize of BKZ used to preprocess the basis
	double* prune,		//!< The bounding function for pruning  (IT IS NOT IMPLEMENTED YET)
	vec_RR& result		//!< Container for the return value
	); 	

/**  Searches for a short vector with pruned enumeration. It performs an enumeration and returns with the FIRST vector found in the domain determined by the boundary function. It is the enumeration algorithm published in this paper: Nicolas Gama, Phong Q. Nguyen and Oded Regev "Lattice Enumeration using Extreme Pruning", Advances in Cryptology – EUROCRYPT 2010 Lecture Notes in Computer Science Volume 6110, 2010, pp 257-278. */ 
void enumerate_epr(
	mat_ZZ& basis, 		//!< The lattice basis 
	int beta,		//!< The blocksize of BKZ used to preprocess the basis
	double* prune,		//!< The bounding function for pruning 
	vec_RR& result		//!< Container for the return value
	); 

#endif
