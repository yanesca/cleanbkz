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

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <NTL/RR.h>

using namespace NTL;

#define RR_PRECISION 120

/** 	@file 
	@brief Contains algorithms for finding the optimal bounding function for a lattice. */

/**  Estimates the running time of pruned enumeration. It performs an exact computation based on the Dirichlet distribution as briefly explained in: Nicolas Gama, Phong Q. Nguyen and Oded Regev "Lattice Enumeration using Extreme Pruning", Advances in Cryptology – EUROCRYPT 2010 Lecture Notes in Computer Science Volume 6110, 2010, pp 257-278. */ 
RR t_extreme(
	RR R[],		 	 //!< The boundary function used for pruning. The odd entries have to be equal to the corresponding even entries (For example.: 1,1,2,2,5,5,... )  
	RR b_star_norm[],	 //!< The lengths of the vectors in the Gram-Schmidt basis
	double t_node,		 //!< The time the enumeration spends with processing a single node 
	double t_reduc,		 //!< The running time of the basis reduction algorithm used 
	int n,			 //!< The dimension of the lattice
	long prec		 //!< The floating point precision of the computations
	);

/**  Estimates the success probability of pruned enumeration. It performs an exact computation based on the Dirichlet distribution as briefly explained in: Nicolas Gama, Phong Q. Nguyen and Oded Regev "Lattice Enumeration using Extreme Pruning", Advances in Cryptology – EUROCRYPT 2010 Lecture Notes in Computer Science Volume 6110, 2010, pp 257-278. */ 
RR p_succ(
	RR Rvec[],		 //!< The boundary function used for pruning. The odd entries have to be equal to the corresponding even entries (For example.: 1,1,2,2,5,5,... ) 
	int n, 			 //!< The dimension of the lattice
	int prec		 //!< The floating point precision of the computations
	);

/**  Generates an optimal boundary function for enumeration with extreme pruning. Uses numerical optimization with random changes as briefly explained in: Nicolas Gama, Phong Q. Nguyen and Oded Regev "Lattice Enumeration using Extreme Pruning", Advances in Cryptology – EUROCRYPT 2010 Lecture Notes in Computer Science Volume 6110, 2010, pp 257-278. */ 
void generate_boundary(	
	NTL::RR b_star_norm[],	//!< The lengths of the vectors in the Gram-Schmidt basis 
	double t_node,		//!< The time the enumeration spends with processing a single node 
	double t_reduc, 	//!< The running time of the basis reduction algorithm used
	int n, 			//!< The dimension of the lattice
	double Rvec[], 		//!< Container with enough space reserved for the result (n double values) 
	NTL::ZZ R,		//!< The square length of the vector the enumeration is looking for 
	ZZ delta,		//!< The step of the random modifications  
	unsigned long iterations,	//!< The number of random modifications to test
	double& p_succ,		//!< Container for the success probability of the resulting function
	double& t_enum,		//!< Container for the estimated running time of the resulting function
	bool quiet		//!< Set this false if you do not want this function to print on the standard output
	); 

/**  Continues a previously started boundary function optimization for enumeration with extreme pruning. Uses numerical optimization with random changes as briefly explained in: Nicolas Gama, Phong Q. Nguyen and Oded Regev "Lattice Enumeration using Extreme Pruning", Advances in Cryptology – EUROCRYPT 2010 Lecture Notes in Computer Science Volume 6110, 2010, pp 257-278. */ 
void continue_boundary_gen(
	std::ifstream& infile,	//!< Reference to the input file stream containing the information about the previous computation 
	unsigned long iterations, 	//!< The number of random modifications to test
	double** Rvec,		//!< This parameter will hold the result 
	int& dim		//!< This parameter will hold the length of the result 
	); 

#endif
