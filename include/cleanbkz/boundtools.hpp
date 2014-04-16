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

#ifndef BOUNDTOOLS_H_
#define BOUNDTOOLS_H_

/** @file Contains tools for constructing boundary functions for pruned enumeration. */

/**  Measures the processing time of a node and the basis reduction enumerate_epr function. Generates random lattices, prerocesses the bases with the prescribed blocksize and terminates the enumeration after processig \p nodes nodes. Computes the sample means of the runing times after /p samples experiments. */ 

void measure_epr(
	int dimension, 			//!< The dimension of the random input lattices 
	int beta,			//!< The blocksize of BKZ used to preprocess the basis
	unsigned long samples,		//!< The number of samples to consider			
	unsigned long nodes,		//!< The number of nodes to process with the enumeration			
	double& t_node,			//!< Container for the sample mean of the processing time of a node 
	double& t_reduc,		//!< Container for the sample mean of the basis reduction 
	unsigned long& zeroes,		//!< Container for the number of samples with zero t_node 
	unsigned long& bias		//!< Container for the number of samples with fewer nodes processed
	); 

#endif
