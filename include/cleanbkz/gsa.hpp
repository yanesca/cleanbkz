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

#ifndef CLEANBKZ_GSA_H_
#define CLEANBKZ_GSA_H_


/** 	@file 
	@brief  Contains constant definitions required to compute vector lengths according to the Geometric Series Assumption. */

const double gsa_slope[]= { 	0, 0, 							// block sizes 0 and 1 make no sense
				-0.08570,  						// LLL
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  	// don't have parameters for block sizes 3-19
				-0.05322, -0.05292 ,-0.05259 ,-0.05218 ,-0.05171, 	// from BKZ-20 to BKZ-24
				-0.05121, -0.05071, -0.05022, -0.04976, -0.04934,	// from BKZ-25 to BKZ-29
				-0.04893, -0.04857, -0.04821, -0.04785, -0.04752,	// from BKZ-30 to BKZ-34
				-0.04724						// BKZ-35
			};

const double gsa_convexity[]= { 0, 0, 							// block sizes 0 and 1 make no sense
				2.06e-05,  						// LLL
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  	// don't have parameters for block sizes 3-19
				-2.71e-05, -3.11e-05, -3.74e-05, -4.58e-05, -5.45e-05,	// from BKZ-20 to BKZ-24
				-6.28e-05, -6.99e-05, -7.61e-05, -8.15e-05, -8.64e-05,	// from BKZ-25 to BKZ-29
				-9.11e-05, -9.54e-05, -9.94e-05, -1.034e-04, -1.068e-04,// from BKZ-30 to BKZ-34
				-1.105e-04						// BKZ-35
			};

#endif
