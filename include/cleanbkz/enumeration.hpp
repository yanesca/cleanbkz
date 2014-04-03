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

/** @file */
/** Searches for shortest vector with pruned enumeration. The inputs are the lattice basis and the pruning function and a container for the result (short vector int the lattice defined by the basis). */
void enumerate(mat_ZZ& basis, double* pruning, vec_RR& result); 	

#endif
