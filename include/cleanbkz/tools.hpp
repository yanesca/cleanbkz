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

#ifndef CLEANBKZ_TOOLS_H_
#define CLEANBKZ_TOOLS_H_

#include <NTL/RR.h>

using namespace NTL;

/** 	@file 
	@brief Contains a ceiling function that was somehow missing from my NTL library. */

/** Funcion for ceiling NTL::RR values */
void ceilPrec(RR& x, const RR& a, long p);

#endif
