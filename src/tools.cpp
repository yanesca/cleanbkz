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

#include <string>
#include <algorithm>
#include <cleanbkz/tools.hpp>

using namespace std;

char* get_cmd_option(char** begin, char** end, const string& option) {
    char** itr= find(begin, end, option);

    if (itr != end && ++itr != end) {
        return *itr;
    }

    return 0;
}

bool cmd_option_exists(char** begin, char** end, const string& option) {
    return find(begin, end, option) != end;
}

void ceilPrec(RR& x, const RR& a, long p) {
	if (p < 1 || NTL_OVERFLOW(p, 1, 0))
		Error("CeilPrec: bad precision");

 	ConvPrec(x,a,p);

	cout << "x = " << x << endl;
	
	RR epsilon;
	epsilon= 10;
	RR exp;
	exp= (x.exponent()+15);
	cout << "exponent = " << x.exponent() << endl;
	cout << "precision = " << x.precision() << endl;
	epsilon= pow(epsilon, exp); 
	cout << "epsilon = " << epsilon << endl;

	x+= epsilon;
}

