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

/** 	@file 
	@brief Contains a functions for parsing the command line, used by the applications. */

#include <string>

using namespace std;

//TODO: documentation

char* get_cmd_option(char** begin, char** end, const string& option);

bool cmd_option_exists(char** begin, char** end, const string& option);

#endif
