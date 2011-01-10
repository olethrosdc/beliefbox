/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef READ_FILE_H
#define READ_FILE_H

#include <vector>
#include "Matrix.h"

int FileToIntVector(std::vector<int>& data, const char* fname, int tmpT);
int ReadClassData(Matrix& data, std::vector<int>& labels, const char* fname);
int ReadFloatDataASCII(Matrix& data, const char* fname);

#endif
