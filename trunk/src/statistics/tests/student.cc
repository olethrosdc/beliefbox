/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// Test the student distribuition

#include "Student.h"

int main ()
{
  Student student(2);

  Matrix T(2,2);
  T(0,0)=2;
  T(1,1)=100;
  T(0,1)=-2;
  T(1,0)=-2;
  Vector x(2);
  student.setDegrees(10);
  student.setPrecision(T);
  for (real u=0; u<=10; u+=1.0) {
	x(0) = u;
	x(1) = u;
	printf("%f %f %f\n", u, student.log_pdf(x), student.pdf(x));
  }
  return 0;
}

