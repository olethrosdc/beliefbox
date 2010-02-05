/* -*- Mode: c++ -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifdef MAKE_MAIN

#include "Ring.h"
#include "real.h"
#include "Random.h"
#include <iostream>
#include <list>
using namespace std;

bool TestRing (int number_of_entries, int ring_size);

int main(void)
{
    int n_tests = 10;
    int errors = 0;

    cout << "\nRunning ordered fixed list test\n";
    cout << "n_tests: " << n_tests << endl;

    for (int test = 0; test<n_tests; ++test) {
        bool res = TestRing(10+(int) ceil(urandom(0,20)), (int) ceil(urandom(0,10)));
        if (!res) {
            
            errors ++;
        }
    }
    fflush(stdout);
    fflush(stderr);
    if (errors) {
        std::cerr << "test failed with " << errors << "/" << n_tests << " failed tests\n";
    } else {
        std::cout << "test complete with no errors in " << n_tests << " tests\n";
    }
    return errors;
}

bool TestRing (int number_of_entries, int ring_size)
{
    cout << "entries: " << number_of_entries << " - size: " << ring_size << endl;
    Ring<int> ring(ring_size);
    list<int> X;
    for (int i=0; i<number_of_entries; ++i) {
        X.push_back(i);//(int) floor(10*urandom()));
        ring.push_back(i);
    }
    Ring<int>::iterator rit = ring.begin();
    list<int>::reverse_iterator lit = X.rbegin();

    bool flag = true;
    for (int i=0; i<ring_size; ++i, ++rit, ++lit) {
        cout << i << " " << *rit << " " << *lit << endl;
        if (*rit != *lit) {
            flag = false;
            //break;
        }
    }
        
    return flag;
}

#endif
