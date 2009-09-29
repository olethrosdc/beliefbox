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

#include "OrderedFixedList.h"
#include "Random.h"
#include "real.h"
#include <iostream>
#include <list>
using namespace std;

bool TestList (int N, int K);

int main(void)
{
    int n_tests = 1000;
    int errors = 0;

    cout << "\nRunning ordered fixed list test\n";
    cout << "n_tests: " << n_tests << endl;

    for (int test = 0; test<n_tests; ++test) {
        bool res = TestList((int) ceil(urandom(0,1000)), (int) ceil(urandom(0,10)));
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

bool TestList (int N, int K)
{
    OrderedFixedList<int> L(K);
    list<real> X;
    vector<int> number(N);
    for (int i=0; i<N; ++i) {
        real x = floor(10*urandom());
        number[i] = i;
        X.push_back(x);
        L.AddPerhaps(x, &number[i]);
    }
    X.sort();
    list<real>::iterator it; 
    list<pair<real,int*> >::iterator oit; 

    it = X.begin();
    oit = L.S.begin();
    bool flag = true;
    for (int i=0; i<K; ++i, ++it, ++oit) {
        if (it == X.end() || oit == L.S.end()) {
            break;
        }
        if (*it != oit->first) {
            flag = false;
            break;
        }
    }
        
    if (!flag) {
        it = X.begin();
        oit = L.S.begin();
        for (int i=0; i<K; ++i, ++it, ++oit) {
            real x = *it;
            real y = *it;
            if (approx_eq(x, y)) {
                cout <<  x << "=" << y << " ";
            } else {
                cout <<  x << "!" << y << " ";
            }
        }   
        cout << endl;
    }

    return flag;
}

#endif
