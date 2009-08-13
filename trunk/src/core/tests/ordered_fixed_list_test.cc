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
    OrderedFixedList L(K);
    list<real> X;
    for (int i=0; i<N; ++i) {
        real x = floor(10*urandom());
        X.push_back(x);
        L.AddPerhaps(x);
    }
    X.sort();
    list<real>::iterator it; 
    list<real>::iterator oit; 

    it = X.begin();
    oit = L.S.begin();
    bool flag = true;
    for (int i=0; i<K; ++i, ++it, ++oit) {
        if (it == X.end() || oit == L.S.end()) {
            break;
        }
        if (*it != *oit) {
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
            if (approx_eq(x, y,0.1)) {
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
