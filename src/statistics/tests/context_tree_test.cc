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

#ifdef MAKE_MAIN
#include "ContextTree.h"
#include "ContextTreeCTW.h"
#include "ContextTreePPM.h"
#include "ContextTreeBMC.h"
#include "Random.h"
#include "RandomNumberGenerator.h"
#include "MersenneTwister.h"
#include "ReadFile.h"
//#include "ResourceUse.h"
#include "EasyClock.h"
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>

class AbstractContextTreeInterface
{                                               
public:
    AbstractContextTreeInterface()
    {}
    virtual ~AbstractContextTreeInterface()
    {}
    virtual real Observe(int x, int y) = 0;
};

template <class T>
class ContextTreeInterface : public AbstractContextTreeInterface
{
protected:
    T* tree;
public:
    ContextTreeInterface(int n_branches_, int n_symbols_, int max_depth_ ) 
    //    : tree(n_branches_, n_symbols_, max_depth_)
    {
        tree = new T(n_branches_, n_symbols_, max_depth_);
    }
    virtual ~ContextTreeInterface()
    {
        delete tree;
    }
    virtual real Observe(int x, int y) {
        return tree->Observe(x, y);
    }
};

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: context_tree_test method depth output [input_file [T]]\n";
        exit(-1);
    }

    std::string method_name(argv[1]);
	int depth = atoi(argv[2]);
    if (depth <= 0) {
        std::cerr << "Depth should be > 0\n";
        exit(-1);
    }

    std::string out_fname(argv[3]);
    std::ofstream out_file;
    out_file.open(out_fname.c_str());
    
	int n_symbols = 2;
	MersenneTwisterRNG mt;
	RandomNumberGenerator* rng = &mt;
	rng->manualSeed(123456791);
	int T = 10000;
	std::vector<int> data(T);
	if (argc==4) {
		for (int t=0; t<T; ++t) {
			data[t] = rng->discrete_uniform(n_symbols);
		} 
	} else if (argc>=5) {
		if (argc==6) {
			T = atoi(argv[5]);
		} else {
			T = 0;
		}
		n_symbols = FileToIntVector(data, argv[4], T);
		T = data.size();
	} 
	

	AbstractContextTreeInterface* tree;
    if (!method_name.compare("BVMM")) {
        tree = new ContextTreeInterface<ContextTree> (n_symbols, n_symbols, depth);
    } else  if (!method_name.compare("CTW")) {
        tree = new ContextTreeInterface<ContextTreeCTW> (n_symbols, n_symbols, depth);
    } else  if (!method_name.compare("PPM")) {
        tree = new ContextTreeInterface<ContextTreePPM> (n_symbols, n_symbols, depth);
    } else  if (!method_name.compare("BMC")) {
        tree = new ContextTreeInterface<ContextTreeBMC> (n_symbols, n_symbols, depth);
    } else {
        std::cerr << "Unknown method " << method_name << std::endl;
        exit(-1);
    }

	std::cout << std::endl;
    double start_time = GetCPU();
    std::vector<real> p(T);
    int x = 0;
    real logsum = 0;
    real log2 = log(2);
	for (int t=0; t<T; ++t) {
		int y = data[t];
        p[t] = tree->Observe(x, y);
        x = y;
        logsum += log(p[t]) / log2;
	}
    double end_time = GetCPU();
	//tree.Show();

    delete tree;

    std::cout << "Depth: " << depth;
    std::cout << ", log loss: " << logsum / (real) T;
    std::cout << ", time: " << end_time - start_time << std::endl;
#if 0
    std::cout << ", RSS: " << getMaxResidentSetSize()
              << ", SHR: " << getSharedMemorySize() 
              << ", DAT: " << getUnsharedDataSetSize ()
              << ", STACK: " << getUnsharedStackSize()
              << std::endl;
#endif
    for (int t=0; t<T; ++t) {
        out_file << p[t] << std::endl;
    }
    out_file.close();


    return 0;
}

#endif
