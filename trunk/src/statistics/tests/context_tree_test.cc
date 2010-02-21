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
#include "Random.h"
#include "RandomNumberGenerator.h"
#include "MersenneTwister.h"
#include "ReadFile.h"
#include "ResourceUse.h"
#include "EasyClock.h"
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>

/// TODO: Make method work with arguments
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: context_tree_test depth output [input_file [T]]\n";
        exit(-1);
    }
	int depth = atoi(argv[1]);
    if (depth <= 0) {
        std::cerr << "Depth should be > 0\n";
        exit(-1);
    }

    std::string out_fname(argv[2]);
    std::ofstream out_file;
    out_file.open(out_fname.c_str());
    
	int n_symbols = 2;
	MersenneTwisterRNG mt;
	RandomNumberGenerator* rng = &mt;
	rng->manualSeed(123456791);
	int T = 10000;
	std::vector<int> data(T);
	if (argc==3) {
		for (int t=0; t<T; ++t) {
			data[t] = rng->discrete_uniform(n_symbols);
		} 
	} else if (argc>=4) {
		if (argc==5) {
			T = atoi(argv[4]);
		} else {
			T = 0;
		}
		n_symbols = FileToIntVector(data, argv[3], T);
		T = data.size();
	} 
	

	ContextTree tree(n_symbols, n_symbols, depth);
	std::cout << std::endl;
    double start_time = GetCPU();
    std::vector<real> p(T);
    int x = 0;
	for (int t=0; t<T; ++t) {
		int y = data[t];
        p[t] = tree.Observe(x, y);
        x = y;
	}
    double end_time = GetCPU();
	//tree.Show();
    std::cout << "Depth: " << depth;
    std::cout << " time: " << end_time - start_time;
    std::cout << " RSS: " << getMaxResidentSetSize()
              << " SHR: " << getSharedMemorySize() 
              << " DAT: " << getUnsharedDataSetSize ()
              << " STACK: " << getUnsharedStackSize()
              << std::endl;
    for (int t=0; t<T; ++t) {
        out_file << p[t] << std::endl;
    }
    out_file.close();


    return 0;
}

#endif
