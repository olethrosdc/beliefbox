#ifdef MAKE_MAIN
#include "ContextTreeKDTree.h"
#include "ContextTreeRealLine.h"
#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "MultivariateNormal.h"
#include "MultivariateNormalUnknownMeanPrecision.h"
#include "BetaDistribution.h"
#include "iWishart.h"


int main (int argc, char** argv)
{
	int n_dim = 2;
	Matrix S0 = Matrix::Unity(n_dim, n_dim)*10;
	int N =  n_dim + 1;
    iWishart iwish = iWishart(N, S0, true);
    Wishart wish = Wishart(N, S0, true);
	int n = 100000;
	Matrix AV(n_dim, n_dim);
	Matrix iAV(n_dim, n_dim);
	for(int i = 0; i<n; ++i) {
        Matrix iV = iwish.generate();
		Matrix V = wish.generate();
	    //printf("iWish\n"); iV.print(stdout);
		//printf(" Wish\n"); V.print(stdout);
		iAV += iV;
		AV += V;
	}
	iAV *= 1.0/(real) n;
	AV *= 1.0/(real) n;
	printf("A iWish\n"); iAV.print(stdout);
	printf("A  Wish\n"); AV.print(stdout);

    return 0;
}

#endif
