#include "ConnectivityMatrix.h"
#include <cstdlib>
#include <cstdio>

ConnectivityMatrix::ConnectivityMatrix(int N, bool directional, int** c, bool invert, real** d) : Graph (N, directional)
{
#if GRAPH_DEBUG_LEVEL > 90
	printf ("Creating Graph as connectivity matrix with %d nodes\n", N);
#endif
	C = new int* [N];
	D = new real* [N];
	for (int i=0; i<N; i++) {
		C[i] = new int[N];
		D[i] = new real[N];
		for (int j=0; j<N; j++) {
			if (c) { //copy c, otherwise make empty matrix
				if (invert) {
					C[i][j] = c[j][i];
					if (d) {
						D[i][j] = d[j][i];
					} else {
						D[i][j] = 1.0;
					}
				} else {
					C[i][j] = c[i][j];
					if (d) {
						D[i][j] = d[i][j];
					} else {
						D[i][j] = 1.0;
					}
				}
			} else {
				C[i][j] = 0;
				D[i][j] = 0.0;
			}
		}
	}

	// make matrix symettric
	if (directional==false) {
		for (int i=0; i<N; i++) {
			for (int j=i+1; j<N; j++) {
				int cij = (C[i][j] + C[j][i]);
				if (cij>1) cij=1;
				C[j][i] = C[i][j] = cij;

			}
		}
	}
}



ConnectivityMatrix::ConnectivityMatrix(SparseGraph& G) : Graph (G.n_nodes(), G.is_directional())
{
#if GRAPH_DEBUG_LEVEL > 90
	printf ("Creating Graph as connectivity matrix from sparse matrix with %d nodes\n", N);
#endif
	C = new int* [N];
	D = new real* [N];
	for (int i=0; i<N; i++) {
		C[i] = new int[N];
		D[i] = new real[N];
		for (int j=0; j<N; j++) {
            C[i][j] = 0;
            D[i][j] = 0.0;
        }
	}
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
            if (G.edge(i,j)) {
                C[i][j] = 1;
                D[i][j] = G.distance(i,j);
            }
        }
    }

	// make matrix symettric
	if (directional==false) {
		for (int i=0; i<N; i++) {
			for (int j=i+1; j<N; j++) {
				int cij = (C[i][j] + C[j][i]);
				if (cij>1) cij=1;
				C[j][i] = C[i][j] = cij;

			}
		}
	}
}


bool ConnectivityMatrix::CheckSymmetry()
{
	if (!directional) {
		for (int i=0; i<N; i++) {
			for (int j=i+1; j<N; j++) {
				if (C[i][j] != C[j][i]) {
					return false;
				}
			}
		}
	}
	return true;
}

ConnectivityMatrix::~ConnectivityMatrix()
{
	for (int i=0; i<N; i++) {
		delete [] C[i];
		delete [] D[i];
	}
	delete [] C;
	delete [] D;
}

bool ConnectivityMatrix::edge (int src, int dst)
{
#if GRAPH_DEBUG_LEVEL >= 95
	if ((directional==false)&&(CheckSymmetry()==false)) {
		fprintf (stderr, "Warning: Undirected graph should have symmetric connectivity matrix.\n");
	}
#endif

	if (C[src][dst])
		return true;
	return false;
}


real ConnectivityMatrix::distance (int src, int dst)
{
	if (!edge(src, dst)) {
		fprintf (stderr, "Warning: Calculate edge length for non-existent edge (%d %d)", src, dst);
	}
	return D[src][dst];
}
