#include "ReadFile.h"
#include <cstdlib>
#include <cstdio>
#include "debug.h"

int FileToIntVector(std::vector<int>& data, char* fname, int tmpT) {
    FILE* file = fopen(fname, "r");
    if (!file) {
        fprintf (stderr, "Error: Could not open file %s\n", fname);
        exit(-1);
    }

	int T = 0;
    int success = fscanf(file, "%d", &T);
	if (success <= 0) {
		Serror("Could not scan file\n");
	}

	if (tmpT > 0 && tmpT < T) {
		T = tmpT;
	}

    printf("horizon: %d\n", T);
    data.resize(T);
    int n_observations = 0;
    for (int t=0; t<T; ++t) {
		int success = fscanf(file, "%d", &data[t]);
		if (success <=0) {
			Serror("Could not scan file\n");
		}
		if (data[t] > n_observations) {
            n_observations = data[t];
        }
        data[t] -= 1;
    }
    fclose(file);
	return n_observations;
}
