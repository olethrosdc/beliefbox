#include "ReadFile.h"
#include <cstdlib>
#include <cstdio>
#include "debug.h"
#include <cerrno>

/** Read a file in.

	The format is:
	number_of_lines
 */
int FileToIntVector(std::vector<int>& data, const char* fname, int tmpT)
{
    FILE* file = fopen(fname, "r");
    if (!file) {
        fprintf (stderr, "Error: Could not open file %s\n", fname);
        exit(-1);
    }

	int T = 0;
    int success = fscanf(file, "%d", &T);
	if (success <= 0) {
            Serror("Could not scan file %s - T =%d - retval: %d - errno %d\n", fname, T, success, errno);
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

/** Read a list of class records.

	Format:
	number_of_lines number_of_colums
	<data>
 */
int ReadClassData(Matrix& data, std::vector<int>& labels, const char* fname) 
{
    FILE* file = fopen(fname, "r");
    if (!file) {
        fprintf (stderr, "Error: Could not open file %s\n", fname);
        exit(-1);
    }

	int T = 0;
    int columns;
    int success = fscanf(file, "%d %d", &T, &columns);
	if (success <= 0) {
            Serror("Could not scan file %s - T =%d - retval: %d - errno %d\n", fname, T, success, errno);
            exit(-1);
	}
    
    printf("# horizon: %d, columns: %d\n", T, columns);
    data.Resize(T, columns - 1);
    labels.resize(T);

    for (int t=0; t<T; ++t) {
        for (int i=0; i<columns - 1; ++i) {
            int success = fscanf(file, "%lf ", &data(t,i));
            //printf ("%f ", data(t,i));
            if (success <=0) {
                Serror("Could not scan file, line %d, column %d, suc: %d, errno: %d\n", t, i, success, errno);
                exit(-1);
            }
        }
        
        int success = fscanf(file, "%d", &labels[t]);
        if (success <=0) {
            Serror("Could not scan file, line %d\n", t);
            exit(-1);
        }
        //printf("%d\n", labels[t]);
    }
    int min_label = Min(labels);
    for (int t=0; t<T; ++t) {
        labels[t] -= min_label;
    }
	return T;
}

/** Read a matrix of data

	Format:
	number_of_lines number_of_colums
	<data>
    
    returns:
    number of lines read
 */
int ReadFloatDataASCII(Matrix& data, const char* fname) 
{
    FILE* file = fopen(fname, "r");
    if (!file) {
        fprintf (stderr, "Error: Could not open file %s\n", fname);
        exit(-1);
    }

	int T = 0;
    int columns;
    int success = fscanf(file, "%d %d", &T, &columns);
	if (success <= 0) {
            Serror("Could not scan file %s - T =%d - retval: %d - errno %d\n", fname, T, success, errno);
            exit(-1);
	}
    
    printf("# horizon: %d, columns: %d\n", T, columns);
    data.Resize(T, columns);


    for (int t=0; t<T; ++t) {
        for (int i=0; i<columns; ++i) {
            int success = fscanf(file, "%lf ", &data(t,i));
            //printf ("%f ", data(t,i));
            if (success <=0) {
                Serror("Could not scan file, line %d, column %d, suc: %d, errno: %d\n", t, i, success, errno);
                exit(-1);
            }
        }
    }
	return T;
}
