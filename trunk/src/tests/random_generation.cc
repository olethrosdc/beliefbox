#include "Distribution.h"
#include <vector>
#include "real.h"
int main()
{
	int T=10000;
	int N=10000;
	std::vector<real> s(T);
	
	for (int j=0; j<T; j++) {
		s[j] = 0.0f;
	}

	real prior = 0.5;
	srand48(time(NULL));
	for (int i=0; i<N; i++) {
		bool c = false;
		real p = 0.75;
		if (urandom()<prior) {
			c = true;
		} 

		for (int j=0; j<T; j++) {
			real r = 0.0f;
			if (c) {
				r = 0.5f;
			} else {
				r = 0.75f;
				//if (urandom()<p) {
				//	r = 1.0f;
				//}
			}
			s[j] += r;
		}
	}
	for (int j=0; j<T; j++) {
		s[j] /= (real) N;
		printf ("%f\n", s[j]);
	}
	return 0;
}
