/* Pretty printing of a uBLAS matrix.
 *
 * Tim Bailey 2005.
 */

#include <cstdio>

namespace ulapack {
	namespace ublas = boost::numeric::ublas;

	template<class F, class A>
	void print_matrix(const ublas::matrix<double,F,A> &m)
	{
		using namespace std;

		for (size_t i = 0; i < m.size1(); ++i) {
			for (size_t j = 0; j < m.size2(); ++j)
				printf("%9.3g\t", m(i,j));
			printf("\n");
		}
		printf("\n");
	}

} // namespace ulapack
