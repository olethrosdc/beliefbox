/* Matrix symmetry functions.
 *
 * The function is_symmetric() returns true if the matrix is symmetric.
 * The function force_symmetry() enforces symmetry of a square matrix by
 * copying either the upper triangle to the lower or vice-versa (upper to
 * lower by default).
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_SYMMETRY_HPP_
#define ULAPACK_SYMMETRY_HPP_

#include "lapack_exception.hpp"
#include "errormacros.hpp"

namespace ulapack {
	namespace ublas = boost::numeric::ublas;

	template<class F, class A>
	bool is_symmetric(const ublas::matrix<double, F, A> &m)
	{
		if (m.size1() != m.size2())
			return false;

		for (size_t i = 0; i < m.size1(); ++i)
			for (size_t j = i+1; j < m.size2(); ++j)
				if (m(i,j) != m(j,i))
					return false;

		return true;
	}

	template<class F, class A>
	void force_symmetry(ublas::matrix<double, F, A> &m, const bool upperToLower=true)
	// Make matrix symmetric by copying upper-to-lower (default) or lower-to-upper.
	{
		if (m.size1() != m.size2())
			throw LogicalError(ERROR_INFO("Matrix is not square"));

		if (upperToLower) {
            for (size_t i = 0; i < m.size1(); ++i)
				for (size_t j = i+1; j < m.size2(); ++j)
					m(j,i) = m(i,j);
		}
		else {
            for (size_t i = 0; i < m.size1(); ++i)
				for (size_t j = i+1; j < m.size2(); ++j)
					m(i,j) = m(j,i);
		}
	}
}

#endif
