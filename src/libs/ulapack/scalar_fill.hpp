/* Scalar fill routines.
 * 
 * Fill all of part of a matrix with a single scalar value. The possible
 * parts are: all, upper-triangular, lower-triangular, strict-upper-triangular, 
 * and strict-lower-triangular.
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_SCALAR_FILL_HPP_
#define ULAPACK_SCALAR_FILL_HPP_

namespace ulapack {

	namespace ublas = boost::numeric::ublas;

	namespace detail {
		
		// Define a series of indexing operations for various triangular forms
		template <class T>
		class TriOp;

		template <>
		class TriOp<ublas::upper> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j <= i; }
			size_t left_idx(size_t i, size_t j) { return j; }
			size_t right_idx(size_t i, size_t j) { return i; }
		};

		template <>
		class TriOp<ublas::lower> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j <= i; }
			size_t left_idx(size_t i, size_t j) { return i; }
			size_t right_idx(size_t i, size_t j) { return j; }
		};

		template <>
		class TriOp<ublas::strict_upper> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j < i; }
			size_t left_idx(size_t i, size_t j) { return j; }
			size_t right_idx(size_t i, size_t j) { return i; }
		};

		template <>
		class TriOp<ublas::strict_lower> {
		public:
			bool cmp_ij(size_t i, size_t j) { return j < i; }
			size_t left_idx(size_t i, size_t j) { return i; }
			size_t right_idx(size_t i, size_t j) { return j; }
		};

	} // namespace detail
    
	template<class F, class A>
	void scalar_fill(ublas::matrix<double, F, A> &m, double x)
	// Fill all m with scalar x.
	{
		double *p = &m.data()[0];
		const double * const end = p + m.size1() * m.size2();
		for (; p != end; ++p) 
			*p = x;
	}

	// Fill  triangular parts of m with scalar x.
	template<class TRI, class F, class A>
	void scalar_fill(ublas::matrix<double, F, A> &m, double x)
	{
		detail::TriOp<TRI> t;
		for (size_t i = 0; i < m.size1(); ++i)
			for (size_t j = 0; t.cmp_ij(i,j); ++j)
				m(t.left_idx(i,j), t.right_idx(i,j)) = x;
	} 

} // namespace ulapack

#endif
