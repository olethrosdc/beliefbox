/* LU decomposition.
 *
 * Perform LU factorisation, and LU-based inverse and determinant.
 *
 * Notes:
 *	(1) When computing the LU decomposition for row-major matrices, it is necessary
 * to transpose the matrix before and after the operation (to be compatible with the
 * LAPACK column-major routines). However, when computing the matrix inverse, *neither* 
 * transposition is required. This is because trans(inv(A)) is equal to inv(trans(A)).
 *
 * TODOs:
 *	(1) Confirm that the second statement of Note (1) is correct. The results seem
 * reasonable, but there is some discrepancy between a row-major and a column-major
 * answer. This is probably just differing roundoff errors.
 *	(2) Implement reciprocal condition estimator: dgecon_
 *	(3) Implement an LU-based solver: dgetrs_
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_ULFACTOR_HPP_
#define ULAPACK_ULFACTOR_HPP_

#include "lapack_exception.hpp"
#include "errormacros.hpp"

namespace ulapack {
	namespace ublas = boost::numeric::ublas;

	namespace detail {

		extern "C" 
		{
		// LAPACK function declarations

		// LU decomposition of a real square matrix
		void dgetrf_( const int *m, const int *n, double *a, 
			const int *lda, const int *ipiv, int *info );

		// Inverse of LU factorised matrix
		void dgetri_(const int *n, double *a, const int *lda, const int *ipiv,
			double *work, const int *lwork, int *info );

		} // extern "C" 

		// LU inversion uses a block algorithm, which is more efficient over 
		// a certain matrix size, CUTOFF.
		enum { CUTOFF=64, BLOCKSIZE=64 }; 

		template<class F, class A>
		int lufactor_basic(ublas::matrix<double,F,A> &m, ublas::vector<int> &pivot)
		// Perform LU decomposition (in-place). Returns 0 if ok, LAPACK info if singular.
		// This "basic" function produces correct LU decompose for column-major matrices only.
		{
			// Call LAPACK routine
			const int M = static_cast<int>(m.size1());
			const int N = static_cast<int>(m.size2());
			int info;
			detail::dgetrf_( &M, &N, &m.data()[0], &M, &pivot.data()[0], &info );

			// Check validity of result
			if (info < 0) 
				throw LogicalError(ERROR_INFO("Invalid argument"), info);
			return info;
		}

		template<class F, class A>
		void inverse_factored(ublas::matrix<double,F,A> &m, 
			const ublas::vector<int> &pivot, ublas::vector<double> &work)
		// Compute inverse (in-place) of a matrix that has already been decomposed by lufactor_basic.
		// The pivot vector values are assumed to correspond to m.
		{
			// Set work vector size
			size_t bsize = m.size1();
			if (m.size1() > detail::CUTOFF) // use block-size > 1 if m large
				bsize *= detail::BLOCKSIZE; 
			if (work.size() < bsize)
				work.resize(bsize);

			// Call LAPACK inversion routine
			const int N = static_cast<int>(m.size1());
			const int wsize = static_cast<int>(work.size());
			int info;
			detail::dgetri_(&N, &m.data()[0], &N, &pivot.data()[0],
				&work.data()[0], &wsize, &info);

			// Check validity of result
			if (info < 0) 
				throw LogicalError(ERROR_INFO("Invalid argument"), info);
			else if (info > 0) 
				throw NumericalError(ERROR_INFO("Matrix is singular"), info);
		}

		template<class F, class A>
		double determinant_factored(const ublas::matrix<double,F,A> &m,
			const ublas::vector<int> &pivot)
		// Compute matrix determinant of a matrix that has already been decomposed by lufactor_basic.
		// The pivot vector values are assumed to correspond to m.
		{
			double det = 1.;
			for (size_t i=0; i < m.size1(); ++i) {
				if (pivot(i) != i+1) 
					det = -det;
				det *= m(i,i);
			}

			return det;
		}

	} // namespace detail

	// --------------------------------------------------------------------------
	// --------------------------- In-place algorithms --------------------------
	// --------------------------------------------------------------------------

	template<class A>
	bool lufactor(ublas::matrix<double, ublas::column_major, A> &m, 
		ublas::vector<int> &pivot)
	// Perform LU decomposition for column-major matrices (in-place).
	// Returns false if matrix is singular.
	{
		const size_t N = std::min(m.size1(), m.size2());
		if (pivot.size() < N)
			pivot.resize(N);
		return detail::lufactor_basic(m, pivot) == 0;
	}

	template<class A>
	bool lufactor(ublas::matrix<double, ublas::row_major, A> &m, 
		ublas::vector<int> &pivot)
	// Perform LU decomposition for row-major matrices (in-place).
	// Returns false if matrix is singular.
	// Warning: row-major lu-decomposition is less efficient than column-major.
	{
		ublas::matrix<double, ublas::column_major, A> tmp(m);
		bool status = lufactor(tmp, pivot); // call col-major version
		m = tmp;
		return status;
		// Note, we copy to column-major here rather than perform transposes explicitly,
		// as this gives us the correct size1() and size2() values for lufactor_basic()
		// whereas the m=trans(m) operations would not. Also, it gives us the right pivot.
	} 

	template<class F, class A>
	void inv_inplace(ublas::matrix<double,F,A> &m)
	// Compute matrix inverse (in-place).
	{
		if (m.size1() != m.size2())
			throw LogicalError(ERROR_INFO("Matrix is not square"));

		// LU decompose
		ublas::vector<int> pivot(m.size1());
		int info = detail::lufactor_basic(m, pivot);
		if (info != 0) 
			throw NumericalError(ERROR_INFO("Matrix is singular"), info);

		// Invert
		ublas::vector<double> work;
		detail::inverse_factored(m, pivot, work);

		// Check that we are using an optimal block size (useful check in pre-release builds)
		assert(work.size() >= work.data()[0]); 
	}

	// --------------------------------------------------------------------------
	// ------------------------- Non-in-place algorithms ------------------------
	// --------------------------------------------------------------------------

	template<class F, class A>
	ublas::matrix<double,F,A> inv(const ublas::matrix<double,F,A> &m)
	// Compute matrix inverse.
	{
		ublas::matrix<double,F,A> tmp(m);
		inv_inplace(tmp);
		return tmp;
	}

	template<class F, class A>
	double det(const ublas::matrix<double,F,A> &m)
	// Compute matrix determinant
	{
		if (m.size1() != m.size2())
			throw LogicalError(ERROR_INFO("Matrix is not square"));

		ublas::matrix<double,F,A> lu(m);
		ublas::vector<int> pivot(lu.size1());
		if (detail::lufactor_basic(lu, pivot) != 0)
			return 0.; // matrix singular implies det == 0

		return detail::determinant_factored(lu, pivot);
	}

	// Wrap everything up as a class, templatised by M - the matrix type.
	// Class caches the essential storage: LU matrix, pivot, and working.
	template<class M>
	class LU {
	public:
		LU(const M &m) : lu_(m), pivot_(m.size1()), 
			info(detail::lufactor_basic(lu_, pivot_)) {}

// Unsafe operations: (they give incorrect results if M is row-major)
//		const M &lu() const { return lu_; }	// Warning, only correct for M column-major
//		const ublas::vector<int> &pivot() const { return pivot_; } // Warning, only correct for M column-major

		bool is_singular() const { return info != 0; }

		M inv() const {
			if (lu_.size1() != lu_.size2())
				throw LogicalError(ERROR_INFO("Matrix is not square"));
			if (is_singular()) 
				throw NumericalError(ERROR_INFO("Upper matrix is singular"), info);

			M invmat(lu_);
			ublas::vector<double> work;
			detail::inverse_factored(invmat, pivot_, work);
			return invmat;
		}

		double det() const 
		{ 
			if (lu_.size1() != lu_.size2())
				throw LogicalError(ERROR_INFO("Matrix is not square"));
			if (is_singular()) 
				return 0.; 
			return detail::determinant_factored(lu_, pivot_);
		}

	private:
		M lu_;
		ublas::vector<int> pivot_;
		const int info;
	};

} // namespace ulapack

#endif
