/* Perform Cholesky decomposition.
 *
 * For positive definite matrices, these functions compute:
 *	- cholesky decomposition (upper or lower triangular as specified)
 *	- inverse
 *	- determinant
 % Each operation comes in an in-place form and a copying form.
 * 
 * Notes: 
 * (1) The cholesky decomposition is for symmetric matrices, therefore it does not
 * matter whether the matrix is row-major or column major. However, for row-major
 * matrices, the result format (lower/upper triangular) is reversed.
 % (2) The LAPACK function dpotrf_() stores the result in the requested lower/upper
 * triangle, but does not zero-out the opposite triangular part.
 *
 * TODOs:
 *	(1) Implement reciprocal condition estimator: dpocon_
 *	(2) Implement a solver: dpotrs_
 *
 * Tim Bailey 2005.
 */

#ifndef ULAPACK_CHOLESKY_HPP_
#define ULAPACK_CHOLESKY_HPP_

#include "scalar_fill.hpp"
#include "symmetry.hpp"
#include "lapack_exception.hpp"
#include "errormacros.hpp"

namespace ulapack {
    namespace ublas = boost::numeric::ublas;

    namespace detail {

        extern "C" 
        {
            // LAPACK function declarations

            // Cholesky factorization of a real symmetric positive definite matrix
            void dpotrf_( const char *uplo, const int *n, double *a, 
                          const int *lda, int *info );

            // Inverse of Cholesky factorised matrix
            void dpotri_( const char *uplo, const int *n, double *a, 
                          const int *lda, int *info );
		
        } // extern "C" 

        template <class A>
        char uplo_flag(const ublas::matrix<double, ublas::row_major, A> &m, const bool upper)
            // Compute "reversed" upper/lower flags for row-major matrices
        {
            return (upper) ? 'L' : 'U';
        }

        template <class A>
        char uplo_flag(const ublas::matrix<double, ublas::column_major, A> &m, const bool upper)
            // Compute normal upper/lower flags for column-major matrices
        {
            return (upper) ? 'U' : 'L';
        }

        template<class F, class A>
        int cholesky_basic_checked(ublas::matrix<double,F,A> &m, const bool upper)
            // Perform Cholesky decomposition, but do not zero out other-triangular part.
            // Returns 0 if matrix was actual positive-definite, otherwise it returns a
            // LAPACK info value.
        {
            if (m.size1() != m.size2())
                throw LogicalError(ERROR_INFO("Matrix is not square"));
            assert(is_symmetric(m)); 

            // Call LAPACK routine
            int info;
            char uplo = detail::uplo_flag(m, upper);
            int size = static_cast<int>(m.size1());
            detail::dpotrf_( &uplo, &size, &m.data()[0], &size, &info );

            // Check validity of result
            if (info < 0) 
                throw LogicalError(ERROR_INFO("Invalid argument"), info);

            return info;
        }

        template<class F, class A>
        void cholesky_basic(ublas::matrix<double,F,A> &m, const bool upper)
            // Perform Cholesky decomposition, but do not zero out other-triangular part.
        {
            int info = cholesky_basic_checked(m, upper);
            if (info > 0) 
                throw NumericalError(ERROR_INFO("Matrix is not positive definite"), info);
        }

        template<class F, class A>
        void zero_strict_triangular(ublas::matrix<double,F,A> &m, const bool upper)
            // Zero out strict-other triangular part of matrix.
        {
            if (upper) 
                scalar_fill<ublas::strict_lower>(m, 0.);
            else 
                scalar_fill<ublas::strict_upper>(m, 0.);
        }

    } // namespace detail

    // --------------------------------------------------------------------------
    // --------------------------- In-place algorithms --------------------------
    // --------------------------------------------------------------------------

    template<class F, class A>
    bool chol_checked_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition (in-place).
	// Returns false if matrix is not positive-definite.
    {
        if (detail::cholesky_basic_checked(m, upper))
            return false;
        detail::zero_strict_triangular(m, upper);
        return true;
    }

    template<class F, class A>
    void chol_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition (in-place).
    {
        detail::cholesky_basic(m, upper);
        detail::zero_strict_triangular(m, upper);
    }

    template<class F, class A>
    void inv_chol_inplace(ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute inverse (in-place) of a pos. def. matrix that has previously been 
	// Cholesky decomposed (ie, "m" is already a Cholesky matrix). 
	// WARNINGS: 
	//	(1) The value of "upper" must match the actual form of "m". This function does not check.
	//	(2) This function does *not* return the inverse of a Cholesky matrix.
    {
        // Call LAPACK routine
        int info;
        char uplo = detail::uplo_flag(m, upper);
        int size = static_cast<int>(m.size1());
        detail::dpotri_( &uplo, &size, m.data().begin(), &size, &info );

        // Check validity of result
        if (info < 0) 
            throw LogicalError(ERROR_INFO("Invalid argument"), info);
        else if (info > 0) 
            throw NumericalError(ERROR_INFO("Inverse does not exist"), info);

        // Copy result to other triangular part (ie, make a symmetric inverse matrix)
        force_symmetry(m, upper);
    }

    template<class F, class A>
    void inv_pd_inplace(ublas::matrix<double,F,A> &m)
	// Compute inverse of a positive definite matrix (in-place).
    {
        detail::cholesky_basic(m, true);
        inv_chol_inplace(m, true);
    }

    // --------------------------------------------------------------------------
    // ------------------------- Non-in-place algorithms ------------------------
    // --------------------------------------------------------------------------

    template<class F, class A>
    ublas::matrix<double,F,A> 
    chol(const ublas::matrix<double,F,A> &m, const bool upper=true)
	// Compute Cholesky decomposition, and return result.
    {
        ublas::matrix<double,F,A> c(m);
        chol_inplace(c, upper);
        return c;
    }

    template<class F, class A>
    ublas::matrix<double,F,A> inv_pd(const ublas::matrix<double,F,A> &m)
	// Compute inverse of a positive definite matrix, and return the result.
    {
        ublas::matrix<double,F,A> inv(m);
        inv_pd_inplace(inv);
        return inv;
    }

    template<class F, class A>
    double det_chol(const ublas::matrix<double,F,A> &m)
	// Compute determinant of Cholesky matrix. 
    {
        assert(m.size1() == m.size2());
        double d = 1.;
        for (size_t i=0; i < m.size1(); ++i)
            d *= m(i,i);
        return d;
    }

    template<class F, class A>
    double det_pd_sqrt(const ublas::matrix<double,F,A> &m)
	// Compute sqrt of determinant of pos. def. matrix. 
    {
        ublas::matrix<double,F,A> c(m);
        detail::cholesky_basic(c, true);
        double dsqrt = det_chol(c);
        assert(dsqrt > 0.);
        return dsqrt;
    }

    template<class F, class A>
    double det_pd(const ublas::matrix<double,F,A> &m)
	// Compute determinant of pos. def. matrix.
    {
        double d = det_pd_sqrt(m);
        return d*d;
    }

    // Wrap everything up as a class, templatised by M - the matrix type.
    template<class M>
    class Cholesky {
    public:
        Cholesky(const M &m, const bool upper=true) :
            cholmat(m), first(true), is_upper(upper), 
            info(detail::cholesky_basic_checked(cholmat, is_upper))
        { }

        bool is_posdef() const { return info == 0; }

        const M& chol() const {
            check_posdef();
            if (first) { 
                first = false;
                detail::zero_strict_triangular(cholmat, is_upper);
            }
            return cholmat;
        }

        M inv() const {
            check_posdef();
            M invmat(cholmat);
            inv_chol_inplace(invmat, is_upper);
            return invmat;
        }

        double det_sqrt() const {
            check_posdef();
            return det_chol(cholmat);
        }

        double det() const {
            double d = det_sqrt();
            return d*d;
        }

    private:
        void check_posdef() const {
            if (info != 0)
                throw NumericalError(ERROR_INFO("Matrix is not positive definite"), info);
        }

        mutable M cholmat;
        mutable bool first;
        const bool is_upper;
        const int info;
    }; 

} // namespace ulapack

#endif
