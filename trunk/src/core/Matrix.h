/* -*- Mode: c++ -*- */
// copyright (c) 2006-2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include "real.h"
#include "Vector.h"
#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
//#include <atlas/cblas.h>



#undef REFERENCE_ACCESS

#define ACCURACY_LIMIT 1e-12

/** \brief An n-by-m dimensional matrix.

    TODO: Use the BLAS interface for some / most routines.
 */
class Matrix
{
public:
    /// Check bounds?
    enum BoundsCheckingStatus {NO_CHECK_BOUNDS=0, CHECK_BOUNDS=1};
#ifdef NDEBUG
    static Matrix Null(int rows, int columns, enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
    static Matrix Unity(int rows, int columns, enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
	
    Matrix (int rows = 1, int columns = 1,  enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
    Matrix (int rows, int columns, real* y, enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
    explicit Matrix (const Vector& v, enum BoundsCheckingStatus check_ = NO_CHECK_BOUNDS);
#else
    static Matrix Null(int rows, int columns, enum BoundsCheckingStatus check = CHECK_BOUNDS);
    static Matrix Unity(int rows, int columns, enum BoundsCheckingStatus check = CHECK_BOUNDS);

    Matrix (int rows = 1, int columns = 1,  enum BoundsCheckingStatus check = CHECK_BOUNDS);
    Matrix (int rows, int columns, real* y, enum BoundsCheckingStatus check = CHECK_BOUNDS);
    explicit Matrix (const Vector& v, enum BoundsCheckingStatus check_ = CHECK_BOUNDS);
#endif
    Matrix (const Matrix& rhs, bool clone = true);
    ~Matrix();
    void Resize(int rows_, int columns_);
	Matrix AddRow(const Vector& rhs);
	Matrix AddColumn(const Vector& rhs);
    Matrix& operator= (const Matrix& rhs);
    bool operator== (const Matrix& rhs) const;
    bool operator!= (const Matrix& rhs) const;
    Matrix operator+ (const Matrix& rhs);
    Matrix& operator+= (const Matrix& rhs);
    Matrix operator- (const Matrix& rhs);
    Matrix& operator-= (const Matrix& rhs);
    Matrix& operator*= (const real& rhs);
    Matrix operator* (const Matrix& rhs);
    Matrix operator* (const real& rhs);
    Matrix operator+ (const real& rhs);
    Matrix operator- (const real& rhs);
    Matrix operator/ (const real& rhs);
    /// Matrix inversion (defaults to GSL with LU)
    Matrix Inverse(real epsilon = ACCURACY_LIMIT) const
    {
		return GSL_Inverse();
		//return Inverse_LU();
	}
	Vector SVD_Solve(const Vector& b) const;

	Matrix GSL_Inverse() const;

    /** Matrix inversion using the Cholesky decomposition.
        
        First call Cholesky with accuracy epsilon to obtain $\fL,
        L'\f$ such that \f$A = LL'\f$.  Then solve for \f$X = A^{-1}\f$:
        \f[
        LL'X = I,
        \f]
        by dynamic programming.
    */
    Matrix Inverse_Cholesky(real epsilon = ACCURACY_LIMIT) const
    {
        Matrix U = Cholesky(epsilon);
        Matrix L(U, false);
        L.Transpose();
        return Inverse(L, U);
    }
    /** Matrix inversion using the LU decomposition.
        
        First call LUDecomposition with accuracy epsilon to obtain $\fL,
        U\f$ such that \f$A = LU\f$.  Then solve for \f$X = A^{-1}\f$:
        \f[
        LUX = I,
        \f]
        by dynamic programming.
    */
    Matrix Inverse_LU(real epsilon = ACCURACY_LIMIT) const
    {
        real det;
        Matrix tmp(*this);
        std::vector<Matrix> A = tmp.LUDecomposition(det, epsilon);
        return Inverse(A[0], A[1]);
    }
    
    Matrix Inverse(const Matrix& L, const Matrix& U) const;
    bool isSymmetric() const;
    bool isTriangular() const;
    bool isUpperTriangular() const;
    bool isLowerTriangular() const;
    real det() const;
	real tr() const;
    real L1Norm() const;
    real L2Norm() const;
	Matrix Multiple(const Matrix& rhs) const;
	Matrix Kron(const Matrix& rhs) const;
    real ColumnSum(int c) const;
    real RowSum(int r) const;
    real compute_det_triangular() const;
    real gaussian_elimination_forward(real epsilon = ACCURACY_LIMIT);
	std::vector<Matrix> QRDecomposition() const;
    std::vector<Matrix> LUDecomposition(real& determinant, real epsilon = ACCURACY_LIMIT);
    Matrix Cholesky(real epsilon = ACCURACY_LIMIT) const;
    void Cholesky(Matrix& chol, real epsilon = ACCURACY_LIMIT) const;
    void Clear();
    void Transpose();
    Vector getColumn(int c) const;
    Vector getRow(int r) const;
    void setColumn(int c, const Vector& x);
    void setRow(int r, const Vector& x);
    void SortRow(int r);
    void SortColumn(int r);
    int Rows() const;
    int Columns() const;
	Vector Vec() const;
	void Vec(const Vector& x);
	Vector RowMax() const;
	Vector ColumnMax() const;
	Vector RowMin() const;
	Vector ColumnMin() const;
    real& operator() (int i, int j);
    const real& operator() (int i, int j) const;
    void print(FILE* f) const;
    friend Matrix operator* (const real& lhs, const Matrix& rhs);
    friend Matrix operator* (const Vector& lhs, const Matrix& rhs);
    friend Vector operator* (const Matrix& lhs, const Vector& rhs);
protected:
    int rows; ///< number of rows in the matrix
    int columns; ///< number of columns in the matrix
    real* x; ///< data
#ifdef REFERENCE_ACCESS
    real** x_list; ///< data pointers
    void MakeReferences();
#endif
    const enum BoundsCheckingStatus checking_bounds;
    bool transposed;
    bool clear_data;
    const real& qGet(int i, int j);
    void qSet(int i, int j, real v);
};

Matrix operator* (const real& lhs, const Matrix& rhs);
Matrix operator* (const Vector& lhs, const Matrix& rhs);
Vector operator* (const Matrix& lhs, const Vector& rhs);
real Mahalanobis2 (const Vector& x, const Matrix& A, const Vector& y);

/// Kronecker product
inline Matrix Kron (const Matrix& lhs, const Matrix& rhs)
{
	return lhs.Kron(rhs);
}

inline Matrix OuterProduct (const Vector& lhs, const Vector& rhs)
{
    Matrix R(lhs.n, rhs.n);
    MatrixProduct(lhs, rhs, R);
    return R;
}

Matrix Transpose (const Matrix& rhs);

/// In-place transpose matrix
inline
void Matrix::Transpose()
{
    transposed = !transposed;
}

inline
int Matrix::Rows() const
{
    if (transposed) {
        return columns;
    }
    return rows;
}

inline
int Matrix::Columns() const
{
    if (transposed) {
        return rows;
    } 
    return columns;
}

inline
real& Matrix::operator() (int i, int j)
{
    if (transposed) {
        int tmp = i;
        i = j;
        j = tmp;
    }
                 
    if (checking_bounds) {
        if (i<0 || j<0 || i>=rows || j>=columns) {
	    fprintf(stderr, "bad access: (%d, %d) on a %d x %d matrix\n",
		    i, j, rows, columns);
	    fflush (stderr);
        throw std::out_of_range("matrix index out of range");
        }
    }
#ifdef REFERENCE_ACCESS
    return x_list[i][j];
#else
    return x[i*columns + j];
#endif
}

inline
const real& Matrix::operator() (int i, int j) const
{
    if (transposed) {
        int tmp = i;
        i = j;
        j = tmp;
    }
                 
    if (checking_bounds) {
        if (i<0 || j<0 || i>=rows || j>=columns) {
	    fprintf(stderr, "bad access: (%d, %d) on a %d x %d matrix\n",
		    i, j, rows, columns);
	    fflush (stderr);
            throw std::out_of_range("matrix index out of range");
        }
    }
#ifdef REFERENCE_ACCESS
    return x_list[i][j];
#else
    return x[i*columns + j];
#endif
}

#if 0
inline
const real& Matrix::qGet(int i, int j)
{
#ifdef REFERENCE_ACCESS
    return x_list[i][j];
#else
    return x[i*columns + j];
#endif
}

inline
void Matrix::qSet(int i, int j, real v)
{
#ifdef REFERENCE_ACCESS
    x_list[i][j] = v;
#else
    x[i*columns + j] = v;
#endif
}
#endif

#endif

