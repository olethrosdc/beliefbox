/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Matrix.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <exception>
#include <stdexcept>
//#include <lapackpp/laslv.h>

Matrix Matrix::Unity(int rows, int columns, enum BoundsCheckingStatus check)
{
    int N = rows*columns;
    real* y = (real*) calloc(N, sizeof(real));
    int K = std::min(rows, columns);
    int step = columns + 1;
    int j = 0;
    for (int i=0; i<K; ++i) {
        y[j] = 1.0;
        j += step;
    }
    return Matrix(rows, columns, y, check);
}

Matrix Matrix::Null(int rows, int columns, enum BoundsCheckingStatus check)
{
    return Matrix(rows, columns, check);
}

Matrix Transpose(const Matrix& rhs)
{
    Matrix tmp = Matrix(rhs, false);
    tmp.Transpose();
    return tmp;
}

Matrix::Matrix (int rows_, int columns_,  enum BoundsCheckingStatus check_)
    : rows(rows_),
      columns(columns_),
      checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
    int N = rows*columns;
    x = (real*) calloc(N, sizeof(real));
}

Matrix::Matrix (int rows_, int columns_, real* y, enum BoundsCheckingStatus check_)
    : rows(rows_),
      columns(columns_),
      x(y),
      checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
}

Matrix::Matrix (const Vector& v, enum BoundsCheckingStatus check_)
    : checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
    rows = v.Size();
    columns = 1;
    x = (real*) malloc(rows * sizeof(real));
    for (int i=0; i<rows; ++i) {
        x[i] = v[i];
    }
}



/// Copy constructor
Matrix::Matrix (const Matrix& rhs, bool clone)
    : rows(rhs.Rows()),
      columns(rhs.Columns()),
      checking_bounds(rhs.checking_bounds),
      transposed(false)
{

    const int M = rows;
    const int N = columns;
    const int K = M*N;

    if (clone) {
        x = (real*) malloc (sizeof(real)*K);
		
        for (int m=0; m<M; ++m) {
            for (int n=0; n<N; ++n) {
                (*this)(m,n) = rhs(m,n);
            }
        }
        clear_data = true;
    } else {
        x = rhs.x;
        clear_data = false;
    }
	
}

Matrix::~Matrix()
{
    if (clear_data) {
        free(x);
    }
}



/// Assign another matrix to this
Matrix& Matrix::operator= (const Matrix& rhs)
{
    if (this == &rhs) return *this;

    transposed = false;
    columns = rhs.Columns();
    rows = rhs.Rows();

    const int M = rows;
    const int N = columns;
    const int K = M*N;
    
    x = (real*) realloc (x, sizeof(real)*K);

    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) = rhs(m,n);
        }
    }
    return *this;
}

/// Resize matrix
void Matrix::Resize (int rows_, int columns_)
{
    transposed = false;
    columns = columns_;
    rows = rows_;

    const int M = rows;
    const int N = columns;
    const int K = M*N;
    
    x = (real*) realloc (x, sizeof(real)*K);

    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) = 0.0;
        }
    }
}


/// Boolean equality operator
bool Matrix::operator== (const Matrix& rhs) const
{
    if (Rows()!=rhs.Rows() || Columns()!=rhs.Columns()) {
        throw std::domain_error("Matrix equality operator");
    }
    for (int i=0; i<Rows(); ++i) {
        for (int j=0; j<Columns(); ++j) {
            if ((*this)(i, j) != rhs(i, j)) {
                printf ("%f - %f=%f\n", (*this)(i, j), rhs(i, j), (*this)(i, j) - rhs(i, j));
                return false;
            }
        }
    }
    return true;
}

/// Boolean inequality operator
bool Matrix::operator!= (const Matrix& rhs) const
{
    return !(*this==rhs);
}

/// Create a matrix through the addition of two other matrices.
Matrix Matrix::operator+ (const Matrix& rhs)
{
    if (Columns() != rhs.Columns() 
        || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix addition error\n");
    }
    int M = Rows();
    int N = Columns();
    
    Matrix lhs(M, N);
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            lhs(m, n) = (*this)(m,n) + rhs(m,n);
        }
    }
    return lhs;
}

/// Add another matrix to this
Matrix& Matrix::operator+= (const Matrix& rhs)
{
    if (Columns() != rhs.Columns() 
        || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix addition error\n");
    }
    int M = Rows();
    int N = Columns();
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) += rhs(m,n);
        }
    }
    return *this;
}

/// Create a matrix through the subtraction of two other matrices.
Matrix Matrix::operator- (const Matrix& rhs)
{
    if (Columns() != rhs.Columns() 
        || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix addition error\n");
    }
    int M = Rows();
    int N = Columns();
    
    Matrix lhs(M, N);
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            lhs(m, n) = (*this)(m,n) - rhs(m,n);
        }
    }
    return lhs;
}


/// Subtract another matrix from this
Matrix& Matrix::operator-= (const Matrix& rhs)
{
    if (Columns() != rhs.Columns() 
        || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix addition error\n");
    }
    int M = Rows();
    int N = Columns();
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) -= rhs(m,n);
        }
    }
    return *this;
}



/// Create a matrix through the multiplication of two other matrices.
Matrix Matrix::operator* (const Matrix& rhs)
{
    if (Columns() != rhs.Rows()) {
        throw std::domain_error("Matrix multiplication error\n");
    }
    int M = Rows();
    int K = Columns();
    int N = rhs.Columns();
    
    Matrix lhs(M, N);
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            real sum = 0.0;
            for (int k=0; k<K; ++k) {
                sum += (*this)(m, k)*rhs(k, n);
            }
            lhs(m, n) = sum;
        }
    }
    return lhs;
}

/// Multiply a matrix with a scalar
Matrix& Matrix::operator*= (const real& rhs)
{
    real* y = x;
    int N = rows * columns;
    for (int n=0; n<N; ++n) {
        (*y++) *= rhs;
    }
    return *this;
}

/// Multiply a matrix with a scalar, creating a new matrix
Matrix Matrix::operator* (const real& rhs)
{
    Matrix lhs(rows, columns);
    if (transposed) {
        lhs.Transpose();
    }
    real* y = x;
    real* z = lhs.x;
    int N = rows * columns;
    for (int n=0; n<N; ++n) {
        (*z++) = (*y++) * rhs;
    }
    return lhs;
}


void Matrix::print(FILE* f)
{
    for (int i=0; i<Rows(); ++i) {
        for (int j=0; j<Columns(); ++j) {
            fprintf (f, "%f ", (*this)(i,j));
        }
        fprintf (f, "\n");
    }
}

/// Multiply a matrix with a scalar, creating a new matrix
Matrix operator* (const real& lhs, const Matrix& rhs)
{
    Matrix v(rhs.Rows(), rhs.Columns());
    if (rhs.transposed) {
        v.Transpose();
    }
    real* y = rhs.x;
    real* z = v.x;
    int N = rhs.Rows() * rhs.Columns();
    for (int n=0; n<N; ++n) {
        (*z++) = (*y++) * lhs;
    }
    return v;
}


/// Multiply a vector with a matrix, creating a new matrix
Matrix operator* (const Vector& lhs, const Matrix& rhs)
{
    if (rhs.Rows() != 1) {
        throw std::domain_error("Vector-matrix multiplication error\n");
    }
    Matrix v(lhs.Size(), rhs.Columns());
    if (rhs.transposed) {
        v.Transpose();
    }

    int N = rhs.Columns();
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            v(i,j) = lhs[i] * rhs(0, j);
        }
    }
    return v;
}

real Matrix::det()
{
    if (isTriangular()) {
        return compute_det_triangular();
    } else {
        Matrix A = *this;
        real swap = A.gaussian_elimination_forward();
        return swap * A.compute_det_triangular();
    }

}

bool Matrix::isTriangular()
{
    if (rows != columns) {
        throw std::domain_error("only square matrices can be triangular");
    }
    if (isLowerTriangular() || isUpperTriangular()) {
        return true;
    }
    return false;
}

bool Matrix::isUpperTriangular()
{
    if (rows != columns) {
        return false;
        throw std::domain_error("only square matrices can be triangular");
    }
    for (int j=0; j<Columns(); ++j) {
        for (int i=j; i<Rows(); ++i) {
            if ((*this)(i,j) != 0.0)
                return false;
        }
    }
    return true;
        
}

bool Matrix::isLowerTriangular()
{
    if (rows != columns) {
        throw std::domain_error("only square matrices can be triangular");
    }
    for (int i=0; i<Rows(); ++i) {
        for (int j=i; j<Columns(); ++j) {
            if ((*this)(i,j) != 0.0)
                return false;
        }
    }
    return true;
}

real Matrix::compute_det_triangular()
{
    Matrix& A = *this;
    real det = A(0,0);
    for (int i=1; i<rows; ++i) {
        det *= A(i,i);
    }
    return det;
}



/** Convert a matrix to an upper triangular matrix in-place.
 *
  * The method does no attempt to optimise the order of rows - it simply
  * selects rows such that the pivot element is non-zero.
  * All elements with norm < epsilon are considered to be zero.  
  *
  * This method throws std::runtime_error when the matrix appears to
  * be singular, i.e. when there is no way to find a useful pivot.
  *
  * It returns -1 if the number of row exchanges is odd, otherwise 1.
  * This is useful for calculating the determinant.
  *
  */
real Matrix::gaussian_elimination_forward(real epsilon)
{
    int n = rows;
    Matrix& A = *this;
    int det = 1;
    for (int k=0; k<n-1; ++k) {

        real d = A(k,k);
        int m = k;
        while (fabs(d)<epsilon) {
            m++;
            if (m>=n) {
                A.print(stderr);
                throw std::runtime_error("Could not do the elimination, matrix singular");
            }
            d = A(m,k);
        }
        if (m!=k) {
            for (int j=k; j<n; ++j) {
                real tmp = A(k,j);
                A(k,j) = A(m,j);
                A(m,j) = tmp;
            }
            det *= -1;
        }

        real inv_d = 1.0 / d;

        for (int i=k+1; i<n; ++i) {
            for (int j=k+1; j<n; ++j) {
                A(i, j) -= A(i, k)*A(k, j)*inv_d;
            }
        }
        // for all the subsequent rows
        for (int i=k+1; i<n; ++i) {
            // clear columns up to and including pivot column
            for (int j=k; j<=k; ++j) {
                A(i, j) = 0.0;
            }
        }
    }
    return (real) det;
}

/** For matrix X=AB, find matrices A, B.
  
  The matrices A, B are returned in a std::vector<Matrix>.

  \f\[
  AB = X
  \f\]

  \f\[
  a_{ii} = 1 \for i \in \{1, \ldots, n\}
  \f\]

  \f\[
  b_{ij} = x_{ij} - \sum_{k=1}^{i-1} a_{ik}b_{kj}
  \f\]
 

  */
std::vector<Matrix> Matrix::LUDecomposition(real& determinant, real epsilon)
{
    if (rows!=columns) {
        throw std::domain_error("LU Decomposition cannot be performed for non-square matrices");
    }
    const int n = rows;
    std::vector<Matrix> retval;

    retval.push_back(Matrix::Unity(n,n));
    retval.push_back(Matrix(n,n));
    Matrix& A = retval[0];
    Matrix& B = retval[1];
    Matrix& X = (*this);

    determinant = 1.0;
    for (int k=0; k<n-1; ++k) {
        real d = X(k,k);
        int m = k;
        while (fabs(d)<epsilon) {
            m++;
            if (m>=n) {
                A.print(stderr);
                throw std::runtime_error("Could not do the elimination, matrix singular");
            }
            d = X(m,k);
        }
        if (m!=k) {
            for (int j=k; j<n; ++j) {
                real tmp = X(k,j);
                X(k,j) = X(m,j);
                X(m,j) = tmp;

            }
            determinant *= -1.0;
        }
    }

    // for each j \in [1,N]
    for (int j=0; j<n; ++j) {
        // find b_ij
        for (int i=0; i<=j; ++i) {
            real sum = 0.0;
            for (int k=0; k<i; ++k) {
                sum += A(i,k)*B(k,j);
            }
            B(i,j) = X(i,j) - sum;
        }

        real invd = 1.0/B(j,j);        
        for (int i=j+1; i<n; ++i) {
            real sum = 0.0;
            for (int k=0; k<j; ++k) {
                sum += A(i,k)*B(k,j);
            }
            A(i,j) = invd * (X(i,j) - sum);
        }
    }
    return retval;
}

real Matrix::ColumnSum(int c)
{
    real sum = 0.0;
    for (int i=0; i<rows; i++) {
        sum += (*this)(i,c);
    }
    return sum;
}
real Matrix::RowSum(int r)
{
    real sum = 0.0;
    for (int i=0; i<columns; i++) {
        sum += (*this)(r,i);
    }
    return sum;
}
