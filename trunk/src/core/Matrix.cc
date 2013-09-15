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
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <cstring>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
//#include <gsl/gsl_cblas.h>

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

/// Return a matrix that is the transpose of rhs.
Matrix Transpose(const Matrix& rhs)
{
    Matrix tmp = Matrix(rhs, false);
    tmp.Transpose();
    return tmp;
}

bool Matrix::isSymmetric() const
{
    if (rows != columns)
        return false;

    for (int i=0; i<rows; ++i)  {
        for (int j=i+1; j<columns; ++j) {
            if ((*this)(i,j) != (*this)(j,i)) {
                return false;
            }
        }
    }
    return true;
}

#ifdef REFERENCE_ACCESS
void Matrix::MakeReferences() {
    x_list = (real**) malloc(rows * sizeof(real*));
    for (int i=0; i<rows; ++i) {
        x_list[i] = &x[i*columns];
    }
}
#endif

Matrix::Matrix (int rows_, int columns_,  enum BoundsCheckingStatus check_)
    : rows(rows_),
      columns(columns_),
      checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
    int N = rows*columns;
    x = (real*) calloc(N, sizeof(real));
#ifdef REFERENCE_ACCESS
    MakeReferences();
#endif
}

Matrix::Matrix (int rows_, int columns_, real* y, enum BoundsCheckingStatus check_)
    : rows(rows_),
      columns(columns_),
      x(y),
      checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
#ifdef REFERENCE_ACCESS
    MakeReferences();
#endif

}

Matrix::Matrix (const Vector& v, enum BoundsCheckingStatus check_)
    : checking_bounds(check_),
      transposed(false),
      clear_data(true)
{
    rows = v.Size();
    columns = 1;
    size_t  N = rows * sizeof(real);
    x = (real*) malloc(N);
    memcpy(x, v.x, N);
#ifdef REFERENCE_ACCESS
    MakeReferences();
#endif

}

/// Copy constructor
Matrix::Matrix (const Matrix& rhs, bool clone)
    : rows(rhs.Rows()),
      columns(rhs.Columns()),
      checking_bounds(rhs.checking_bounds),
      transposed(rhs.transposed)
{

    const int M = rows;
    const int N = columns;
    const int K = M*N;

    if (clone) {
        x = (real*) malloc (sizeof(real)*K);
#ifdef REFERENCE_ACCESS
        MakeReferences();
#endif
        for (int m=0; m<M; ++m) {
            for (int n=0; n<N; ++n) {
                (*this)(m,n) = rhs(m,n);
            }
        }
        clear_data = true;
    } else {
        x = rhs.x;
#ifdef REFERENCE_ACCESS
        x_list = rhs.x_list;
#endif
        clear_data = false;
    }

}

Matrix::~Matrix()
{
    if (clear_data) {
        free(x);
#ifdef REFERENCE_ACCESS
        free(x_list);
#endif
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
    assert(x);
#ifdef REFERENCE_ACCESS
    x_list = (real**) realloc(x_list, rows * sizeof(real*));
    for (int i=0; i<rows; ++i) {
        x_list[i] = &x[i*columns];
    }
#endif    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) = rhs(m,n);
        }
    }


    return *this;
}

void Matrix::Clear ()
{
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            (*this)(i,j) = 0.0;
        }
    }
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
    assert(x);

    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            (*this)(m,n) = 0.0;
        }
    }

#ifdef REFERENCE_ACCESS
    x_list = (real**) malloc(rows * sizeof(real*));
    for (int i=0; i<rows; ++i) {
        x_list[i] = &x[i*columns];
    }
#endif

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

/** Create a matrix through the addition of two other matrices.
	
	We distinguish four cases.

	1. \f$C = A + B\f$. Then \f$C = A; C += B\f$
	2. \f$C = A^\top + B\f$. Then \f$C = A^\top; C += B\f$.
	3. \f$C = A + B^\top\f$. Then \f$C = B^\top; C += A\f$.
	4. \f$C = A^\top + B^\top\f$. Then \f$C^\top = A; C^\top = B\f$.

 */
Matrix Matrix::operator+ (const Matrix& rhs)
{
    if (Columns() != rhs.Columns() 
        || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix addition error\n");
    }
    int M = Rows();
    int N = Columns();
    
    Matrix lhs(M, N);
#if 1
	gsl_matrix_const_view A_view
		= gsl_matrix_const_view_array(x, rows, columns);
	gsl_matrix_const_view B_view
		= gsl_matrix_const_view_array(rhs.x, rhs.rows, rhs.columns);
	gsl_matrix_view C_view
		= gsl_matrix_view_array(lhs.x, lhs.rows, lhs.columns);
	
	if (!transposed && !lhs.transposed) {
		gsl_matrix_memcpy(&C_view.matrix, &A_view.matrix);
		gsl_matrix_add(&C_view.matrix, &B_view.matrix);
	} else if (transposed && !lhs.transposed) {
		gsl_matrix_transpose_memcpy(&C_view.matrix, &A_view.matrix);
		gsl_matrix_add(&C_view.matrix, &B_view.matrix);
	} else if (!transposed && lhs.transposed) {
		gsl_matrix_transpose_memcpy(&C_view.matrix, &B_view.matrix);
		gsl_matrix_add(&C_view.matrix, &A_view.matrix);
	} else {
		gsl_matrix_memcpy(&C_view.matrix, &A_view.matrix);
		gsl_matrix_add(&C_view.matrix, &B_view.matrix);
		lhs.transposed = true;
	}
#else
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            lhs(m, n) = (*this)(m,n) + rhs(m,n);
        }
    }
#endif
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
    int N = rhs.Columns();
    
    Matrix C(M, N);
	
	CBLAS_TRANSPOSE Trans_A = transposed ? CblasTrans : CblasNoTrans;
	CBLAS_TRANSPOSE Trans_B = rhs.transposed ? CblasTrans : CblasNoTrans;

	gsl_matrix_const_view A_view
		= gsl_matrix_const_view_array(x, rows, columns);
	gsl_matrix_const_view B_view
		= gsl_matrix_const_view_array(rhs.x, rhs.rows, rhs.columns);
	gsl_matrix_view C_view = gsl_matrix_view_array(C.x, M, N);

	gsl_blas_dgemm(Trans_A, Trans_B,
				   1.0, &A_view.matrix, &B_view.matrix,
				   0.0, &C_view.matrix);
	return C;

#if 0    
    //int K = Columns();
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
#endif
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

/// Divide a matrix by a scalar, creating a new matrix
Matrix Matrix::operator/ (const real& rhs)
{
	real inv = 1.0 / rhs;
	return (*this) * inv;
}

/// Add a matrix to a scalar, creating a new matrix
Matrix Matrix::operator+ (const real& rhs)
{
    Matrix lhs(rows, columns);
    if (transposed) {
        lhs.Transpose();
    }
    real* y = x;
    real* z = lhs.x;
    int N = rows * columns;
    for (int n=0; n<N; ++n) {
        (*z++) = (*y++) + rhs;
    }
    return lhs;
}

/// Subtract a scalar from a matrix, creating a new matrix
Matrix Matrix::operator- (const real& rhs)
{
    Matrix lhs(rows, columns);
    if (transposed) {
        lhs.Transpose();
    }
    real* y = x;
    real* z = lhs.x;
    int N = rows * columns;
    for (int n=0; n<N; ++n) {
        (*z++) = (*y++) - rhs;
    }
    return lhs;
}

void Matrix::print(FILE* f) const
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
///
/// The vector is Kx1, and the rhs is 1xN, by necessity,
/// giving a KxN matrix
Matrix operator* (const Vector& lhs, const Matrix& rhs)
{
    if (rhs.Rows() != 1) {
        throw std::domain_error("Vector-matrix multiplication error\n");
    }
    int K = lhs.Size();
    int N = rhs.Columns();
    Matrix v(K, N);

    for (int i=0; i<K; ++i) {
        for (int j=0; j<N; ++j) {
            v(i,j) = lhs[i] * rhs(0, j);
        }
    }
    return v;
}


/// Multiply a vector with a matrix, creating a new vector
///
/// The lhs is NxK, and the rhs is Kx1, giving a Nx1 Vector
Vector operator* (const Matrix& lhs, const Vector& rhs)
{
    if (rhs.Size() != lhs.Columns()) {
        fprintf (stderr, "Multiplication: (%d x %d) * (%d x 1)\n", lhs.Rows(), lhs.Columns(), rhs.Size());
        throw std::domain_error("matrix-vector multiplication error\n");
    }

    int K = lhs.Rows();
    int N = lhs.Columns();

    Vector v(K);

    for (int i=0; i<K; ++i) {
        real vi = 0;
        for (int j=0; j<N; ++j) {
            vi += lhs(i,j) * rhs(j);
        }
        v(i) = vi;
    }
    return v;
}

real Matrix::det() const
{
    if (isTriangular()) {
        return compute_det_triangular();
    } else {
        Matrix A = *this;
        real swap = A.gaussian_elimination_forward();
        return swap * A.compute_det_triangular();
    }

}

/// Matrix trace.
real Matrix::tr() const
{
	assert(rows == columns);
	if(rows != columns){
		throw std::domain_error("Trace cannot be computed for non-square matrices");
	}
	real trace = 0.0;
	for(int i = 0; i<rows; ++i)
	{
		trace += (*this)(i,i);
	}
	return trace;
}

/// Calculate the element-by-element product of two arrays.
Matrix Matrix::Multiple(const Matrix& rhs) const
{
	if (Columns() != rhs.Columns() || Rows() != rhs.Rows()) {
        throw std::domain_error("Matrix element-by-element product\n");
    }
	int M = Rows();
    int N = Columns();
    
    Matrix lhs(M, N);
    
    for (int m=0; m<M; ++m) {
        for (int n=0; n<N; ++n) {
            lhs(m, n) = (*this)(m,n) * rhs(m,n);
        }
    }
	return lhs;
}

/// Kronecker product.
Matrix Matrix::Kron(const Matrix& rhs) const
{
	int r_rows = rhs.Rows();
	int r_cols = rhs.Columns();
	
	Matrix K(rows*r_rows,columns*r_cols);
	
	for(int i = 0; i < rows; ++i)
	{
		for(int j = 0; j < columns; ++j)
		{
			for(int k = 0; k < r_rows; ++k)
			{
				for(int l = 0; l < r_cols; ++l)
				{
					K(i*r_rows + k, j*r_cols + l) = (*this)(i,j)*rhs(k,l);
				}
			}
		}
	}
	return K;
}

bool Matrix::isTriangular() const
{
    if (rows != columns) {
        throw std::domain_error("only square matrices can be triangular");
    }
    if (isLowerTriangular() || isUpperTriangular()) {
        return true;
    }
    return false;
}

bool Matrix::isUpperTriangular() const
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

bool Matrix::isLowerTriangular() const
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

real Matrix::compute_det_triangular() const
{
    const Matrix& A = *this;
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

/** For matrix X=AB, find matrices A, B, using in-place elimination.
  
  The matrices A, B are returned in a std::vector<Matrix>.

  \f[
  AB = X
  \f]

  \f[
  a_{ii} = 1, \quad \textrm{for} ~ i \in \{1, \ldots, n\}
  \f]

  \f[
  b_{ij} = x_{ij} - \sum_{k=1}^{i-1} a_{ik}b_{kj}
  \f]
 

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


    // pivot non-zero elements
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

    determinant = 1.0;
    for (int i=0; i<n; ++i) {
        determinant *= A(i,i)*B(i,i);
    }
    return retval;
}

/// QR decomposition using Householder reflections.
std::vector<Matrix> Matrix::QRDecomposition() const
{
	
	const int m = rows;
	const int n = columns;
	std::vector<Matrix> retval;

	retval.push_back(Matrix(m,n));
    retval.push_back(Matrix(n,n));
    Matrix& Q = retval[0];
    Matrix& R = retval[1];
	Vector d(n);
	Matrix  A = (*this);
	real f;
	
	for(int k = 0; k < n; ++k)
	{
		real s = 0.0;
		for(int j = k; j < m; ++j) {
			s += A(j,k)*A(j,k);
		}
		s = sqrt(s);
		if(s == 0){
			fprintf(stderr,"\nERROR: Matrix rank is smaller than the number of columns\n");
			throw std::runtime_error("Could not do QR Decomposition, matrix not positive definite");
		}
		d(k) = (A(k,k) > 0) ? (-s) : s;
		f = sqrt(s*(s + fabs(A(k,k))));
		A(k,k) = A(k,k) - d(k);
		for(int j = k; j < m; ++j) {
			A(j,k) = A(j,k) / f;
		}
		for(int i = (k+1); i < n; ++i)
		{
			s = 0.0;
			for(int j = k; j < m; ++j) {
				s += A(j,k)*A(j,i);
			}
			for(int j = k; j < m; ++j) {
				A(j,i) = A(j,i) - A(j,k)*s;
			}
		}
	}
	for (int i = 0;	i < n; ++i) {
		Vector y(m);
		y(i) = 1.0;
		for( int j = (n-1); j >= 0; --j)
		{
			real s = 0.0;
			for( int k = j; k < m; ++k) {
				s += A(k,j)*y(k);
			}
			for( int k = j; k < m; ++k) {
				y(k) = y(k) - A(k,j)*s;
			}
		}
		Q.setColumn(i,y);
	}
	for(int i=0; i<n; ++i) {
		R(i,i) = d(i);
		for(int j = (i + 1); j<n; ++j) {
			R(i,j) = A(i,j);
		}
	}
	return retval;
}	

/** Cholesky decomposition.
    
   
    Do a Cholesky decomposition, after adding epsilon * I to the
    matrix. The original matrix is unchanged.
*/
Matrix Matrix::Cholesky(real epsilon) const
{
    int n = Rows();
    assert (n == Columns());
    Matrix chol(n, n);
    Cholesky(chol, epsilon);
    return chol;
}


/** Cholesky decomposition.
    
    Based on code by Catherine Bourgain (2006).
    Original code at:
    http://www.stat.uchicago.edu/~mcpeek/software/CCQLSpackage1.3/

    Can be safely called with chol = *this; however then the lower triangular
    part of the matrix must be cleared by the user. (TODO?).
*/
void Matrix::Cholesky(Matrix& chol, real epsilon) const
{
    int n = Rows();
    assert (n == Columns());
    for (int i=0; i<n; i++) {
        // do diagonal first
        chol(i,i) = (*this)(i,i) + epsilon;
        for (int k=0; k<i; k++)  {
            chol(i, i) -= chol(k, i)*chol(k, i);
        }
        if (chol(i, i) <= 0) {
            fprintf(stderr,"\nERROR: non-positive definite matrix!\n");
            fprintf(stderr, "\nproblem from %d %f\n",i,chol(i, i));
            throw std::runtime_error("Could not do Cholesky, matrix not positive definite");
        }
        chol(i, i) = sqrt(chol(i, i));
        
        
        for (int j=i+1; j<n; j++) {
            chol(i, j) = (*this)(i, j);
            for (int k=0; k<i; k++)
                chol(i, j) -= chol(k, i)*chol(k, j);
            chol(i, j) /= chol(i, i);
        }
    }
}

/** Invert matrix using GSL LU Decomp */
Matrix Matrix::GSL_Inverse() const
{
	int N = Rows();
	assert(N==Columns());
	Matrix A(*this);
	Matrix R(N, N);
	R.transposed = transposed;
	gsl_matrix_view M_view = gsl_matrix_view_array(A.x, N, N);
	gsl_matrix_view R_view = gsl_matrix_view_array(R.x, N, N);
	gsl_permutation * perm = gsl_permutation_alloc (N);
	int s;
	gsl_linalg_LU_decomp (&M_view.matrix, perm, &s);
	gsl_linalg_LU_invert (&M_view.matrix, perm, &R_view.matrix);
	gsl_permutation_free(perm);
	return R;
}
/** Matrix inversion using a factorised matrix A = LU.
    
    Here L is lower triangluar and U is an upper triangular matrix.

    Given \f$A = LU\f$,  solve for \f$X = A^{-1}\f$: 
    \f[
    LUX = I,
    \f]
    by dynamic programming.
 */
Matrix Matrix::Inverse(const Matrix& L, const Matrix& U) const
{
    int n = L.Rows();
    assert(L.Columns() == n);
    assert(U.Rows() == n);
    assert(U.Columns() == n);
    
    // Build this column by column
    Matrix B = Unity(n,n);


    // Forward substitution - find Y : LY=B;
    Matrix Y(n, n);
    for (int c=0; c<n; ++c) {
        Y(0,c) = B(0, c) / L(0,0);
        for (int i=1; i<n; ++i) {
            real s = 0;
            for (int j=0; j<i; ++j) {
                s += Y(j,c) * L(i,j);
            }
            Y(i,c) = (B(i, c) - s) / L(i,i);
        }
    }

    // Backward substitution
    Matrix X(n, n);
    for (int c=0; c<n; ++c) {
        X(n - 1, c) = Y(n - 1, c) / U(n - 1, n - 1);
        for (int i = n-2; i>=0; --i) {
            real s = 0;
            for (int j=i+1; j<n; ++j) {
                s += X(j, c) * U (i, j);
            }
            X(i, c) = (Y(i, c) - s) / U(i, i);
        }
    }

    return X;
}

/// Forward substitution - find Y : LY=B;
/// 
/// Implementation for square matrices only.
Matrix forward_substitution(const Matrix& L, const Matrix& B)
{
    int n = L.Rows();

    Matrix Y(n, n);
    for (int c=0; c<n; ++c) {
        Y(0,c) = B(0, c) / L(0,0);
        for (int i=1; i<n; ++i) {
            real s = 0;
            for (int j=0; j<i; ++j) {
                s += Y(j,c) * L(i,j);
            }
            Y(i,c) = (B(i, c) - s) / L(i,i);
        }
    }
    return Y;
}

/// Get the sum of column c
real Matrix::ColumnSum(int c) const
{
    real sum = 0.0;
    for (int i=0; i<rows; i++) {
        sum += (*this)(i,c);
    }
    return sum;
}

/// Get the sum of column r
real Matrix::RowSum(int r) const
{
    real sum = 0.0;
    for (int i=0; i<columns; i++) {
        sum += (*this)(r,i);
    }
    return sum;
}




/// Return a column vector that is the column-wise maximum
Vector Matrix::ColumnMax() const
{
	Vector C(Rows());
	for (int i=0; i<Rows(); ++i) {
		C(i) = (*this)(i, 0);
		for (int j=1; j<Columns(); ++j) {
			real x = (*this)(i, j);
			if (x > C(j)) {
				C(j) = x;
			}
		}
	}
    return C;
}

void Matrix::Vec(const Vector& x)
{
	assert(!(x.Size()%Rows()));
	for(int i = 0; i<x.Size(); ++i)
	{
		div_t divresult = div(i,Rows());
		(*this)(divresult.rem,divresult.quot) = x(i);
	}
}

/// Return a column vector with the matrix columns stacked.
Vector Matrix::Vec() const
{
	Vector R(Columns()*Rows());
	for(int j=0; j<Columns(); ++j)
	{
		for(int i=0; i<Rows(); ++i)
		{
			R(i + j*Rows()) = (*this)(i,j);
		}
	}
	return R;
}

/// Return a row vector that is the row-wise maximum
Vector Matrix::RowMax() const
{
	Vector R(Columns());
    for (int j=0; j<Columns(); ++j) {
		R(j) = (*this)(0,j);
		for (int i=1; i<Rows(); ++i) {
			real x = (*this)(i, j);
			if (x > R(j)) {
				R(j) = x;
			}
		}
	}
    return R;
}


/// Return a column vector that is the column-wise minimum
Vector Matrix::ColumnMin() const
{
	Vector C(Rows());
	for (int i=0; i<Rows(); ++i) {
		C(i) = (*this)(i, 0);
		for (int j=1; j<Columns(); ++j) {
			real x = (*this)(i, j);
			if (x < C(j)) {
				C(j) = x;
			}
		}
	}
    return C;
}

/// Return a row vector that is the row-wise minimum
Vector Matrix::RowMin() const
{
	Vector R(Columns());
    for (int j=0; j<Columns(); ++j) {
		R(j) = (*this)(0,j);
		for (int i=1; i<Rows(); ++i) {
			real x = (*this)(i, j);
			if (x < R(j)) {
				R(j) = x;
			}
		}
	}
    return R;
}


Vector Matrix::getColumn(int c) const
{
    Vector column(rows);
    for (int i=0; i<rows; i++) {
        column[i] = (*this)(i,c);
    }
    return column;
}

Vector Matrix::getRow(int r) const
{
    Vector row(columns);
    for (int i=0; i<columns; i++) {
        row[i] = (*this)(r,i);
    }
    return row;
}


void Matrix::setColumn(int c, const Vector& x)
{
    assert(rows == x.Size());
    for (int i=0; i<rows; i++) {
        (*this)(i,c) = x[i];
    }
}

void Matrix::setRow(int r, const Vector& x)
{
    assert(columns == x.Size());
    for (int i=0; i<columns; i++) {
        (*this)(r,i) = x[i];
    }

}


void Matrix::SortRow(int r)
{
    Vector row(columns);
    for (int i=0; i<columns; i++) {
        row[i] = (*this)(r,i);
    }
    std::sort(&row[0], &row[columns]);
    for (int i=0; i<columns; i++) {
        (*this)(r,i) = row[i];
    }
}


void Matrix::SortColumn(int c)
{
    Vector column(rows);
    for (int i=0; i<rows; i++) {
        column[i] = (*this)(i,c);
    }

    std::sort(&column[0], &column[rows]);
    for (int i=0; i<rows; i++) {
        (*this)(i,c) = column[i];
    }
}

real Matrix::L1Norm() const
{
    real s = 0.0;
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            s += abs((*this)(i, j));
        }
    }                                           
    return s;
}


real Matrix::L2Norm() const
{
    real s = 0.0;
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            real d = (*this)(i, j);
            s += d * d;
        }
    }
    return s;
}

/// Quick estimation of \f$d = x'Ayf\f$.
///
/// This is much faster than writing d * x * A *  y explicitly.
real Mahalanobis2 (const Vector& x, const Matrix& A, const Vector& y)
{
    assert(x.Size() == y.Size());
    assert(A.Rows() == A.Columns());
    assert(A.Rows() == x.Size());

    int n = x.Size();
    real d = 0;
    for (int i=0; i<n; ++i) {
        real x_i = x(i);
        for (int j=0; j<n; ++j) {
            d += x_i * A(i,j) * y(j);
        }
    }
    return d;
}

