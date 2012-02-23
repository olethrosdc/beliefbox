// -*- Mode: C++ -*-
// $Id: Vector.c,v 1.1 2006/11/07 18:03:35 cdimitrakakis Exp cdimitrakakis $

#include "Vector.h"
#include "Matrix.h"

#include <exception>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cassert>

Vector Vector::Unity(int N_, enum BoundsCheckingStatus check)
{
    Vector v(N_);
    for (int i=0; i<N_; ++i) {
        v(i) = 1.0;
    }
    return v;
}


Vector* NewVector (int n)
{
    Vector* v;
    v = new Vector (n);
    return v;
}

int DeleteVector (Vector* vector)
{
    delete vector;
    return 0;
}

void CopyVector (Vector* const lhs, const Vector * const rhs)
{
    int n = lhs->n;
    assert(lhs && rhs);
    assert(n == rhs->n);
    for (int i=0; i<n; i++) {
        lhs->x[i] = rhs->x[i];
    }
}


#if 1
Vector::Vector()
{
    n = 0;
    maxN = 0;
    x = NULL;
    checking_bounds = NO_CHECK_BOUNDS;
}

Vector::Vector(int N_, enum BoundsCheckingStatus check)
{
    n = N_;
    maxN = n;
    if (n==0) {
        x = NULL;
    } else {
        x = (real*) calloc(n, sizeof(real));
    }
    checking_bounds = check;
}

/** Copy from an array.
    TODO replace assignment with memcopy
*/
Vector::Vector (int N_, real* y, enum BoundsCheckingStatus check)
{
    n = N_;
    maxN = n;
    if (n==0) {
        x = NULL;
    } else {
        x = (real*) malloc(sizeof(real)*n);
        for (int i=0; i<n; i++) {
            x[i] = y[i];
        }
    }
    checking_bounds = check;
}

/** Copy constructor
    TODO replace assignment with memcopy
*/
Vector::Vector (const Vector& rhs)
{
    n = rhs.n;
    maxN = n;
    if (n==0) {
        x = NULL;
    } else {
        x = (real*) malloc(sizeof(real)*n);
        for (int i=0; i<n; i++) {
            x[i] = rhs[i];
        }
    }
    checking_bounds = rhs.checking_bounds;
}

/// Destructor
Vector::~Vector()
{
    if (x) {
        free(x);
    }
}

/// Print vector out at file f, make a new line
void Vector::print(FILE* f) const
{
    printf(f);
    fprintf(f, "\n");
}

/// Print a vector out without newline.
void Vector::printf(FILE* f) const
{
    for (int i=0; i<n; ++i) {
        fprintf (f, "%f ", x[i]);
    }
}


/// Inequality operator.
///
/// Return true iff this[i] > rhs[i] for all i.
const bool Vector::operator> (const Vector& rhs) const
{
    if (this == &rhs) return false;
	assert(n == rhs.n);
    for (int i=0; i<n; i++) {
		if (x[i] <= rhs[i]) {
			return false;
		}
    }
    return true;
}

/// Inequality operator.
///
/// Return true iff this[i] < rhs[i] for all i.
const bool Vector::operator< (const Vector& rhs) const
{
    if (this == &rhs) return false;
	assert(n == rhs.n);
    for (int i=0; i<n; i++) {
		if (x[i] >= rhs[i]) {
			return false;
		}
    }
    return true;
}

/// Equality operator.
///
/// Return true iff this[i] == rhs[i] for all i.
const bool Vector::operator== (const Vector& rhs) const
{
    if (this == &rhs) return true;
	assert(n == rhs.n);
    for (int i=0; i<n; i++) {
		if (x[i] != rhs[i]) {
			return false;
		}
    }
    return true;
}

/// Assignment operator
Vector& Vector::operator= (const Vector& rhs)
{
    if (this == &rhs) return *this;
    Resize(rhs.n);
    for (int i=0; i<n; i++) {
        x[i] = rhs[i];
    }
    return *this;
}

/// Assignment operator
Vector& Vector::operator= (const real& rhs)
{
    Resize(1);
    x[0] = rhs;
    return *this;
}

/// log Sum of all in vector
real Vector::logSum() const
{
    real log_sum = LOG_ZERO;
    for (int i=0; i<n; ++i) {
        log_sum = logAdd (log_sum, x[i]);
    }
    return log_sum;
}

/// Sum of all in vector
real Vector::Sum() const
{
    real sum = 0;
    for (int i=0; i<n; ++i) {
        sum += x[i];
    }
    return sum;
}

/// P-norm of a vector
real Vector::Norm(real p) const
{
    assert(p > 0);
    real log_sum = LOG_ZERO;
    for (int i=0; i<n; ++i) {
        log_sum = logAdd(log_sum, fabs(x[i]) * p);
    }
    return exp(log_sum / p);
}

/// L1norm of a vector
real Vector::L1Norm() const
{
    real sum = 0.0;
    for (int i=0; i<n; ++i) {
        sum += fabs(x[i]);
    }
    return sum;
}

/// L2norm of a vector
real Vector::L2Norm() const
{
    return sqrt(SquareNorm());
}

/// L2norm of a vector
real Vector::SquareNorm() const
{
    real sum = 0.0;
    for (int i=0; i<n; ++i) {
        sum += x[i] * x[i];
    }
    return sum;
}

/// Sum a range of a vector
real Vector::Sum(const int start, const int end) const
{
    if (checking_bounds) {
        if (start < 0 || end >= n) {
            throw std::out_of_range("index out of range");
        }
    }
    real sum = 0;
    for (int i=start; i<=end; ++i) {
        sum += x[i];
    }
    return sum;
}

/// Addition
const Vector Vector::operator+ (const Vector& rhs) const
{
    assert (rhs.n==n);
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] + rhs[i];
    }
    return lhs;
}

/// self-addition
Vector& Vector::operator+= (const Vector& rhs)
{
    assert (rhs.n==n);
    for (int i=0; i<n; i++) {
        x[i] += rhs[i];
    }
    return *this;
}


/// Subtraction
const Vector Vector::operator- (const Vector& rhs) const
{
    assert (rhs.n==n);
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] - rhs[i];
    }
    return lhs;
}

/// Self-substraction
Vector& Vector::operator-= (const Vector& rhs)
{
    assert (rhs.n==n);
    for (int i=0; i<n; i++) {
        x[i] -= rhs[i];
    }
    return *this;
}


/// Per-element multiplication
const Vector Vector::operator* (const Vector& rhs) const
{
    assert (rhs.n==n);
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] * rhs[i];
    }
    return lhs;
}

/// Per-element self-multiplication
Vector& Vector::operator*= (const Vector& rhs)
{
    assert (rhs.n==n);
    for (int i=0; i<n; i++) {
        x[i] *= rhs[i];
    }
    return *this;
}

/// Per-element division
const Vector Vector::operator/ (const Vector& rhs) const
{
    assert (rhs.n==n);
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] / rhs[i];
    }
    return lhs;
}

/// Per-element self-division
Vector& Vector::operator/= (const Vector& rhs)
{
    assert (rhs.n==n);
    for (int i=0; i<n; i++) {
        x[i] /= rhs[i];
    }
    return *this;
}


/* ----------- SCALAR OPERATORS -------------------------*/
/// Scalar addition
const Vector Vector::operator+ (const real& rhs) const
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] + rhs;
    }
    return lhs;
}
/// Scalar subtraction
const Vector Vector::operator- (const real& rhs) const
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] - rhs;
    }
    return lhs;
}
/// Scalar multiplication
const Vector Vector::operator* (const real& rhs) const
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i]*rhs;
    }
    return lhs;
}
/// Scalar division
const Vector Vector::operator/ (const real& rhs) const
{
    Vector lhs (n);
    real inv = 1.0 / rhs;
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] * inv;
    }
    return lhs;
}


/* ----------- SELF SCALAR OPERATORS -------------------------*/

/// Self scalar addition
Vector& Vector::operator+= (const real& rhs)
{
    for (int i=0; i<n; i++) {
        x[i] += rhs;
    }
    return *this;
}

/// Self scalar subtraction
Vector& Vector::operator-= (const real& rhs)
{
    for (int i=0; i<n; i++) {
        x[i] -= rhs;
    }
    return *this;
}


/// Self scalar multiplication
Vector& Vector::operator*= (const real& rhs)
{
    for (int i=0; i<n; i++) {
        x[i] *= rhs;
    }
    return *this;
}

/// Self scalar multiplication
Vector& Vector::operator/= (const real& rhs)
{
    real inv = 1.0 / rhs;
    for (int i=0; i<n; i++) {
        x[i] *= inv;
    }
    return *this;
}

/// Exponentiation


/* --------------- INDEXING ------------------- */
/// index
real& Vector::operator[] (int index)
{
    if (checking_bounds) {
        if ((index<0) || (index>=n)) {
            throw std::out_of_range("index out of range");
        }
    }
    return x[index];
}

/// index
const real& Vector::operator[] (int index) const
{
    if (checking_bounds) {
        if ((index<0) || (index>=n)) {
            throw std::out_of_range("index out of range");
        }
    }
    return x[index];
}

/// Change size
void Vector::Clear()
{ 
    memset(x, 0, sizeof(real) * n);
}

/// Change size
void Vector::Resize(int N_)
{ 
    n = N_;
    if (n > maxN) {
        if (maxN == 0) {
            x = (real*) malloc (sizeof(real)*n);
        } else {
            x = (real*) realloc(x, sizeof(real)*n);
        }
        maxN = n;
    }
}
#endif

/// \brief Add vector lhs to rhs, save result to res.
/// It is safe to use identical handles for all arguments.
void Add (const Vector* lhs, const Vector* rhs, Vector* res)
{
    assert(lhs && rhs  && res);
    int n=lhs->n;
    assert((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] + rhs->x[i];
    }
}  

/// \brief Sub vector rhs from lhs, save result to res.
/// It is safe to use identical handles for all arguments.
void Sub (const Vector* lhs, const Vector* rhs, Vector* res)
{
    assert(lhs && rhs  && res);
    int n=lhs->n;
    assert((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] - rhs->x[i];
    }
}  


/// \brief Multiply lhs and rhs element by element, save result to res.
/// It is safe to use identical handles for all arguments.
void Mul (const Vector* lhs, const Vector* rhs, Vector* res)
{
    assert(lhs && rhs  && res);
    int n=lhs->n;
    assert((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] * rhs->x[i];
    }
}  

/// \brief Divide lhs and rhs element by element, save result to res.
/// It is safe to use identical handles for all arguments.
void Div (const Vector* lhs, const Vector* rhs, Vector* res)
{
    assert(lhs && rhs  && res);
    int n=lhs->n;
    assert((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] / rhs->x[i];
    }
}  

/// \brief Return the inner product of two vectors.
/// If all vectors are column vectors, this is \f$x=a' b\f$.
real Product (const Vector* const lhs, const Vector* const rhs)
{
    assert(lhs && rhs);
    int n=lhs->n;
    assert(n==rhs->n);
    real res = 0.0;
    for (int i=0; i<n; i++) {
        res += lhs->x[i] * rhs->x[i];
    }
    return res;
}


/// \brief Save product of two vectors onto a matrix.
/// If all vectors are column vectors, this is \f$x=a b'\f$.
void Product (const Vector* lhs, const Vector* rhs, Matrix* res)
{
    assert(lhs && rhs  && res);
    int n=lhs->n;
    assert((n==rhs->n)&&(n==res->Columns()));
    assert (n==res->Rows());
    for (int i=0; i<n; i++) { // columns
        for (int j=0; j<n; j++) { //rows
            (*res)(i,j) = lhs->x[j] * rhs->x[i];
        }
    }
}

/// \brief Return the inner product of two vectors.
/// If all vectors are column vectors, this is \f$x=a' b\f$.
real Product (const Vector& lhs, const Vector& rhs)
{
    int n=lhs.n;
    assert(n==rhs.n);
    real res = 0.0;
    for (int i=0; i<n; i++) {
        res += lhs.x[i] * rhs.x[i];
    }
    return res;
}


/// \brief Save product of two vectors onto a matrix.
/// If all vectors are column vectors, this is \f$x=a b'\f$.
void Product (const Vector& lhs, const Vector& rhs, Matrix& res)
{
    int n=lhs.n;
    assert((n==rhs.n)&&(n==res.Columns()));
    assert (n==res.Rows());
    for (int i=0; i<n; i++) { // columns
        for (int j=0; j<n; j++) { //rows
            res(i,j) = lhs.x[j] * rhs.x[i];
        }
    }
}



