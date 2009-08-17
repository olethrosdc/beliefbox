// -*- Mode: C++ -*-
// $Id: Vector.c,v 1.1 2006/11/07 18:03:35 cdimitrakakis Exp cdimitrakakis $

#include "Vector.h"
#include "SmartAssert.h"
#include "Matrix.h"

#include <exception>
#include <stdexcept>

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
    SMART_ASSERT(lhs && rhs)(lhs)(rhs);
    SMART_ASSERT(n == rhs->n)(lhs->n)(rhs->n);
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
        x = (real*) malloc(sizeof(real)*n);
        for (int i=0; i<n; i++) {
            x[i] = 0.0;
        }
    }
    checking_bounds = check;
}

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

/// Copy constructor
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

Vector::~Vector()
{
    if (x) {
        free(x);
    }
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

/// Sum of all in vector
real Vector::Sum() const
{
    real sum = 0;
    for (int i=0; i<n; ++i) {
	sum += x[i];
    }
    return sum;
}

/// Addition
Vector Vector::operator+ (const Vector& rhs)
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
Vector Vector::operator- (const Vector& rhs)
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
Vector Vector::operator* (const Vector& rhs)
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
Vector Vector::operator/ (const Vector& rhs)
{
    assert (rhs.n==n);
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] * rhs[i];
    }
    return lhs;
}

/// Per-element self-multiplication
Vector& Vector::operator/= (const Vector& rhs)
{
    assert (rhs.n==n);
    for (int i=0; i<n; i++) {
        x[i] *= rhs[i];
    }
    return *this;
}


/* ----------- SCALAR OPERATORS -------------------------*/
/// Scalar addition
Vector Vector::operator+ (const real& rhs)
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] + rhs;
    }
    return lhs;
}
/// Scalar subtraction
Vector Vector::operator- (const real& rhs)
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i] - rhs;
    }
    return lhs;
}
/// Scalar multiplication
Vector Vector::operator* (const real& rhs)
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i]*rhs;
    }
    return lhs;
}
/// Scalar division
Vector Vector::operator/ (const real& rhs)
{
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = x[i]/rhs;
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
    for (int i=0; i<n; i++) {
        x[i] /= rhs;
    }
    return *this;
}


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
void Vector::Resize(int N_)
{ 
    n = N_;
    if (n>maxN) {
        if (n==0) {
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
    SMART_ASSERT(lhs && rhs  && res);
    int n=lhs->n;
    SMART_ASSERT((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] + rhs->x[i];
    }
}  

/// \brief Sub vector rhs from lhs, save result to res.
/// It is safe to use identical handles for all arguments.
void Sub (const Vector* lhs, const Vector* rhs, Vector* res)
{
    SMART_ASSERT(lhs && rhs  && res);
    int n=lhs->n;
    SMART_ASSERT((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] - rhs->x[i];
    }
}  


/// \brief Multiply lhs and rhs element by element, save result to res.
/// It is safe to use identical handles for all arguments.
void Mul (const Vector* lhs, const Vector* rhs, Vector* res)
{
    SMART_ASSERT(lhs && rhs  && res);
    int n=lhs->n;
    SMART_ASSERT((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] * rhs->x[i];
    }
}  

/// \brief Divide lhs and rhs element by element, save result to res.
/// It is safe to use identical handles for all arguments.
void Div (const Vector* lhs, const Vector* rhs, Vector* res)
{
    SMART_ASSERT(lhs && rhs  && res);
    int n=lhs->n;
    SMART_ASSERT((n==rhs->n)&&(n==res->n));
    for (int i=0; i<n; i++) {
        res->x[i] = lhs->x[i] / rhs->x[i];
    }
}  

/// \brief Return the inner product of two vectors.
/// If all vectors are column vectors, this is \f$x=a' b\f$.
real Product (const Vector* const lhs, const Vector* const rhs)
{
    SMART_ASSERT(lhs && rhs);
    int n=lhs->n;
    SMART_ASSERT(n==rhs->n);
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
    SMART_ASSERT(lhs && rhs  && res);
    int n=lhs->n;
    SMART_ASSERT((n==rhs->n)&&(n==res->Columns()))(n)(rhs->n)(res->Columns());
    SMART_ASSERT (n==res->Rows())(n)(res->Rows());
    for (int i=0; i<n; i++) { // columns
        for (int j=0; j<n; j++) { //rows
            (*res)(i,j) = lhs->x[j] * rhs->x[i];
        }
    }
}

#if 0
/// Return \f$\|a-b\|_1\f$, the L1 norm between two vectors.
real L1Norm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return L1Norm (lhs->x, rhs->x, lhs->n);
}

/// Return \f$\|a-b\|\f$, the euclidean norm between two vectors.
real EuclideanNorm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return EuclideanNorm (lhs->x, rhs->x, lhs->n);
}

/// Return \f$\|a-b\|^2\f$, the square of the euclidean norm between two
/// vectors.
real SquareNorm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return SquareNorm (lhs->x, rhs->x, lhs->n);
}

/// Get maximum element
real Max(const Vector* v)
{
    return Max(v->Size(), v->x);
}

/// Get minimum element
real Min(const Vector* v)
{
    return Min(v->Size(), v->x);
}

/// Get maximum element
int ArgMax(const Vector* v)
{
    return ArgMax(v->Size(), v->x);
}

/// Get minimum element
int ArgMin(const Vector* v)
{
    return ArgMin(v->Size(), v->x);
}

real Span(const Vector* v)
{
    return Max(v) - Min(v);
}
#endif
