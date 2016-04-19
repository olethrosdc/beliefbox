#include "svgp.h"
#include "real.h"
#include "math.h"
#include <random>
#include "MultivariateNormal.h"
#include <gsl/gsl_linalg.h>

SVGP::SVGP(Matrix& X_, Vector& Y_, Matrix& Z_, real noise_var_, real sig_var_, Vector scale_length_)
	: X(X_),
	  Y(Y_),
	  Z(Z_),
	  noise_var(noise_var_),
	  sig_var(sig_var_),
	  scale_length(scale_length_)
{
    init();
}
SVGP::SVGP(Matrix &X_, Vector& Y_, int num_inducing, real noise_var_, real sig_var_, Vector scale_length_)
    : X(X_),
      Y(Y_),
      noise_var(noise_var_),
      sig_var(sig_var_),
      scale_length(scale_length_)
{
    //create Z with inducing points placed decided by k-means on X

    srand(time(NULL));
    Z = Matrix(num_inducing,X.Columns());
    for(int i=0;i<num_inducing;i++) {
        Z.setRow(i,rand()%X.Rows());
    }
    //std::vector<std::vector<Vector>> kmeans(num_inducing);
    //the assignment + update step here...
    init();
}
void SVGP::init() {
    Beta = 1/(noise_var*noise_var);
    N = X.Rows();
    d = X.Columns();
    m = Z.Rows();
    assert(d==Z.Columns());
    assert(N==Y.Size());
    assert(d==scale_length.Size());
    FullUpdateGaussianProcess();
}
Matrix SVGP::Kernel(const Matrix& A, const Matrix& B) {
	assert(A.Columns()==B.Columns());

	int N_ = A.Rows();
	int M_ = B.Rows();

	Matrix cov(N_,M_);

	real extra_noise = 0;
	if(N_ == M_)
		extra_noise = noise_var*noise_var;

	real sig_noise = sig_var*sig_var;
	for(int i=0; i<M_; ++i) {
		Vector x = B.getRow(i);
		for(int j=0; j<N_; j++) {
			real delta = ((x - A.getRow(j))/scale_length).SquareNorm();
			delta = sig_noise*exp(-0.5*delta);
			if(i == j) 
				delta += extra_noise;
			cov(j,i) = delta;
			//cov(i,j) = delta;
		}
	}
	return cov;
}

void SVGP::UpdateGaussianProcess() {
	Kmm = Kernel(Z,Z);
	Knm = Kernel(X,Z);
	Kmn = Kernel(Z,X);

	invKmm = Kmm.Inverse();
	K_tilde = Knn - Knm * invKmm * Kmn;

	u = MultivariateNormal(m,S).generate();
}
void SVGP::FullUpdateGaussianProcess() {
	Kmm = Kernel(Z,Z);
	Knm = Kernel(X,Z);
	Kmn = Kernel(Z,X);

	//expensive stuff

	Knn = Kernel(X,X);

	invKmm = Kmm.Inverse();
	K_tilde = Knn - Knm * invKmm * Kmn;

	Matrix q_prec = Beta * invKmm * Kmn * Knm * invKmm + invKmm;
	Matrix q_var = q_prec.Inverse();
	Matrix q_mean_tmp = Beta * q_prec.Inverse() * invKmm * Kmn;
	Vector q_mean = q_mean_tmp * Y;

	S = q_var;
	m = q_mean;

	u = MultivariateNormal(m,S).generate();
}
void SVGP::Prediction(const Vector& x, real& mean, real& var) {
	Matrix tmp(x);
	Vector Kx = Kernel(tmp,Z).getRow(0);

	mean = Product(Kx,invKmm * u);
	var = sig_var*sig_var - Product(Kx,invKmm * Kx);
}
real SVGP::LogLikelihood() {
	//D_KL(q(u)||p(u))

	Matrix q_prec = Beta * invKmm * Kmn * Knm * invKmm + invKmm;
	Matrix q_var = q_prec.Inverse();

	Matrix q_mean_tmp = Beta * q_prec.Inverse() * invKmm * Kmn;
	Vector q_mean = q_mean_tmp * Y;

	p_var = Kmm;
	p_mean = Vector::Null(q_mean.Size());

	Matrix inv_p_var = p_var.Inverse();
	real lhs = (inv_p_var * q_prec).tr();
	Vector A = inv_p_var * (p_mean - q_mean);
	real B = Product((p_mean - q_mean),A);
	KL = lhs + B - N + log(p_var.det()/q_prec.det());
	KL *= 0.5;

	L = 0;
	for(int i=0;i<N;i++) {
		Vector k_i = Kmn.getColumn(i);
		Matrix outer = OuterProduct(k_i,k_i);

		Matrix lambda_i = Beta * invKmm * outer * invKmm;
		real k_ii = K_tilde(i,i);
		Vector rhs = invKmm * q_mean;
		real pdfmean = Product(k_i,rhs);
		real pdfvar = 1/Beta;
		real pdfcond = Y(i);

		real pdf = (real) 1/(pdfvar * sqrt(2 * M_PI));
		real expon = exp(-(pdfcond-pdfmean)*(pdfcond-pdfmean)/(2*pdfvar*pdfvar));
		pdf *= expon;
		pdf = (real) log(pdf);
		
		L += pdf - 1/2 * Beta * k_ii - 1/2 * (S * lambda_i).tr();
	}
	L -= KL;
	return L;
}
//OK, here is where I'm a bit confused on how to update the inducing inputs Z
void SVGP::local_update(const Matrix& X_samples, const Vector& Y_samples) {
    //here Z should be updated somehow..
    //first, transform Z matrix into double array

    double *data;
}
void SVGP::global_update(const Matrix& X_samples, const Vector& Y_samples) {

	Matrix Kmn_samples = Kernel(Z,X_samples);
	Matrix Knm_samples = Kernel(X_samples,Z);

	Matrix q_prec = Beta * invKmm * Kmn_samples * Knm_samples * invKmm + invKmm;

	Vector theta_1 = S.Inverse() * m;
	Matrix theta_2 = -0.5 * S.Inverse();

	Vector tmpResult = Beta * invKmm * Kmn_samples * Y_samples - theta_1;
	tmpResult *= l;
	tmpResult += theta_1;
	theta_1 = tmpResult;
	theta_2 = theta_2 + l * (-1/2 * q_prec + 1/2 * theta_2);


	S = theta_2.Inverse();
	S *= -2.0;

	m = S * theta_1;
}

Vector SVGP::optimize_Z(int max_iters) {


	//eval loglikelihood

}
void SVGP::AddObservation(const std::vector<Vector>& x, const std::vector<real>& y) {
	//assert(x[0].Size()==d);
	for(int i = 0;i<x.size();i++) {
		Vector x_ = x[i];
		real y_ = y[i];
		AddObservation(x_,y_);
	}
	FullUpdateGaussianProcess();
}
void SVGP::AddObservation(const Vector& x, const real& y) {
	X.AddRow(x);
	Y.AddElement(y);
}
void SVGP::Clear() {
	//TODO
}
