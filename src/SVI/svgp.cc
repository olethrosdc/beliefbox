#include "svgp.h"
#include "real.h"
#include "math.h"
#include <random>
#include "MultivariateNormal.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

SVGP::SVGP(Matrix& X_, Vector& Y_, Matrix& Z_, real noise_var_, real sig_var_, Vector scale_length_, int samples)
	: X(X_),
	  Y(Y_),
	  Z(Z_),
	  noise_var(noise_var_),
	  sig_var(sig_var_),
	  scale_length(scale_length_),
	  samples(samples)
{
    init();
}
SVGP::SVGP(Matrix &X_, Vector& Y_, int num_inducing, real noise_var_, real sig_var_, Vector scale_length_, int samples)
    : X(X_),
      Y(Y_),
      noise_var(noise_var_),
      sig_var(sig_var_),
      scale_length(scale_length_),
	  samples(samples)
{
    //create Z with inducing points placed decided by k-means on X

    srand(time(NULL));
    Z = Matrix(num_inducing,X.Columns());
    for(int i=0;i<num_inducing;i++) {
        Z.setRow(i,X.getRow(rand()%X.Rows()));
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
	Matrix X_samples;
	Vector Y_samples;
	getSamples(X_samples,Y_samples);

	Matrix Kmm_samples = Kernel(Z,Z);
	Matrix Knm_samples = Kernel(X_samples,Z);
	Matrix Kmn_samples = Kernel(Z,X_samples);
	Matrix Knn_samples = Kernel(X_samples,X_samples);

	Matrix invKmm_samples = Kmm_samples.Inverse();
	K_tilde = Knn_samples - Knm_samples * invKmm_samples * Kmn_samples;

	global_update(X_samples,Y_samples);

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

	global_update(X,Y);

	u = MultivariateNormal(m,S).generate();
}
void SVGP::Prediction(const Vector& x, real& mean, real& var) {
	Matrix tmp(1,d);
	tmp.setRow(0,x);
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
void SVGP::local_update() {
	/*
    //here Z should be updated somehow..
    //first, transform Z matrix into double array

    double *data = new double[num_inducing*d];
    for(int i=0;i<num_inducing;i++) {
        for(int j=0;j<d;j++) {
            data[i*m+d] = Z(i,j);
        }
    }
    Matrix sampleMatrix;
    Vector observationVector;
    getSamples(sampleMatrix,observationVector);

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    double size;

    gsl_vector_view x = gsl_vector_view_array(data,num_inducing*d);

    ss = gsl_vector_alloc(num_inducing*d);
    gsl_vector_set_all(ss,1.0);
    
    SVGP* ptr2 = this;
    auto ptr = [=](gsl_vector * x)->real{return ptr2->LogLikelihood(x);};
    gsl_function_pp<decltype(ptr)> Fp(ptr);

    void *par;
    minex_func.n = num_inducing*d;
    minex_func.f = &LogLikelihood;
    minex_func.f = static_cast<gsl_function*>(&Fp);
    minex_func.params = par;

    s = gsl_multimin_fminimizer_alloc (T, num_inducing*d);
    gsl_multimin_fiminimizer_set(s, &minex_func, x, ss);

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if(status)
            break;
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-2);

        if (status = GSL_SUCCESS) {
            printf("converged");
        }
    } while(status == GSL_CONTINUE && iter < 10);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //TODO save best simplex to Z
    delete[] data;
	*/
}
/*
void SVGP::LogLikelihood(double* data, Matrix& sampleMatrix, Vector& observationVector) {
=======
real SVGP::LogLikelihood(const gsl_vector *v) {
>>>>>>> 5e630cfe343b7114924c51f3f9d7b8a71a952741

    Matrix Z_ = Matrix(d,num_inducing);
    for(int i=0;i<num_inducing;i++) {
        for(int j=0;j<d;j++) {
            Z_(i,j) = gsl_vector_get(v,i*d+j);
        }
    }
    Matrix Kmn_samples = Kernel(Z_,currentSample);
    Matrix Knm_samples = Kernel(currentSample,Z_);
    Matrix Kmm_ = Kernel(Z_,Z_);
    Matrix invKmm_ = Kmm_.Inverse();
    Matrix Knn_ = Kernel(currentSample,currentSample);
	Matrix K_tilde_ = Knn_ - Knm_samples * invKmm_ * Kmn_samples;


    Matrix q_prec = Beta * invKmm_ * Kmn_samples * Knm_samples * invKmm_ + invKmm_;
    Matrix q_var = q_prec.Inverse();

	Matrix q_mean_tmp = Beta * q_prec.Inverse() * invKmm_ * Kmn_samples;
	Vector q_mean = q_mean_tmp * currentObservation;

	p_var = Kmm_;
	p_mean = Vector::Null(q_mean.Size());

	Matrix inv_p_var = p_var.Inverse();
	real lhs = (inv_p_var * q_prec).tr();
	Vector A = inv_p_var * (p_mean - q_mean);
	real B = Product((p_mean - q_mean),A);
	KL = lhs + B - N + log(p_var.det()/q_prec.det());
	KL *= 0.5;

	L = 0;
    
	for(int i=0;i<samples;i++) {
		Vector k_i = Kmn_samples.getColumn(i);
		Matrix outer = OuterProduct(k_i,k_i);

		Matrix lambda_i = Beta * invKmm_ * outer * invKmm_;
		real k_ii = K_tilde_(i,i);
		Vector rhs = invKmm_ * q_mean;
		real pdfmean = Product(k_i,rhs);
		real pdfvar = 1/Beta;
		real pdfcond = currentSample(i);

		real pdf = (real) 1/(pdfvar * sqrt(2 * M_PI));
		real expon = exp(-(pdfcond-pdfmean)*(pdfcond-pdfmean)/(2*pdfvar*pdfvar));
		pdf *= expon;
		pdf = (real) log(pdf);
		
		L += pdf - 1/2 * Beta * k_ii - 1/2 * (S * lambda_i).tr();
	}
	L -= KL;
	//return L;
    return 1/L; //to make it a minimzation problem
}
*/
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

void SVGP::getSamples(Matrix& sampleMatrix,Vector& observationVector) {
    sampleMatrix = Matrix(samples,d);
	observationVector = Vector(samples);

    srand(time(NULL));
    for(int i=0;i<samples;i++) {
        int idx = rand()%X.Rows();
        sampleMatrix.setRow(i,X.getRow(idx));
        observationVector(i) = Y(idx); 
    }
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
