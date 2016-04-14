#include "svgp.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

SVGP::SVGP(const Matrix& X_, const Vector& Y_, const Matrix& Z_, double noise_var_, double sig_var_, const Vector& scale_length_)
{
	Beta = 1/(noise_var*noise_var);
	N = X_.Rows();
	d = X_.Columns();
	num_inducing = Z_.Rows();

	assert(d==Z_.Columns());
	assert(N==Y_.Size());
	assert(d==scale_length_.Size());

	X = gsl_matrix_alloc (N,d);
	Y = gsl_vector_alloc (N);
	Z = gsl_matrix_alloc (num_inducing,d);
	scale_length = gsl_vector_alloc (d);

	for(int i = 0;i<N;i++) {
		for(int j = 0;j<d;j++) {
			gsl_matrix_set(X,i,j,X_(i,j));
			gsl_matrix_set(Z,i,j,Z_(i,j));
		}
		gsl_vector_set(Y,i,Y_(i));
		gsl_vector_set(scale_length,i,scale_length_(i));
	}
	FullUpdateGaussianProcess();
}
SVGP::SVGP(gsl_matrix * X_, gsl_vector * Y_, gsl_matrix * Z_, double noise_var_, double sig_var_, gsl_vector * scale_length_)
	//TODO make sure data is copied
	: X(X_),
	  Y(Y_),
	  Z(Z_),
	  noise_var(noise_var_),
	  sig_var(sig_var_),
	  scale_length(scale_length_)
{
	Beta = 1/(noise_var*noise_var);
	N = X->size1;
	d = X->size2;
	num_inducing = Z->size1;

	assert(d==Z->size2);
	assert(N==Y->size);
	assert(d==scale_length->size);

	FullUpdateGaussianProcess();
}
void SVGP::init() {
	Kmm = gsl_matrix_alloc(num_inducing,num_inducing);
	Kmn = gsl_matrix_alloc(num_inducing,N);
	Knm = gsl_matrix_alloc(N,num_inducing);
	Knn = gsl_matrix_alloc(N,N);
	invKmm = gsl_matrix_alloc(num_inducing,num_inducing);
	K_tilde = gsl_matrix_alloc(N,N);

	q_var = gsl_matrix_alloc(num_inducing,num_inducing);
	q_mean = gsl_vector_alloc(num_inducing);
	p_var = gsl_matrix_alloc(num_inducing,num_inducing);
	p_mean = gsl_vector_alloc(num_inducing);

	S = gsl_matrix_alloc(num_inducing,num_inducing);
	m = gsl_vector_alloc(num_inducing);

	u = gsl_vector_alloc(num_inducing);
}
void SVGP::Clear() {
	gsl_matrix_free(Kmm);
	gsl_matrix_free(Kmn);
	gsl_matrix_free(Knm);
	gsl_matrix_free(Knn);
	gsl_matrix_free(invKmm);
	gsl_matrix_free(K_tilde);

	gsl_matrix_free(q_var);
	gsl_vector_free(q_mean);
	gsl_matrix_free(p_var);
	gsl_vector_free(p_mean);

	gsl_matrix_free(S);
	gsl_vector_free(m);

	gsl_vector_free(u);

	gsl_matrix_free(X);
	gsl_vector_free(Y);
	gsl_matrix_free(Z);
	gsl_vector_free(scale_length);
}
void SVGP::Kernel(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * cov) {
	assert(A->size2==B->size2);

	int d_ = A->size2;
	int N_ = A->size1;
	int M_ = B->size1;

	cov = gsl_matrix_alloc (N_, M_);

	double extra_noise = 0;
	if(N_ == M_)
		extra_noise = noise_var*noise_var;

	double sig_noise = sig_var*sig_var;
	for(int i=0; i<M_; ++i) {
		gsl_vector * x = gsl_vector_alloc (d_);
		gsl_matrix_get_row(x,B,i);
		for(int j=0; j<N_; j++) {
			gsl_vector * A_j = gsl_vector_alloc (d_);
			gsl_matrix_get_row(A_j,A,j);
			gsl_vector * x_ = gsl_vector_alloc (d_);
			gsl_vector_memcpy(x_,x);
			gsl_vector_sub(x_,A_j);
			gsl_vector_div(x_,scale_length);
			double delta = gsl_blas_dnrm2(x_);
			delta = sig_noise*exp(-0.5*delta);
			if(i == j) 
				delta += extra_noise;
			gsl_matrix_set(cov,j,i,delta);
			gsl_vector_free(x_);
			//cov(i,j) = delta;
		}
		gsl_vector_free(x);
	}
}

void SVGP::UpdateGaussianProcess() {

	Kernel(Z,Z,Kmm);
	Kernel(X,Z,Knm);
	Kernel(Z,X,Kmn);

	/* K_tilde = Knn - Knm * invKmm * Kmn; */

	gsl_matrix_memcpy(invKmm,Kmm);
	gsl_linalg_cholesky_decomp(invKmm);
	gsl_linalg_cholesky_invert(invKmm);

	gsl_matrix * tmpmn = gsl_matrix_alloc(num_inducing,N);
	gsl_matrix * tmpnn = gsl_matrix_alloc(N,N);

	gsl_matrix_view A = gsl_matrix_view_array(invKmm->data, num_inducing, num_inducing);
	gsl_matrix_view B = gsl_matrix_view_array(Kmn->data, num_inducing, N);
	gsl_matrix_view C = gsl_matrix_view_array(tmpmn->data, num_inducing, N);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

	A = gsl_matrix_view_array(Knm->data, N, num_inducing);
	B = gsl_matrix_view_array(tmpmn->data,num_inducing,N);
	C = gsl_matrix_view_array(tmpnn->data,N,N);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

	gsl_matrix_memcpy(K_tilde,Knn);
	gsl_matrix_sub(K_tilde,&C.matrix);

	gsl_matrix_free(tmpmn);
	gsl_matrix_free(tmpnn);

	// generate m multivariate samples and save in u

	gsl_matrix * sigmachol = gsl_matrix_alloc(num_inducing,num_inducing);
	gsl_matrix_memcpy(sigmachol,S);
	gsl_linalg_cholesky_decomp(sigmachol);

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	for(int i = 0;i<num_inducing;i++) {
		gsl_vector_set(u,i,gsl_ran_ugaussian(r));
	}

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, sigmachol, u);
	gsl_vector_add(u,m);
	
	gsl_matrix_free(sigmachol);
	gsl_rng_free(r);
}
void SVGP::FullUpdateGaussianProcess() {

	//expensive stuff

	Kernel(X,X,Knn);

	// local_update
	// global_update

	UpdateGaussianProcess();
}
void SVGP::Prediction(const gsl_vector * x, double& mean, double& var) {
	gsl_matrix * tmp = gsl_matrix_alloc(1,d);
	gsl_matrix_set_row(tmp,0,x);
	gsl_matrix * Kx = gsl_matrix_alloc(1,num_inducing);
	gsl_vector * Kxv = gsl_vector_alloc(num_inducing);
	Kernel(tmp,Z,Kx);

	gsl_matrix_get_row(Kxv,Kx,0);
	gsl_vector * tmpv = gsl_vector_alloc(d);
	gsl_blas_dgemv( CblasNoTrans, 1.0, invKmm, u, 0.0, tmpv);

	double * res;
	gsl_blas_ddot(Kxv,tmpv,res);

	mean = *res;

	gsl_blas_dgemv( CblasNoTrans, 1.0, invKmm, Kxv, 0.0, tmpv);
	gsl_blas_ddot(Kxv,tmpv,res);

	var = *res;
	var = sig_var*sig_var - var;

	gsl_matrix_free(tmp);
	gsl_matrix_free(Kx);
	gsl_vector_free(Kxv);
}
double SVGP::LogLikelihood() {
	//D_KL(q(u)||p(u))

	q_prec = Beta * invKmm * Kmn * Knm * invKmm + invKmm;
	q_var = q_prec.Inverse();

	gsl_matrix q_mean_tmp = Beta * q_prec.Inverse() * invKmm * Kmn;
	q_mean = q_mean_tmp * Y;

	p_var = Kmm;
	p_mean = gsl_vector::Null(q_mean.Size());

	gsl_matrix inv_p_var = p_var.Inverse();
	double lhs = (inv_p_var * q_prec).tr();
	gsl_vector A = inv_p_var * (p_mean - q_mean);
	double B = Product((p_mean - q_mean),A);
	KL = lhs + B - N + log(p_var.det()/q_prec.det());
	KL *= 0.5;

	L = 0;
	for(int i=0;i<N;i++) {
		gsl_vector k_i = Kmn.getColumn(i);
		gsl_matrix outer = OuterProduct(k_i,k_i);

		gsl_matrix lambda_i = Beta * invKmm * outer * invKmm;
		double k_ii = K_tilde(i,i);
		gsl_vector rhs = invKmm * q_mean;
		double pdfmean = Product(k_i,rhs);
		double pdfvar = 1/Beta;
		double pdfcond = Y(i);

		double pdf = (double) 1/(pdfvar * sqrt(2 * M_PI));
		double expon = exp(-(pdfcond-pdfmean)*(pdfcond-pdfmean)/(2*pdfvar*pdfvar));
		pdf *= expon;
		pdf = (double) log(pdf);
		
		L += pdf - 1/2 * Beta * k_ii - 1/2 * (S * lambda_i).tr();
	}
	L -= KL;
	return L;
}

double SVGP::min_Likelihood(double* data) {
	gsl_matrix * Z_ = gsl_matrix_alloc (num_inducing, d);
	gsl_matrix * q_prec_ = gsl_matrix_alloc (num_inducing, num_inducing);
	gsl_matrix * invKmm_ = gsl_matrix_alloc (num_inducing, num_inducing);
	gsl_matrix * Kmn_ = gsl_matrix_alloc (num_inducing, num_inducing);
	Z_->data = data;
}
//OK, here is where I'm a bit confused on how to update the inducing inputs Z
void SVGP::local_update(const gsl_matrix& X_samples, const gsl_vector& Y_samples) {
  /*
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  double par[2] = { 1.0, 2.0 };

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = 2;
  my_func.params = &par;

  Starting point, x = (5,7) 
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);

  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fdfminimizer_alloc (T, 2);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, 1e-3);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

      printf ("%5d %.5f %.5f %10.5f\n", iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1), 
              s->f);

    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return 0;
  */
}
void SVGP::global_update(const gsl_matrix& X_samples, const gsl_vector& Y_samples) {
    gsl_matrix Kmn_samples = Kernel(Z,X_samples);

	gsl_vector theta_1 = S.Inverse() * m;
	gsl_matrix theta_2 = -0.5 * S.Inverse();

	gsl_vector tmpResult = Beta * invKmm * Kmn_samples * Y_samples - theta_1;
	tmpResult *= l;
	tmpResult += theta_1;
	theta_1 = tmpResult;
	theta_2 = theta_2 + l * (-1/2 * q_prec + 1/2 * theta_2);


	S = theta_2.Inverse();
	S *= -2.0;

	//LUDecomp and LUSolve

	gsl_matrix A_ = S.Inverse();
	int size = A_.Rows();
	double* a_data = new double[size*size];
	double* b_data = new double[size];
	
	for(int i = 0;i<size;i++) {
		int offset = i*size;
		for(int j = 0;j<size;j++) {
			a_data[offset+j] = A_(i,j);
		}
		b_data[i] = theta_1(i);
	}
	gsl_matrix_view m_ 
    = gsl_matrix_view_array (a_data, size, size);

	gsl_vector_view b
    = gsl_vector_view_array (b_data, size);

	gsl_vector *x = gsl_vector_alloc (size);
  
	int s;

	gsl_permutation * p = gsl_permutation_alloc (size);

	gsl_linalg_LU_decomp (&m_.gsl_matrix, p, &s);

	gsl_linalg_LU_solve (&m_.gsl_matrix, p, &b.gsl_vector, x);

	gsl_permutation_free (p);

	gsl_vector m = gsl_vector::Null(size);;
	m.x = x->data;

	gsl_vector_free (x);
	delete [] a_data;
	delete [] b_data;

}
void SVGP::AddObservation(const std::gsl_vector<gsl_vector>& x, const std::gsl_vector<double>& y) {
	//assert(x[0].Size()==d);
	for(int i = 0;i<x.size();i++) {
		gsl_vector x_ = x[i];
		double y_ = y[i];
		AddObservation(x_,y_);
	}
	FullUpdateGaussianProcess();
}
void SVGP::AddObservation(const gsl_vector& x, const double& y) {
	X.AddRow(x);
	Y.AddElement(y);
}
