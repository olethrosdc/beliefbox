#ifndef SVGP_H
#define SVGP_H

#define M_PI 3.14159265358979323846

#include "math.h"
#include "Vector.h"
#include "Matrix.h"

//TODO 
//Change so that svgp uses gsl internally

class SVGP
	{
	protected:
		gsl_matrix * X; ///< Samples (N x d)
		gsl_vector * Y; ///< Output (N x 1)
		int N; ///< Number of samples
		int d;

		// Kernel parameters
		double noise_var;
		double sig_var;
		gsl_vector * scale_length;
		
		/// SVI parameters
		int num_inducing; ///< inducing inputs (m)
		gsl_matrix * Z; ///< Hidden variables (m x d)
		gsl_vector * u; ///< Global variables (m x 1) (targets)

		/// SVI fields
		gsl_vector * q_mean;
		gsl_matrix * q_var;
		gsl_matrix * q_prec;

		// variational distribution parametrization
		gsl_matrix * S;
		gsl_vector * m;

		double l; //step length
		int subsamples; //not currently used.. but should be

		gsl_vector * p_mean;
		gsl_matrix * p_var;

		gsl_matrix * Kmm;
		gsl_matrix * Knm;
		gsl_matrix * Kmn;
		gsl_matrix * Knn;
		gsl_matrix * invKmm;
		gsl_matrix * K_tilde;

		double Beta; //noise precision
		double KL; //KL divergence
		double L; //likelihood

        //preferably these could be done in minibatches but only with single examples for now
        //the example is repeated <subsamples> times as in Hoffman et al [2013]
		virtual void init();
        virtual void local_update(const gsl_matrix *& X_samples, const gsl_vector *& Y_samples); //updating Z (local/latent variables) unsure how to do this 
        virtual void global_update(const gsl_matrix *& X_samples, const gsl_vector *& Y_samples); //updating m, S (which parametrizes q(u) and in turn gives global param u)
		virtual double min_Likelihood(double* data); //to be used as an objective function to minimize likelihood given some Z*
		virtual void Kernel(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * cov); 

	public:
		SVGP(const Matrix& X, const Vector& Y, const Matrix& Z, double noise_var, double sig_var, const Vector& scale_length);
		SVGP(gsl_matrix * X, gsl_vector * Y, gsl_matrix * Z, double noise_var, double sig_var, gsl_vector * scale_length);
		//SVGP(gsl_matrix *& X, gsl_vector *& Y, gsl_vector *& Z);
		virtual void Prediction(const gsl_vector * x, double mean, double var);
		virtual void UpdateGaussianProcess(); //update 
		virtual void FullUpdateGaussianProcess();
		virtual double LogLikelihood();
		virtual void AddObservation(const gsl_vector * x, const double& y);
		virtual void AddObservation(const std::vector<gsl_vector *>& x, const std::vector<double>& y);
		virtual void Clear();

};

#endif
