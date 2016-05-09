#ifndef SVGP_H
#define SVGP_H

#define M_PI 3.14159265358979323846

#include "Vector.h"
#include "Matrix.h"
#include "GaussianProcess.h"
#include <gsl/gsl_math.h>

class SVGP //: public GaussianProcess
	{
	protected:
		Matrix X; ///< Samples (N x M)
		Vector Y; ///< Output (N x 1)
		int N; ///< Number of samples
		int d;

		// Kernel parameters
		real noise_var;
		real sig_var;
		Vector scale_length;
		
		/// SVI parameters
		int num_inducing; ///< inducing inputs (J)
		Matrix Z; ///< Hidden variables (N x J)
		Vector u; ///< Global variables (N x 1) (targets)

		// variational distribution parametrization
		Matrix S;
		Vector m;

		real l; //step length
		int samples; //not currently used.. but should be

		Vector p_mean;
		Matrix p_var;

		Matrix Kmm;
		Matrix Knm;
		Matrix Kmn;
		Matrix Knn;
		Matrix invKmm;
		Matrix K_tilde;

		real Beta;
		real KL;
		real L;

        Matrix currentSample;
        Vector currentObservation;

        //preferably these could be done in minibatches but only with single examples for now
        //the example is repeated <subsamples> times as in Hoffman et al [2013]
        virtual void local_update(); //updating Z (local/latent variables)
        virtual void global_update(const Matrix& X_samples, const Vector& Y_samples); //updating m, S (which parametrizes q(u) and in turn gives global param u)
		//virtual Vector optimize_Z(int max_iters);
        virtual void init();
        //virtual real LogLikelihood(double* data); //computes the likelihood given a Z*, to use for optimizing the position of Z
        virtual void getSamples(Matrix& X_, Vector& Y_);
	public:
		SVGP(Matrix& X, Vector& Y, Matrix& Z, real noise_var, real sig_var, Vector scale_length, int samples);
        //initializes Z through k-means
		SVGP(Matrix& X, Vector& Y, int num_inducing, real noise_var, real sig_var, Vector scale_length, int samples);
		virtual Matrix Kernel(const Matrix& A, const Matrix& B);
		virtual void Prediction(const Vector& x, real& mean, real& var);
		virtual void UpdateGaussianProcess(); //update 
		virtual void FullUpdateGaussianProcess();
		virtual real LogLikelihood();
        //virtual real LogLikelihood(const gsl_vector *v);
		virtual void AddObservation(const Vector& x, const real& y);
		virtual void AddObservation(const std::vector<Vector>& x, const std::vector<real>& y);
		virtual void Clear();

};
#endif
