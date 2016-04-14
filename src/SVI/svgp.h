#ifndef SVGP_H
#define SVGP_H

#define M_PI 3.14159265358979323846

#include "Vector.h"
#include "Matrix.h"
#include "GaussianProcess.h"

class SVGP : public GaussianProcess
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

		/// SVI fields
		Vector q_mean;
		Matrix q_var;
		Matrix q_prec;

		// variational distribution parametrization
		Matrix S;
		Vector m;

		real l; //step length
		int subsamples; //not currently used.. but should be

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

        //preferably these could be done in minibatches but only with single examples for now
        //the example is repeated <subsamples> times as in Hoffman et al [2013]
        virtual void local_update(const Matrix& X_samples, const Vector& Y_samples); //updating Z (local/latent variables) unsure how to do this 
        virtual void global_update(const Matrix& X_samples, const Vector& Y_samples); //updating m, S (which parametrizes q(u) and in turn gives global param u)
		virtual Vector optimize_Z(int max_iters);

	public:
		SVGP(Matrix& X, Vector& Y, Matrix& Z, real noise_var, real sig_var, Vector scale_length);
		//SVGP(Matrix& X, Vector& Y, Vector& Z);
		virtual Matrix Kernel(const Matrix& A, const Matrix& B);
		virtual void Prediction(const Vector& x, real& mean, real& var);
		virtual void UpdateGaussianProcess(); //update 
		virtual void FullUpdateGaussianProcess();
		virtual real LogLikelihood();
		virtual void AddObservation(const Vector& x, const real& y);
		virtual void AddObservation(const std::vector<Vector>& x, const std::vector<real>& y);
		virtual void Clear();

};

#endif