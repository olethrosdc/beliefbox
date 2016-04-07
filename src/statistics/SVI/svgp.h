#ifndef SVGP_H
#define SVGP_H

#define M_PI 3.14159265358979323846

#include "Vector.h"
#include "Matrix.h"

class SVGP
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

	public:
		SVGP(Matrix& X, Vector& Y, Matrix& Z, real noise_var, real sig_var, Vector scale_length);
		//SVGP(Matrix& X, Vector& Y, Vector& Z);
		virtual void optimize_Z(); //optimize the positioning of the inducing inputs using likelihood
		virtual void optimize_svi(); //optimize the variational distribution parameters S, m
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