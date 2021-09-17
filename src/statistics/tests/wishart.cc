#ifdef MAKE_MAIN
#include "ContextTreeKDTree.h"
#include "ContextTreeRealLine.h"
#include "Random.h"
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "MultivariateNormal.h"
#include "MultivariateNormalUnknownMeanPrecision.h"
#include "BetaDistribution.h"
#include "iWishart.h"


int main (int argc, char** argv)
{
	Vector mean(2); mean[0] = 1; mean[1] = 1;
	Matrix sigma(2,2);
	sigma(0,0) = 100;
	sigma(0,1) = 1;
	sigma(1,0) = 1;
	sigma(1,1) = 1;
//	sigma.print(stdout);sigma.Inverse().print(stdout); Matrix mat = sigma.Inverse() * sigma; mat.print(stdout); 
//    real det; sigma.LUDecomposition(det); printf("det:%f\n",det); 

	int n = 10000;
	Matrix X(n,2);
	int count = 0;
    MultivariateNormal m_normal(mean, sigma.Inverse());
    real log_sum = 0;
	for(int i = 0; i<n; ++i) {
        X.setRow(i,m_normal.generate());
        log_sum += m_normal.log_pdf(X.getRow(i));
//        X.getRow(i).print(stdout);
//        printf("pdf:%f\n",F(i,0));
	}
	printf("log-sum:%f\n",log_sum);
	printf("X\n");
//	X.print(stdout);

    Matrix S0(4,4);
    S0(0,0) = 10; S0(0,1) = 1; S0(0,2) = 1; S0(0,3) = 1;
    S0(1,0) = 1; S0(1,1) = 2; S0(1,2) = 1; S0(1,3) = 1;
    S0(2,0) = 1; S0(2,1) = 1; S0(2,2) = 2; S0(2,3) = 1;
    S0(3,0) = 1; S0(3,1) = 1; S0(3,2) = 1; S0(3,3) = 2;
 //   real det; S0.LUDecomposition(det); printf("det:%f\n",det);
    int N = 100;
    iWishart iwish = iWishart(N, S0, true);
    Wishart wish = Wishart(N, S0, true);
    log_sum = 0; Matrix expected_X = Matrix::Null(4,4);
	for(int i = 0; i<n; ++i) {
        //Matrix V = iwish.generate();
		Matrix V = wish.generate();

        expected_X += V; // For chking iWihart and Wishart r.v.s generated from their own distributions
		//expected_X += V.Inverse(); //For chking expected value of iWishart r.v. generated using Wishart
//        V = S0+Matrix::Unity(4,4)*10;
//        V.print(stdout);
        //real log_pdf = iwish.log_pdf(V.Inverse());

        //Calculating and adding Jacobian in log-scale
        //real det = V.det();
        //real log_pdf_wish = wish.log_pdf(V) + 5*log(det); //Jacobian is det|V|**(p+1). Refer Anderson book in ADP thread

 //       printf("%f %f\n",log_pdf,log_pdf_wish);
  //      if (log_pdf > 0) printf("i:%d val:%f\n",i,log_pdf);
		//        log_sum += log_pdf;
	}
    expected_X = expected_X/(1.0*n*N);
	expected_X.print(stdout);
	//expected_X.Inverse().print(stdout);

	//<< Expected value of Wishart
	//expected_X = ((N-4-1)*expected_X/n).Inverse(); expected_X.print(stdout); //<< Expected value of iWishart
	//printf("log-sum:%f\n",log_sum);

    return 0;
}

#endif
