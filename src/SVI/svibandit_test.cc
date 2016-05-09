#include "svgp.h"
#include "GaussianProcess.h"
#include "SparseGaussianProcess.h"

int main() {
	double x_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };
	double z_data[] = { 0.25, 0.19, 0.40, 0.70,
					  0.41, 0.74, 0.96, 0.18};
	double y_data[] = { 1.0, 2.0, 3.0, 4.0 };

	double x0[] = {0.18,0.6,0.57,0.96};
	double x1[] = {0.41,0.24,0.99,0.58};
	double x2[] = {0.14,0.30,0.97,0.66};
	double x3[] = {0.51,0.13,0.19,0.85};
	double z0[] = {0.25, 0.19, 0.40, 0.70};
	double z1[] = {0.41, 0.74, 0.96, 0.18};
	double y0[] = {1.0,2.0,3.0,4.0};
	std::vector<real> v0(x0,x0+sizeof(x0)/sizeof(real));
	std::vector<real> v1(x1,x1+sizeof(x1)/sizeof(real));
	std::vector<real> v2(x2,x2+sizeof(x2)/sizeof(real));
	std::vector<real> v3(x3,x3+sizeof(x3)/sizeof(real));
	std::vector<real> v4(z0,z0+sizeof(z0)/sizeof(real));
	std::vector<real> v5(z1,z1+sizeof(z1)/sizeof(real));
	Matrix X(4,4);
	X.setRow(0,Vector(v0));
	X.setRow(1,Vector(v1));
	X.setRow(2,Vector(v2));
	X.setRow(3,Vector(v3));
	Matrix Z(2,4);
	Z.setRow(0,Vector(v4));
	Z.setRow(1,Vector(v5));
	Vector Y(4);
	Y(0) = y0[0];
	Y(1) = y0[1];
	Y(2) = y0[2];
	Y(3) = y0[3];

	real mean;
	real var;

	SVGP svgp(X,Y,X.Rows(),1.0,1.0,Vector::Unity(X.Columns()),100);
	//svgp.FullUpdateGaussianProcess();
	svgp.UpdateGaussianProcess();
	Vector gp_scale_length = Vector(X.Columns()) + 1.0;
	SparseGaussianProcess sgp(1.0,gp_scale_length,1.0,1.0,gp_scale_length);
	sgp.Observe(X,Y);
	svgp.Prediction(X.getRow(0),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	sgp.Prediction(X.getRow(0),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	svgp.Prediction(X.getRow(1),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	sgp.Prediction(X.getRow(1),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	svgp.Prediction(X.getRow(2),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	sgp.Prediction(X.getRow(2),mean,var);
	std::cout<<"mean: "<<mean<<", var "<<var<<"\n";
	//GaussianProcess gp(X,Y,1.0,1.0,1.0);
	std::cout<<svgp.LogLikelihood()<<"\n";
	//std::cout<<sgp.LogLikelihood()<<"\n";
	//std::cout<<gp.LogLikelihood()<<"\n";

	return 0;
}
