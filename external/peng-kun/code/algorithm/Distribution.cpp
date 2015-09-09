#include "Distribution.h"
#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include "Node.h"
#include <iomanip>


/**************************************************************************
 * Low level class, deal with specific Gaussian distribution operations
 *It provides some primitive functions to do distribution sampling and distribution mulplication, addition , substraction and division.
 *Also functions that provide MAX/MIN over a set of correlated or independentvariables are included. 
 **************************************************************************
*/


Distribution::Distribution(const double m,const double v){
	mean = m;
	var = v;
}

/*Mulplication of two Gaussian distributions*/
Distribution Distribution::operator*(const Distribution &dis){
	double var1 = dis.getVar();
	double var2 = this->var;
	double mean1 = dis.getMean();
	double mean2 = this->mean;
	/* When one of the distribution is (mean,0)*/
	if(var1==0)
		return  Distribution(mean2,var2);
	else if(var2==0)
		return Distribution(mean1,var1);

	double newVar = 1/(1/var1 + 1/var2);
	double newMean = (mean1/var1 + mean2/var2)*newVar;
	
	return  Distribution(newMean,newVar);


}

Distribution Distribution::operator+(const Distribution &dis){
	double var1 = this->var;
	double mean1 = this->mean;
	
	double var2 = dis.getVar();
	double mean2 = dis.getMean();
	
	return Distribution(mean1+mean2,var1+var2);
}

Distribution Distribution::operator-(const Distribution &dis){
	double var1 = this->var;
	double mean1 = this->mean;
	
	double var2 = dis.getVar();
	double mean2 = dis.getMean();
	
	return Distribution(mean1-mean2,var1-var2);
}

/* Division of two Gaussian distributions. If they are the same, then result   is Distribution(0,0)
*/ 

Distribution Distribution::operator/(const Distribution &dis){
	double var1 = this->var;
	double mean1 = this->mean;

	double var2 = dis.getVar();
	double mean2  =  dis.getMean();	

	/* Consider situations when the variance of two distributions are the same and when variance of the divisor is zero*/
	if(var2==0)
		return Distribution(mean1,var1);
	if(var1 == var2)
		return Distribution(0,0);

	double newVar = 1/(1/var1 - 1/var2);
	double newMean = (mean1/var1 - mean2/var2)*newVar;

	return Distribution(newMean,newVar);

}

/* Two Gaussian distributions are equal if and only if their means and variants are equivalent respectively*/
bool Distribution::operator==(const Distribution& dis){
	double var1 = this->var;
	double mean1 = this->mean;
	
	double var2 = dis.getVar();
	double mean2 = dis.getMean();
	
	if(var1==var2 && mean1==mean2)
		return 1;
	else
		return 0;
}

/*====================================================================
 The following four functions are used for generate value of PDF and CDF.
 Methods normalPdf and normalPhi are static methods, they generate value of  PDF and CDF of standard Gaussian.
 ====================================================================
*/
double Distribution::pdf(double x){
	double a = (x-this->mean)*(x-this->mean);
	double b =  -a /(2*this->var);
	double c =  exp(b);
	double coefficient = pow(M_PI*2*this->var,-0.5);

	return coefficient * c;
}

double Distribution::normalPdf(double x){
	double a = x*x;
	double b = -a/2;
	double c = exp(b);
	double coefficient = pow(M_PI*2,-0.5);
	
	return coefficient*c;
}

double Distribution::normalPhi(double x){
	double a1 =  0.254829592;
	double a2 = -0.284496736;
    	double a3 =  1.421413741;
    	double a4 = -1.453152027;
    	double a5 =  1.061405429;
    	double p  =  0.3275911;
 
    	// Save the sign of x
    	int sign = 1;
    	if(x < 0)
        	sign = -1;
    	x = fabs(x)/sqrt(2.0);
 
 	// A&S formula 7.1.26
    	double t = 1.0/(1.0 + p*x);
    	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    	return 0.5*(1.0 + sign*y);


}

double Distribution::phi(double x){
	x = (x - this->mean)/this->var;
	double a1 =  0.254829592;
	double a2 = -0.284496736;
    	double a3 =  1.421413741;
    	double a4 = -1.453152027;
    	double a5 =  1.061405429;
    	double p  =  0.3275911;
 
    	// Save the sign of x
    	int sign = 1;
    	if(x < 0)
        	sign = -1;
    	x = fabs(x)/sqrt(2.0);
 
 	// A&S formula 7.1.26
    	double t = 1.0/(1.0 + p*x);
    	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    	return 0.5*(1.0 + sign*y);

}

/* Min of two Gaussian distributions; P is the correlated coefficient*/
Distribution Distribution::getMin(Distribution &dis1,Distribution &dis2,double p){
	double var1 = dis1.getVar();
	double mean1 = dis1.getMean();
	double deviation1 = sqrt(var1);

	double var2 = dis2.getVar();
	double mean2 = dis2.getMean();
	double deviation2 = sqrt(var2);
	
	/*procedure to calculate the min*/
	double a = sqrt(var1+var2-2*p*deviation1*deviation2);
	double k;

	/*double value has precision problem, and if p==1 and one==two then Min(one,two) = one*/
	if((dis1 == dis2) && (1-p)<0.0000001){
		return Distribution(mean1,var1);
	}
	else{
		k = (mean1-mean2)/a;	
	}

	double phiK = Distribution::normalPhi(k);
	double phiNegativeK = Distribution::normalPhi(-k);
	double pdfNegativeK = Distribution::normalPdf(-k);
	
	double firstMoment = mean1*phiNegativeK + mean2*phiK - a*pdfNegativeK;
	double secondMoment = (var1+pow(mean1,2))*phiNegativeK + (var2+pow(mean2,2))*phiK - (mean1+mean2)*a*pdfNegativeK;

	return Distribution(firstMoment,secondMoment-pow(firstMoment,2));
}

/* Max of two Gaussian distributions; P is the correlated coefficient*/
Distribution Distribution::getMax(Distribution &dis1,Distribution &dis2,double p){
	double var1 = dis1.getVar();
	double mean1 = dis1.getMean();
	double deviation1 = sqrt(var1);
	
	double var2 = dis2.getVar();
	double mean2 = dis2.getMean();
	double deviation2 = sqrt(var2);
	
	double a = sqrt(var1+var2-2*p*deviation1*deviation2);
	double k;

	/*double value has precision problem, and if p==1 and one==two then Min(one,two) = one*/
	if((dis1 == dis2)&&(1-p)<0.000001){
		return Distribution(mean1,var1);
	}
	else{
		k = (mean1-mean2)/a;	
	}
	double phiK = Distribution::normalPhi(k);
	double phiNegativeK = Distribution::normalPhi(-k);
	double pdfK = Distribution::normalPdf(k);
	
        double firstMoment = mean1*phiK + mean2*phiNegativeK + a*pdfK;
	double secondMoment = (var1+pow(mean1,2))*phiK + (var2+pow(mean2,2))*phiNegativeK + (mean1+mean2)*a*pdfK;
	
	return Distribution(firstMoment,secondMoment-pow(firstMoment,2));
}
	
/* Max of a set of Independent Gaussian variables*/
Distribution Distribution::getMaxOfIndependentSet(std::vector<Distribution> &variables){
	Distribution intermedia = variables[0];
	for(std::vector<Distribution>::iterator it=variables.begin()+1;it<variables.end();it++){
		intermedia = Distribution::getMax(intermedia,(*it));
	}

	return intermedia;

}

/* Min of a set of Independent Gaussian variables*/
Distribution Distribution::getMinOfIndependentSet(std::vector<Distribution> &variables){
	Distribution intermedia = variables[0];
	for(std::vector<Distribution>::iterator it=variables.begin()+1;it<variables.end();it++){
		intermedia = Distribution::getMin(intermedia,(*it));
	}
        return intermedia;
}

/*Max of two correlated Gaussian distribution*/
Distribution Distribution::getMaxOfCorrelatedSet(std::vector<Node*> &variables){
	
	Distribution parentDis =  variables[0]->getParent()->getGDis();
	Distribution intermedia;
	if(variables[0]->getVDis()==Distribution(0,0))
		intermedia = variables[0]->getDelta()+parentDis; 	
	else
		intermedia = variables[0]->getVDis();
	int i=2;
	int j=i-1;
	int size = variables.size();
	double deviations[size-1];
	double k[size-1];
	for(std::vector<Node*>::iterator it=variables.begin()+1;it<variables.end();it++,i++,j++){
		double p = Distribution::getP(i,j,deviations,k,variables);
		//precision problem
		if(p>1)
			p=1.0;
			
		Distribution aR ;
		if((**it).getVDis()==Distribution(0,0))
			aR= (**it).getDelta()+parentDis;
		else
			aR = (**it).getVDis();
		k[j-1] = (intermedia.getMean()-aR.getMean())/sqrt((intermedia.getVar()+aR.getVar()));
		intermedia = getMax(intermedia,aR,p);
		deviations[j-1] = sqrt(intermedia.getVar());
	}
	return intermedia;

}

/*Min of two correlated Gaussian distribution*/
Distribution Distribution::getMinOfCorrelatedSet(std::vector<Node*> &variables){
	Distribution parentDis = variables[0]->getParent()->getGDis();
	Distribution intermedia ;
	if(variables[0]->getVDis()==Distribution(0,0))
		intermedia = variables[0]->getDelta()+parentDis; 	
	else
		intermedia = variables[0]->getVDis();
	int i = 2;
	int j = i-1;
	int size = variables.size();
	double deviations[size-1];
	double k[size-1];
	for(std::vector<Node*>::iterator it=variables.begin()+1;it<variables.end();it++,i++,j++){
		double p = Distribution::getP(i,j,deviations,k,variables);
		//double precision
		if(p>1)
			p=1.0;

		Distribution aR ;
		if((**it).getVDis()==Distribution(0,0))
			aR= (**it).getDelta()+parentDis;
		else
			aR = (**it).getVDis();
		k[j-1] = (intermedia.getMean()-aR.getMean())/sqrt((intermedia.getVar()+aR.getVar()));
		intermedia = getMin(intermedia,aR,p);
		deviations[j-1] = sqrt(intermedia.getVar());
	}
	return intermedia;
}



/* p = var(g)/(var(g)+var(delta1))*(var(g)+var(delta2)))*/
double Distribution::getP(Node &one,Node &two){
	double varParent= one.getParent()->getGDis().getVar();
	
	double devOne= sqrt(one.getDelta().getVar()+varParent);

	double devTwo = sqrt(two.getDelta().getVar()+varParent);


	double coefficient = varParent/(devOne*devTwo);
	return coefficient;
}

double Distribution::getP(int i,int j,double deviation[],double k[],std::vector<Node*> &variables){
	if(j==1)
		return getP(*(variables[i-1]),*(variables[j-1]));
	else{
		double parentV = variables[0]->getParent()->getGDis().getVar();
		double a = 1/deviation[j-2];
		double dev = variables[j-1]->getDelta().getVar()+parentV;
		double b = Distribution(0,1).phi(-k[j-2])*sqrt(dev)*Distribution::getP(*(variables[i-1]),*(variables[j-1]));
		double dev2=0.0; 
		if(j-3<0)
			dev2 = sqrt(variables[j-2]->getDelta().getVar()+parentV);
		else
			dev2 = deviation[j-3];
		
		double c = dev2*Distribution(0,1).phi(k[j-2])*getP(i,j-1,deviation,k,variables);
		return (b+c)*a;
	}
	



}

double Distribution::getSample(){
	std::random_device gen;
    	std::normal_distribution<double> normal(this->mean, sqrt(this->var));
    	return normal(gen);
}
double Distribution::getMean() const{
	return mean;
}

double Distribution::getVar() const{
	return var;
} 

