#pragma once
#include <vector>
class Node;
class Distribution{
	private:
		double mean;
		double var;
	public:
		Distribution(const double m=0.0,const double v=0.0);
		double getMean() const;
		double getVar() const;
		Distribution operator*(const Distribution &dis);
		Distribution operator/(const Distribution &dis); 
		Distribution operator+(const Distribution &dis);
		Distribution operator-(const Distribution &dis);
		bool  operator==(const Distribution &dis);
		double getSample();
		static Distribution getMin(Distribution &one,Distribution &two,double p=0.0);
		static Distribution getMax(Distribution &one,Distribution &Two,double p=0.0);
		static Distribution getMaxOfCorrelatedSet(std::vector<Node*> &variableSet);
		static Distribution getMaxOfIndependentSet(std::vector<Distribution> &variableSet);
		static double getP(Node &one,Node &two);
		static double getP(int i,int j,double deviation[],double k[],std::vector<Node*> &variables);
		static Distribution getMinOfIndependentSet(std::vector<Distribution> &variableSet);
		static Distribution getMinOfCorrelatedSet(std::vector<Node*>&variableSet);
		double phi(double x);
		double pdf(double x);
		static double normalPhi(double x);
		static double normalPdf(double x);
};
