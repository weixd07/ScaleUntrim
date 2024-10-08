#ifndef BSPLINE_BASIS_H
#define BSPLINE_BASIS_H

#include <array>
#include <vector>
//#include <utility>
//#include <string>

using namespace std;

class BSplineBasis
{
public:
	BSplineBasis(int deg,const vector<double>& v);
	void SetEpsilon(double tol){ epsilon=tol; }
	void Set(int deg,const vector<double>& v);
	bool Check(int deg,const vector<double>& v, double eps=1.e-14) const;//valid degree and knot vector
	void SetKnotSpan();
	int FindSpan(double par) const;
//	int Evaluate(double par,int deriv,vector<double>& val);
	void BasisFunction(int pos,double par,int deriv,vector<double>& val) const;
	void BasisFunctionVal(int i, double par, vector<double>& val) const;
	void BasisFunctionDer(int i, double par, vector<vector<double>>& der, int n0=1) const;
	//void test();
    void get_greville_points(vector<double>& grev) const;

    BSplineBasis() : p(0), order(1), nbf(0), epsilon(1.e-14) {}
	
	int p;//polynomial order
	int nbf;//# of basis functions
	int order;
	vector<double> kv;//knot vector
	vector<int> ks; // knot span
	double epsilon;
};

#endif