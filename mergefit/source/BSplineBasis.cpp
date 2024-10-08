#include "BSplineBasis.h"
//#include "KnotInsertion.h"
#include <iostream>
#include <cmath>

using namespace std;

typedef unsigned int uint;

BSplineBasis::BSplineBasis(int deg,const vector<double>& v)
{
	epsilon=1.e-14;
	if(Check(deg,v,epsilon))
	{
		p=deg;
		kv=v;
		nbf=v.size()-deg-1;
	}
	else
	{
		p=3;
		nbf=1;
		double tmp[5]={0.,1.,2.,3.,4.};
		kv.assign(tmp,tmp+5);
	}
	order=p+1;
	SetKnotSpan();
}

void BSplineBasis::Set(int deg,const vector<double>& v)
{
	if(Check(deg,v,epsilon))
	{
		p=deg;
		kv=v;
		nbf=v.size()-deg-1;
		order=p+1;
		SetKnotSpan();
	}
	else
	{
		cerr<<"Fail to set B-spline basis due to the invalid knot vector!\n";
		abort();
	}

}

bool BSplineBasis::Check(int deg,const vector<double>& v, double eps) const
{
	int ord(deg+1);
	if(deg<0)
	{
		cerr<<"Negative degree!\n";
		return false;
	}
	if(v.size()<=ord)
	{
		cerr<<"Length of knot vector should be at least p+2!\n";
		return false;
	}
	else
	{
		int flag(0);
		for(uint i=0;i<v.size()-1;i++)
		{
			if(v[i]-v[i+1]>eps)
			{
				flag=1;
				break;
			}
		}
		if(flag==1)
		{
			cerr<<"Knot vector is not in increasing order!\n";
			return false;
		}
		else
		{
			int rep(0);
			//check first knot
			for(int i=0; i <= ord; i++)
			{
				if(fabs(v[i]-v[0])<eps) rep++;
			}
			if(rep > ord)
			{
				cerr<<"First knot repeated more than "<<deg<<"+1 times\n";
				return false;
			}
			//check last knot
			rep=0;
			for(int i=0; i <= ord; i++)
			{
				if(fabs(v[v.size()-1-i]-v.back())<eps) rep++;
			}
			if(rep > ord)
			{
				cerr<<"Last knot repeated more than "<<deg<<"+1 times\n";
				return false;
			}
			//check interior knots
			for(int i=ord; i<v.size()-ord; i++)
			{
				rep=0;
				for(int j=i; j<ord; j++)
				{
					if(fabs(v[j]-v[i])<eps) rep++;
				}
				if(rep>deg)
				{
					cerr<<"Interior knot at position " << i <<  " repeated more than "<<p<<" times\n";
					return false;
				}
			}
		}
	}
	return true;
}

void BSplineBasis::SetKnotSpan()
{
	ks.clear();
	ks.reserve(kv.size()-2*p-1); // open knot vector
	for(int i=0; i<kv.size()-1; i++)
	{
		if(kv[i+1]-kv[i] > epsilon)
		{
			ks.push_back(i);
		}
	}
}

int BSplineBasis::FindSpan(double par) const
{
	if(par < kv.front() || par > kv.back())
	{
		cerr<<"The parametric value "<< par << " is out of range of knot vector!\n";
		return -1;
	}
	if(fabs(par-kv[0]) < epsilon)
	{
		return ks[0];
	}
	if(fabs(par-kv.back()) < epsilon)
	{
		return ks.back();
	}

//	int ilow(0), ihigh(ks.size()), imid((ilow+ihigh)/2);
//	cout<<"par: "<<par<<"\n";
//	cout<<"kv: ";
//	for(auto kt:kv) cout<<kt<<" ";
//	cout<<"\n";
//	while(par < kv[ks[imid]] || par >= kv[ks[imid]+1])
//	{
//		if(par < kv[ks[imid]]) ihigh = imid;
//		else ilow = imid;
//		imid = (ilow+ihigh)/2;
//		cout<<ilow<<" "<<ihigh<<" "<<imid<<"\n";
//	}
//	cout<<"\n"; getchar();

	int pos(-1);
	for(uint i=0; i<ks.size(); i++)
	{
		if(par>=kv[ks[i]] && par<kv[ks[i]+1])
		{
			pos=ks[i];
			break;
		}
	}
	if(pos<0)
	{
		cerr<<"Can't find the span including " << par << "!\n";
		return -1;
	}
	return pos;
}

void BSplineBasis::BasisFunction(int bid,double par,int deriv,vector<double>& val) const
{
	if(bid<0 || bid>=nbf)
	{
		cerr<<"Wrong basis functions ID: " << bid << "!\n";
		return;
	}
	if(kv.front()-par > epsilon || par-kv.back() > epsilon)
	{
		cerr<<"Wrong parametric value: "<< par << "!\n";
		return;
	}
	if(deriv>1)
	{
		cerr<<"Wrong derivative order: "<< deriv << "!\n";
		return;
	}
	val.clear();
	val.resize(deriv+1,0.);

	vector<double> N0(p+1,0.);
	double left,right,d1;
	int i,j,q;
	if(fabs(par-kv[bid])<epsilon)
	{
		par=kv[bid];
		N0[0]=1.;
	}
	else if(fabs(par-kv[bid+p+1])<epsilon)
	{
		par=kv[bid+p+1];
		N0[p]=1.;
	}
	else
	{
		for(i=bid; i<bid+p+1; i++)
		{
			if(par>=kv[i] && par<kv[i+1])
			{
				N0[i-bid]=1.;
				break;
			}
		}
	}

	//recursive
	for(q=1; q<=p; q++)
	{
		j = p+1-q;
		for(i=0; i<j; i++)
		{
			if(fabs(kv[bid+i+q]-kv[bid+i]) < epsilon) left=0.;
			else left = (par-kv[bid+i])/(kv[bid+i+q]-kv[bid+i]);
			if(fabs(kv[bid+i+q+1]-kv[bid+i+1]) < epsilon) right=0.;
			else right = (kv[bid+i+q+1]-par)/(kv[bid+i+q+1]-kv[bid+i+1]);
			if(q==p)//first derivative
			{
				double left1,right1;
				if(fabs(kv[bid+i+q]-kv[bid+i]) < epsilon) left1=0.;
				else left1 = q/(kv[bid+i+q]-kv[bid+i]);
				if(fabs(kv[bid+i+q+1]-kv[bid+i+1]) < epsilon) right1=0.;
				else right1 = q/(kv[bid+i+q+1]-kv[bid+i+1]);
				d1=left1*N0[0]-right1*N0[1];
			}
			N0[i] = left*N0[i] + right*N0[i+1];
		}
	}
	val[0]=N0[0];
	if(deriv==1) val[1]=d1;
}

void BSplineBasis::BasisFunctionVal(int i, double par, vector<double>& val) const
{
	val.clear();
	val.resize(order,0.);
	val[0]=1.;
	vector<double> left(order,0.), right(order,0.);
	double tmp, saved;
	int j, r;
	for(j=1; j<=p; j++)
	{
		left[j] = par-kv[i+1-j];
		right[j] = kv[i+j] - par;
		saved = 0.;
		for(r=0; r<j; r++)
		{
			tmp = val[r] / (right[r+1] + left[j-r]);
			val[r] = saved + right[r+1] * tmp;
			saved = left[j-r] * tmp;
		}
		val[j] = saved;
	}
}

void BSplineBasis::BasisFunctionDer(int i, double par, vector<vector<double>>& der, int n0) const
{
	if(n0 <= 0)
	{
		cerr<<"Derivative order "<< n0 << " is not larger than 0!\n";
		return;
	}
//	else if(n0 > p)
//	{
//		cerr<<"Derivatives order " << n0 << " is higher than spline degree "<< p << "!\n";
//	}
	int n = n0 < p ? n0 : p;
	der.clear();
	der.resize(n0+1,vector<double>(order,0.));

	vector<double> left(order,0.), right(order,0.);
	double tmp, saved, d;
	vector<vector<double>> ndu(order, vector<double>(order,0.));
	vector<vector<double>> a(2, vector<double>(order,0.));
	ndu[0][0] = 1.;
	int j, k, r, rk, pk, j1, j2, s1, s2;
	for(j=1; j<=p; j++)
	{
		left[j] = par-kv[i+1-j];
		right[j] = kv[i+j] - par;
		saved = 0.;
		for(r=0; r<j; r++)
		{
			ndu[j][r] = right[r+1] + left[j-r];
			tmp = ndu[r][j-1]/ndu[j][r];
			ndu[r][j] = saved + right[r+1] * tmp;
			saved = left[j-r] * tmp;
		}
		ndu[j][j] = saved;
	}
	for(j=0; j<=p; j++)
	{
		der[0][j] = ndu[j][p];
	}
	for(r=0; r<=p; r++)
	{
		s1 = 0;
		s2 = 1;
		a[0][0] = 1.;
		for(k=1; k<=n; k++)
		{
			d = 0.;
			rk = r - k;
			pk = p - k;
			if(r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk+1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if(rk >= -1) j1 = 1;
			else j1 = -rk;
			if(r-1 <= pk) j2= k-1;
			else j2 = p-r;
			for(j=j1; j<=j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
				d += a[s2][j] * ndu[rk+j][pk];
			}
			if(r <= pk)
			{
				a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			der[k][r] = d;
			j = s1;
			s1 = s2;
			s2 = j;
		}
	}

	r = p;
	for(k=1; k<=n; k++)
	{
		for(j=0; j<=p; j++)
		{
			der[k][j] *= r;
		}
		if(p > k) r *= (p - k);//this is different from book
	}
}

void BSplineBasis::get_greville_points(vector<double> &grev) const
{
    grev.clear();
    grev.resize(nbf);
    for (int i=0; i<nbf; i++)
    {
        grev[i] = 0.;
        for (int j=0; j<p; j++)
        {
            grev[i] += kv[i+1+j];
        }
        grev[i] /= p;
    }
}

















