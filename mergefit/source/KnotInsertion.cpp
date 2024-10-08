#include "KnotInsertion.h"
#include <iostream>
#include <algorithm>
#include <cmath>

void InitialT(const vector<double>& k1,const vector<double>& k2,vector<vector<double> >& T)
{
	unsigned int i,j;
	T.clear();
	T.resize(k2.size()-1,vector<double>(k1.size()-1,0.));
	for(i=0;i<k2.size()-1;i++)
		for(j=0;j<k1.size()-1;j++)
			if(k2[i]>=k1[j] && k2[i]<k1[j+1])
				T[i][j]=1.;
}

void KnotInsert(const vector<double>& k1,const vector<double>& k2,const vector<vector<double> >& T1,int q,vector<vector<double> >& T2)
{
	unsigned int i,j;
	T2.clear();
	T2.resize(k2.size()-q-1,vector<double>(k1.size()-q-1,0.));
	//T2=MatrixXd::Zero(k2.size()-q-1,k1.size()-q-1);
	for(i=0;i<k2.size()-q-1;i++)
		for(j=0;j<k1.size()-q-1;j++)
		{
			if(k1[j+q]==k1[j] && k1[j+q+1]!=k1[j+1])
				T2[i][j]=(k1[j+q+1]-k2[i+q])/(k1[j+q+1]-k1[j+1])*T1[i][j+1];
			if(k1[j+q]!=k1[j] && k1[j+q+1]==k1[j+1])
				T2[i][j]=(k2[i+q]-k1[j])/(k1[j+q]-k1[j])*T1[i][j];
			if(k1[j+q]!=k1[j] && k1[j+q+1]!=k1[j+1])
				T2[i][j]=(k2[i+q]-k1[j])/(k1[j+q]-k1[j])*T1[i][j]+(k1[j+q+1]-k2[i+q])/(k1[j+q+1]-k1[j+1])*T1[i][j+1];
		}
}

void TMatrix(const vector<double>& k1,const vector<double>& k2,int p,vector<vector<double> >& T)
{
	vector<vector<double> > T_tmp;
	InitialT(k1,k2,T_tmp);
	for(int n=1;n<=p;n++)
	{
		KnotInsert(k1,k2,T_tmp,n,T);
		T_tmp.clear();
		T_tmp=T;
	}
}

void InitialT(const vector<double>& k1,const vector<double>& k2, SparseMatrix<double>& T)
{
	unsigned int i,j;
	T.resize(k2.size()-1,k1.size()-1);
	vector<Triplet<double>> trip;
	trip.reserve((k2.size()-1)*(0+2));
	for(i=0;i<k2.size()-1;i++)
	{
		for(j=0;j<k1.size()-1;j++)
		{
			if(k2[i]>=k1[j] && k2[i]<k1[j+1])
			{
				trip.push_back(Triplet<double>(i,j,1.));
			}
		}
	}
	T.setFromTriplets(trip.begin(),trip.end());
}

void KnotInsert(const vector<double>& k1,const vector<double>& k2, const SparseMatrix<double>& T1,int q, SparseMatrix<double>& T2)
{
	unsigned int i,j;
	T2.resize(k2.size()-q-1,k1.size()-q-1);
	vector<Triplet<double>> trip;
	trip.reserve((k2.size()-1)*(q+2));
	double tmp;
	for(i=0;i<k2.size()-q-1;i++)
	{
		for(j=0;j<k1.size()-q-1;j++)
		{
			if(k1[j+q]==k1[j] && k1[j+q+1]!=k1[j+1])
			{
				tmp = (k1[j+q+1]-k2[i+q])/(k1[j+q+1]-k1[j+1]) * T1.coeff(i,j+1);
				trip.push_back(Triplet<double>(i,j,tmp));
			}
			if(k1[j+q]!=k1[j] && k1[j+q+1]==k1[j+1])
			{
				tmp = (k2[i+q]-k1[j])/(k1[j+q]-k1[j]) * T1.coeff(i,j);
				trip.push_back(Triplet<double>(i,j,tmp));
			}
			if(k1[j+q]!=k1[j] && k1[j+q+1]!=k1[j+1])
			{
				tmp = (k2[i+q]-k1[j])/(k1[j+q]-k1[j])*T1.coeff(i,j)+(k1[j+q+1]-k2[i+q])/(k1[j+q+1]-k1[j+1]) * T1.coeff(i,j+1);
				trip.push_back(Triplet<double>(i,j,tmp));
			}
		}
	}
	T2.setFromTriplets(trip.begin(),trip.end());
}

void TMatrix(const vector<double>& k1,const vector<double>& k2,int p, SparseMatrix<double>& T)
{
	SparseMatrix<double> T_tmp;
	InitialT(k1,k2,T_tmp);
	for(int n=1; n<=p; n++)
	{
		KnotInsert(k1,k2,T_tmp,n,T);
		T_tmp.resize(T.rows(),T.cols());
		T_tmp = T;
	}
}

void InsertKnots(const vector<double>& kv,const vector<double>& knots,vector<double>& kv1)
{
	kv1.clear();
	if(knots.size()==0)
	{
		kv1=kv;
		return;
	}
	unsigned int i(0),j(0);
	while(i<kv.size() && j<knots.size())
	{
		if(kv[i]<knots[j])
		{
			kv1.push_back(kv[i]);
			i++;
		}
		else
		{
			kv1.push_back(knots[j]);
			j++;
		}
	}
	while(i<kv.size())
	{
		kv1.push_back(kv[i]);
		i++;
	}
	while(j<knots.size())
	{
		kv1.push_back(knots[j]);
		j++;
	}
}

void BezierInsertKnots(vector<double>& kv, std::array<double,2>& knots,vector<double>& kv1)
{
	kv1.clear();
	unsigned int i;
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]<knots[0])
		{
			kv1.push_back(kv[i]);
		}
	}
	for(i=0;i<3;i++)
	{
		kv1.push_back(knots[0]);
	}
	if(knots[0]==kv.front())
	{
		kv1.push_back(knots[0]);
	}
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]>knots[0] && kv[i]<knots[1])
		{
			kv1.push_back(kv[i]);
		}
	}
	for(i=0;i<3;i++)
	{
		kv1.push_back(knots[1]);
	}
	if(knots[1]==kv.back())
	{
		kv1.push_back(knots[1]);
	}
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]>knots[1])
		{
			kv1.push_back(kv[i]);
		}
	}
}

void BezierInsertKnots4(vector<double>& kv, std::array<double,2>& knots,vector<double>& kv1)
{
	kv1.clear();
	unsigned int i;
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]<knots[0])
		{
			kv1.push_back(kv[i]);
		}
	}
	for(i=0;i<4;i++)
	{
		kv1.push_back(knots[0]);
	}
	if(knots[0]==kv.front())
	{
		kv1.push_back(knots[0]);
	}
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]>knots[0] && kv[i]<knots[1])
		{
			kv1.push_back(kv[i]);
		}
	}
	for(i=0;i<4;i++)
	{
		kv1.push_back(knots[1]);
	}
	if(knots[1]==kv.back())
	{
		kv1.push_back(knots[1]);
	}
	for(i=0;i<kv.size();i++)
	{
		if(kv[i]>knots[1])
		{
			kv1.push_back(kv[i]);
		}
	}
}

void InsertKnotsC1(const std::array<double,5>& kv, vector<double>& kv1)
{
	kv1.clear();
	kv1.push_back(kv[0]);
	const double tol(1.e-8);
	for (int i = 1; i < 4; i++)
	{
		kv1.push_back(kv[i]);
		if (fabs(kv[i] - kv[i + 1]) > tol && fabs(kv[i] - kv[i - 1]) > tol)
		{
			kv1.push_back(kv[i]);
		}
	}
	kv1.push_back(kv[4]);
}

void DegreeElevate(vector<vector<double> >& demat)
{
	demat.clear();
	int p(3);
	demat.resize(25,vector<double>(16,0.));
	int loc1(0);
	for(int j=0; j<5; j++)
	{
		for(int i=0; i<5; i++)
		{
			double a(double(i)/double(p+1)), b(double(j)/double(p+1));
			double coef[4]={(1.-a)*(1.-b),a*(1.-b),(1.-a)*b,a*b};
			int loc0[4]={4*j+i,4*j+i-1,4*(j-1)+i,4*(j-1)+i-1};
			if(i==0)
			{
				loc0[1]=-1; loc0[3]=-1;
			}
			if(j==0)
			{
				loc0[2]=-1; loc0[3]=-1;
			}
			if(i==4)
			{
				loc0[0]=-1; loc0[2]=-1;
			}
			if(j==4)
			{
				loc0[0]=-1; loc0[1]=-1;
			}
			for(int k=0; k<4; k++)
			{
				if(loc0[k]!=-1)
				{
					demat[loc1][loc0[k]]=coef[k];
				}
			}
			loc1++;
		}
	}
}

void BezierRefineBi3(const vector<std::array<double, 3> >& bzpt0, vector<std::array<double, 3> >& bzpt1)
{
	if (bzpt0.size() != 16) return;
	bzpt1.clear();
	bzpt1.resize(49);
	double smat1d[7][4] = { { 1.,0.,0.,0. },{ .5,.5,0.,0. },{ .25,.5,.25,0. },{ .125,.375,.375,.125 },
	{ 0.,.25,.5,.25 },{ 0.,0.,.5,.5 },{ 0.,0.,0.,1. } };
	double ctmp;
	int loc(0);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			bzpt1[loc][0] = 0.; bzpt1[loc][1] = 0.; bzpt1[loc][2] = 0.;
			int loc0(0);
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int j0 = 0; j0 < 4; j0++)
				{
					ctmp = smat1d[j][j0] * smat1d[i][i0];
					bzpt1[loc][0] += ctmp*bzpt0[loc0][0];
					bzpt1[loc][1] += ctmp*bzpt0[loc0][1];
					bzpt1[loc][2] += ctmp*bzpt0[loc0][2];
					loc0++;
				}
			}
			loc++;
		}
	}
}