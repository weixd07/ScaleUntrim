#ifndef KNOT_INSERTION_H
#define KNOT_INSERTION_H

#include <vector>
#include <array>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

void InitialT(const vector<double>& k1,const vector<double>& k2,vector<vector<double> >& T);
void KnotInsert(const vector<double>& k1,const vector<double>& k2,const vector<vector<double> >& T1,int q,vector<vector<double> >& T2);
void TMatrix(const vector<double>& k1,const vector<double>& k2,int p,vector<vector<double> >& T);

void InitialT(const vector<double>& k1, const vector<double>& k2, SparseMatrix<double>& T);
void KnotInsert(const vector<double>& k1, const vector<double>& k2, const SparseMatrix<double>& T1, int q, SparseMatrix<double>& T2);
void TMatrix(const vector<double>& k1, const vector<double>& k2, int p, SparseMatrix<double>& T);

void InsertKnots(const vector<double>& kv,const vector<double>& knots,vector<double>& kv1);
void BezierInsertKnots(vector<double>& kv,std::array<double,2>& knots,vector<double>& kv1);
void BezierInsertKnots4(vector<double>& kv, std::array<double,2>& knots,vector<double>& kv1);

void InsertKnotsC1(const std::array<double,5>& kv, vector<double>& kv1);

void DegreeElevate(vector<vector<double> >& demat);//from degree-3 to degree-4

void BezierRefineBi3(const vector<std::array<double,3> >& bzpt0, vector<std::array<double, 3> >& bzpt1);

#endif