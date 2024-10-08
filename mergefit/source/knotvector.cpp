#include "knotvector.h"

void knotvector::setvec(Standard_Integer Degree, vector<Standard_Real> knotvector)
{
	knotDegree = Degree;

	vector<Standard_Real> vector_of_knots;


	vector<int> vector_of_multiplicity;
	int num = 0;
	for (int i = 0; i < knotvector.size(); i = i + num)
	{
		num = count(knotvector.begin() + num, knotvector.end(), knotvector[i]);
		vector_of_multiplicity.push_back(num);
	}

	knotmultiplicity.Resize(1, vector_of_multiplicity.size(), true);

	for (int i = 0; i < vector_of_multiplicity.size(); i++)
	{
		knotmultiplicity.SetValue(i + 1, vector_of_multiplicity[i]);
	}


	auto it = unique(knotvector.begin(), knotvector.end());
	knotvector.erase(it, knotvector.end());
	vector_of_knots = knotvector;//[0 0.5 1] 

	knots.Resize(1, vector_of_knots.size(), true);

	for (Standard_Integer i = 0; i < vector_of_knots.size(); i++)
	{
		knots.SetValue(i + 1, vector_of_knots[i]);
	}
}