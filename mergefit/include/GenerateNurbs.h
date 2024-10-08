#pragma once
#include "controlpoints.h"
#include "knotvector.h"

class GenerateNurbs
{
public:
	controlpoints cp;
	knotvector uvec;
	knotvector vvec;
	TopoDS_Face get();
};

