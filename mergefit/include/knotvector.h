#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <Standard_TypeDef.hxx>

#include <STEPControl_Writer.hxx>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>

#include <BRep_Builder.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_Sewing.hxx>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <TopoDS.hxx>

using namespace std;

class knotvector
{
public:
	Standard_Integer knotDegree;
	TColStd_Array1OfReal knots;//[0 0.5 1]
	TColStd_Array1OfInteger knotmultiplicity;//[000] [111]
	void setvec(Standard_Integer Degree, vector<Standard_Real> knotvector);
};