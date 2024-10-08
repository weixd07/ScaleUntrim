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

#include <fstream>
#include <string>
#include <sstream>


using namespace std;

class controlpoints
{
public:
	TColgp_Array2OfPnt controlpoints_matrix;
	TColStd_Array2OfReal weights_matrix;
};

