#pragma once
#include "GenerateNurbs.h"
#include <gp_Ax2.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopExp.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <TopExp_Explorer.hxx>

class GenerateShell
{
public:
	BRepBuilderAPI_Sewing sewing;
	BRepBuilderAPI_Sewing sewing_mirror;
	void Shell(Standard_CString a);
	void Read(const char* filename);
	void MirrorTransformation(const char* a);
	Standard_Real getArea(TopoDS_Shape myshell);
private:
	Standard_Integer patchnumber_is;
	TopoDS_Face tempFace;
};