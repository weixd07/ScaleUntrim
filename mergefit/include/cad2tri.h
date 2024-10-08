#ifndef CAD2TRI_H
#define CAD2TRI_H
#include <TopoDS_Shape.hxx>
double compute_area(const char* in_file);
void triangular(const char* input_file, const char* output_file, double ratio = 0.0001);
void convert(const char* conv_in, const char* conv_out);
#endif
