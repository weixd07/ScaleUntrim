#include <iostream>
#include "merge_fit.h"

int main(int argc, char *argv[])
{
    std::string cad_in, cad_out, config_file;

    MergeFit mf;
    mf.set_surface_degree(3);
    mf.run(cad_in, cad_out, config_file);
    return 0;
}
