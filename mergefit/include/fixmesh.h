#pragma execution_character_set("utf-8")

#ifndef __FIXMESH_H
#define __FIXMESH_H
#include <stdio.h>
#include "load.h"

using namespace std;



class fIxmesh {
public:

	void Fixmesh(load& mRes, int local_layer, int iteration);

};

#endif
