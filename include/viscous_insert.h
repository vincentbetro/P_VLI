#include <stdio.h>
#include "geometry.h"

#ifndef viscous_insert_h
#define viscous_insert_h

int viscous_insert(int v_layers, geometry *geom, char sname[],
                   int osmoo1, int osmoo2, int lsmoo, int nsmoo, int nsearch,
                   double geo1, double geo2, double geo_min, double geo_max, double aao_fac,
                   int cflag, int restart);
int mesh_deform(geometry *geom, char sname[], char dname[],
                   int osmoo, int lsmoo, int nsearch, int cnvrg);

#endif
