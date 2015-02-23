#include <stdio.h>
#include <stdlib.h>
#include "Point.h"

#ifndef spoint_obj_h
#define spoint_obj_h

class spoint_obj
{
 public:
  spoint_obj()
  { nn=0; node=0; spacing=0;}
  ~spoint_obj()
  {
    delete[] node;
    delete[] spacing;
    nn=0;
    node=0;
    spacing=0;
  }

  int read_spacing(char *fname);
  int write_spacing(char *fname);

  int nn;
  Point *node;
  double *spacing;
};
#endif

