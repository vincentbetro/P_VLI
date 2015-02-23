#include "List.h"

#ifndef Geom_Facet_h
#define Geom_Facet_h
class Geom_Facet
{
 public:

  int nodes[3];
  int mat;
  double space, n_space;
  List *facets;
};
#endif
