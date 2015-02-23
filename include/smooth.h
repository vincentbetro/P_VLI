#include <stdio.h>
#include "Point.h"
#include "Vector.h"
#include "geometry.h"
#include "List.h"
#include "Edge.h"
#include "Bnormal.h"
#include "Chain_node.h"

#ifndef Smesh_obj_h
#define Smesh_obj_h

class EDGE : public Bar_2
{
 public:

  EDGE()
  {
    for (int i=0; i < 6; i++)
      RT[i] = 0.0;
  }

  double RT[6];
};

class Smesh_obj
{
 public:
  Smesh_obj()
  { nn=ntet=npyr=npri=nhex=nb=nc=nv=0;
    nodec=0;
    nodep=0;
    t_n=0;
    q_n=0;
    tet_n=0;
    pyr_n=0;
    pri_n=0;
    hex_n=0;
    nt=0;
    nq=0;
    bname=0;
    pmap = 0;
    tri_map = 0;
    quad_map = 0;
    tet_map = 0;
    pyr_map = 0;
    pri_map = 0;
    hex_map = 0;
  }
  ~Smesh_obj()
  {
    int b,i;
    for (b=0; b < nb; b++)
    {
      if (nt[b] > 0)
      {
        for (i=0; i < nt[b]; i++)
          free(t_n[b][i]);
        free(t_n[b]);
      }
      if (nq[b] > 0)
      {
        for (i=0; i < nq[b]; i++)
          free(q_n[b][i]);
        free(q_n[b]);
      }
	  if (nt[b] > 0 && tri_map != 0)
      {
	    free(tri_map[b]);
      }
	  if (nq[b] > 0 && quad_map != 0)
      {
	    free(quad_map[b]);
      }
	  free(bname[b]);
    }
    free(t_n);
    free(q_n);
    free(bname);
	if (tri_map != 0)
	  free(tri_map);
	if (quad_map != 0)
	  free(quad_map);
	if (tet_map != 0)
	  free(tet_map);
	if (pyr_map != 0)
	  free(pyr_map);
	if (pri_map != 0)
	  free(pri_map);
	if (hex_map != 0)
      free(hex_map);
    free(nt);
    free(nq);
    if (ntet > 0)
    {
      for (i=0; i < ntet; i++)
        free(tet_n[i]);
      free(tet_n);
    }
    if (npyr > 0)
    {
      for (i=0; i < npyr; i++)
        free(pyr_n[i]);
      free(pyr_n);
    }
    if (npri > 0)
    {
      for (i=0; i < npri; i++)
        free(pri_n[i]);
      free(pri_n);
    }
    if (nhex > 0)
    {
      for (i=0; i < nhex; i++)
        free(hex_n[i]);
      free(hex_n);
    }
    free(nodec);
    free(nodep);
    if (pmap != 0)
      {
	for (i = 0; i < nn; i++)
	  free(pmap[i]);
	free(pmap);
      }
    nn=ntet=npyr=npri=nhex=nb=nc=nv=0;
  }
  int check_volumes();
  int check_Jacobians();
  int create_compress_row_storage(int **ia, int **iau, int **ja);
  int create_chain_nodes(int bft[], int tag[], Chain_node **ch_node);
  int create_floating_nodes(int wsmoo, int ftag[], Bnormal **bf_nodes);
  void create_node_to_node(List **nhash);
  int critical_points(geometry *geom, int **cpn);
  int displaced(geometry *geom, char dname[]);
  void elliptic(const int pflag, double u[], double v[], double w[], double Enode[],
                double mat[][3][3], int ia[], int iau[], int ja[],
                List **ltet, List **lpyr, List **lpri, List **lhex);
  void gg_gradients(double f[], double &fx, double &fy, double &fz,
                         List *ltet, List *lpyr, List *lpri, List *lhex);
  int mesh_move(geometry *geom, int osmoo, int lsmoo, int vtag,
                 int cnvrg_order, int nsearch,
                 double geo1, double geo2);
  void output_K(int mdim, double K[][3][3], int ia[], int iau[], int ja[]);
  void output_GMRES(int mdim, double K[][3][3], int ia[], int iau[], int ja[], double u[], double v[], double w[]);
  int read_comp_physical(char sname[]);
  int write_mesh(Point *node, char sname[]);
  void mesh_stats();
  int smooth_io(int mode, geometry *geom, char sname[]);
  void exchange_tag(int *tag, int min_flag);
  void exchange_double(double *fun);
  void exchange_Vector(Vector *vec);

  int nn, ntet, npyr, npri, nhex, nb, nc, nv;
  Point *nodec;
  Point *nodep;
  int ***t_n;
  int ***q_n;
  int **tet_n;
  int **pyr_n;
  int **pri_n;
  int **hex_n;
  int *nt;
  int *nq;
  char **bname;
  int **pmap;
  //these are not needed in P_OPT, but must be generated to be in CGNS file for use with VLI
  int **tri_map, **quad_map, *tet_map, *pyr_map, *pri_map, *hex_map;
};
#endif

