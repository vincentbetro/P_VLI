#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "geometry.h"
#include "Point.h"
#include "Vector.h"
#include "Edge.h"
#include "Util.h"
#include "viscous_insert.h"
#include "smooth.h"
#include "Linked_List.h"
#include "Bnormal.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define CHECKPOINT(file) {fprintf((file),"\nCheckpoint, file %s, line %d\n",__FILE__, __LINE__); fflush((file));}

#ifdef _DEBUG
    int dbg = 1;
#else
    int dbg = 0;
#endif

#define MAXSWEEP 200

class SW_edge : public Bar_2
{
 public:
  int boundary;
};

class FACE
{
 public:
  int node[4];
};

class BN_FLIST
{
 public:
  int nf;
  FACE *face;
};

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/
extern int displace;

int mesh_deform(geometry *geom, char sname[], char dname[],
                   int osmoo, int lsmoo, int nsearch, int cnvrg)
{
  int b, c, e, i, j, k, l, lay, m, n, n0, n1, n2, n3, n4, n5, n6, n7, s, p, q, z, r, flag;
  int b0, b1, b2, b3, ne, nfn;
  int orig_nn, old_nn, nn, ntet, npyr, npri, nhex, maxtag, m0, m1, m2, m3;
  int tot_neg, tot_neg_skw, tot_pos_skw, tot_pos;
  int *tag, *bft, *itmp;
  char filename[32];
  double mag, maxds;
  double fac, tfac, geo;
  double tol;
  double relax = 0.25;
  Point cg, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, pi, pj;
  Vector norm, vec, v0, v1, v2, v3;
  SW_edge *edge;
  BN_FLIST *bn_face;
  Bnormal *bf_node;
  const int bdim = 132; //buffer dim
  char buff[bdim]; //buffer
  //P_VLI: declare MPI vars
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim2, *bdim21;
  char **sbuff, **rbuff;
  //P_VLI: gathering vars
  int *glmax = 0;
  int **new_elem = 0;
  int *nelem_owned = 0;
  int *own = 0;
  int temp_node, temp_node2, test;

  Smesh_obj *mesh = new Smesh_obj();
  
  //have to make file name explicit for given processor, since input deck only takes generic file name
  if (num_procs > 1)
  {
    sprintf(buff,"%s",sname);
    char *ptr = strstr(buff,".cgns");
    if (ptr == NULL)
    {
      fprintf(stderr,"\nCGNS suffix <.cgns> not found in file name!");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    } else
    {
      //reset cursor to overwrite extension with new extension and proc num
      *ptr = '\0';
    }
    sprintf(sname,"%s_%d.cgns",buff,my_rank);
  }

  mesh->smooth_io(-1,geom,sname);
  
  // check for displaced boundary nodes or displaced/rotate bodies
  int local = mesh->displaced(geom,dname);

  MPI_Allreduce(&local,&displace,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if (!displace)
    return(1);

  //
  // perform Linear-Elastic smoothing
  //
  int success1;
  int success=0;
  int tgn = 1;
  double geo1 = 1.25;
  double geo2 = 1.0;

  MPI_Barrier(MPI_COMM_WORLD); //sync

  success1=mesh->mesh_move(geom,osmoo,lsmoo,tgn,cnvrg,nsearch,geo1,geo2);

  MPI_Barrier(MPI_COMM_WORLD); //to assure all return ready before reduce

  MPI_Allreduce(&success1,&success,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  if (!success)
  {
    fprintf(stderr,"\nMesh movement failed!\n");
    fflush(stderr);
    sprintf(buff,"failed.cgns");
    if (num_procs > 1)
      sprintf(buff,"failed_%d.cgns",my_rank);
    mesh->smooth_io(1,geom,buff);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
    
  MPI_Barrier(MPI_COMM_WORLD);

  mesh->smooth_io(1,geom,sname);
  
  delete mesh;

  return(1);
}
