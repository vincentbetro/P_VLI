#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "geometry.h"
#include "smooth.h"
#include "Util.h"
#include "List.h"
#include "nbtri.h"
#include "gmres.h"
#include "Bnormal.h"
#include "Chain_node.h"
#include "svdcmp.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

double pi = 4.0*atan(1.0);

extern double kappa;
extern int E_mode;

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

double Youngs_Modulus(Point p1, Point p2, Point p3, Point p4, int type)
{
  double E=1.0;
  double cn=1.0;
  double vol;
  double wgt[3][3];

  if (E_mode == 1 || E_mode == 3)
  {
    switch(type)
    {
      case 0:
        cn = tetrahedral_condition_number(p1,p2,p3,p4,1);
        break;
      case 1:
        wgt[0][0] = 1.0;
        wgt[1][0] = 0.0;
        wgt[2][0] = 0.0;
        wgt[0][1] = 0.0;
        wgt[1][1] = 1.0;
        wgt[2][1] = 0.0;
        wgt[0][2] = 0.5;
        wgt[1][2] = 0.5;
        //wgt[2][2] = 0.5;
        wgt[2][2] = 1.0/sqrt(2.0);
        cn = condition_number(p1,p2,p3,p4,wgt);
        break;
      case 2:
        wgt[0][0] = 1.0;
        wgt[1][0] = 0.0;
        wgt[2][0] = 0.0;
        wgt[0][1] = 0.5;
        wgt[1][1] = 0.5*sqrt(3);
        wgt[2][1] = 0.0;
        wgt[0][2] = 0.0;
        wgt[1][2] = 0.0;
        wgt[2][2] = 1.0;
        cn = condition_number(p1,p2,p3,p4,wgt);
        break;
      case 3:
        cn = tetrahedral_condition_number(p1,p2,p3,p4,0);
        break;
    }
  }

  double tol = 1.0e-20;
  switch(E_mode)
  {
    case 1:
      E = cn;
      break;
    case 2:
      vol = fabs(tetrahedral_volume(p1,p2,p3,p4));
      E = 1.0/MAX(tol,vol);
      break;
    case 3:
      vol = fabs(tetrahedral_volume(p1,p2,p3,p4));
      E = cn/MAX(tol,vol);
      break;
    default:
      E = 1.0;
      break;
  }

  E = MIN(E,1.0e20);

  return(E);
}

double Poissons_Ratio(Point p1, Point p2, Point p3, Point p4)
{
  double nu = 0.25;

  //nu = 0.0;

  //double ar = tetrahedral_aspect_ratio(p1,p2,p3,p4);
  //nu = MAX(0.001,MIN(0.499,0.5*(1.0-1.0/ar)));

  return(nu);
}

double metric_length(Vector v, double RT[3][3])
{
  int i;
  double mag, tmp[3];

  tmp[0]=tmp[1]=tmp[2]=0.0;
  for (i=0; i < 3; i++)
    tmp[i]=v[0]*RT[0][i]+v[1]*RT[1][i]+v[2]*RT[2][i];
  mag = tmp[0]*v[0]+tmp[1]*v[1]+tmp[2]*v[2];

  if (mag < 1.0e-15)
    return(1.0);
  else
    return(sqrt(mag));
}

double average_metric_length(Vector v, double RT1[3][3], double RT2[3][3])
{
  int i;
  double mag1, mag2, tmp[3];

  tmp[0]=tmp[1]=tmp[2]=0.0;
  for (i=0; i < 3; i++)
    tmp[i]=v[0]*RT1[0][i]+v[1]*RT1[1][i]+v[2]*RT1[2][i];
  mag1 = tmp[0]*v[0]+tmp[1]*v[1]+tmp[2]*v[2];

  tmp[0]=tmp[1]=tmp[2]=0.0;
  for (i=0; i < 3; i++)
    tmp[i]=v[0]*RT2[0][i]+v[1]*RT2[1][i]+v[2]*RT2[2][i];
  mag2 = tmp[0]*v[0]+tmp[1]*v[1]+tmp[2]*v[2];

  if (mag1 < 1.0e-15 || mag2 < 1.0e-15)
    return(1.0);
  else
    return(0.5*(sqrt(mag1)+sqrt(mag2)));
}

void Smesh_obj::create_node_to_node(List **nhash)
{
  int c, n, n0, n1, n2, n3, n4, n5, n6, n7;

  for (n=0; n < nn; n++)
    nhash[n]->Redimension(0);

  for (c=0; c < ntet; c++)
  {
    n0 = tet_n[c][0];
    n1 = tet_n[c][1];
    n2 = tet_n[c][2];
    n3 = tet_n[c][3];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n2);
    nhash[n0]->Check_List(n3);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n3);
    nhash[n2]->Check_List(n0);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n1);
    nhash[n3]->Check_List(n2);
  }
  for (c=0; c < npyr; c++)
  {
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n4);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n4);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n2);
    nhash[n4]->Check_List(n3);
  }
  for (c=0; c < npri; c++)
  {
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n2);
    nhash[n0]->Check_List(n3);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n0);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n5);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n4);
    nhash[n3]->Check_List(n5);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n3);
    nhash[n4]->Check_List(n5);
    nhash[n5]->Check_List(n2);
    nhash[n5]->Check_List(n3);
    nhash[n5]->Check_List(n4);
  }
  for (c=0; c < nhex; c++)
  {
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n5);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n6);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n7);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n5);
    nhash[n4]->Check_List(n7);
    nhash[n5]->Check_List(n1);
    nhash[n5]->Check_List(n4);
    nhash[n5]->Check_List(n6);
    nhash[n6]->Check_List(n2);
    nhash[n6]->Check_List(n5);
    nhash[n6]->Check_List(n7);
    nhash[n7]->Check_List(n3);
    nhash[n7]->Check_List(n4);
    nhash[n7]->Check_List(n6);
  }

}

int Smesh_obj::critical_points(geometry *geom, int **cpn)
{
  int *tag;
  int g, i, j, k, n, ncp, n0, n1, n2, n3;

  tag = new int[nn];

  // identify critical points in mesh
  ncp = 0;
  if (nb == geom->ngb)
  {
    ncp = geom->ncp;
    *cpn = new int[ncp];
    for (k=0; k < ncp; k++)
    {
      for (n=0; n < nn; n++)
        tag[n] = 0;
      for (j=0; j < geom->c_point[k].nmat; j++)
      {
        g = geom->c_point[k].mat[j];
        for (i=0; i < nt[g]; i++)
        {
          n0 = t_n[g][i][0];
          n1 = t_n[g][i][1];
          n2 = t_n[g][i][2];
          if (tag[n0] == j)
            tag[n0] = j+1;
          if (tag[n1] == j)
            tag[n1] = j+1;
          if (tag[n2] == j)
            tag[n2] = j+1;
        }
        for (i=0; i < nq[g]; i++)
        {
          n0 = q_n[g][i][0];
          n1 = q_n[g][i][1];
          n2 = q_n[g][i][2];
          n3 = q_n[g][i][3];
          if (tag[n0] == j)
            tag[n0] = j+1;
          if (tag[n1] == j)
            tag[n1] = j+1;
          if (tag[n2] == j)
            tag[n2] = j+1;
          if (tag[n3] == j)
            tag[n3] = j+1;
        }
      }

      double dsmn = 1.0e20;
      i = -1;
      for (n=0; n < nn; n++)
      {
        if (tag[n] == geom->c_point[k].nmat)
        {
          double ds = distance(geom->g_vert[geom->c_point[k].index],nodep[n]);
          if (ds < dsmn)
          {
            dsmn = ds;
            i = n;
          }
        }
      }
      (*cpn)[k] = i;
    }
  }

  delete[] tag;

  return(ncp);
}

int Smesh_obj::create_chain_nodes(int bft[], int tag[], Chain_node **ch_node)
{
  int b, i, j, m, n, n0, n1, n2, n3;
  int ncn;
  int *itmp;

  itmp = new int[nn];

  ncn = 0;
  (*ch_node) = 0;

  // set critical floating node to line up with chain
  for (b=0; b < nb-1; b++)
  {
    if (!bft[b])
      continue;

    for (n=0; n < nn; n++)
      itmp[n] = 0;
    for (i=0; i < nt[b]; i++)
    {
      n0 = t_n[b][i][0];
      n1 = t_n[b][i][1];
      n2 = t_n[b][i][2];
      itmp[n0] = b+1;
      itmp[n1] = b+1;
      itmp[n2] = b+1;
    }
    for (i=0; i < nq[b]; i++)
    {
      n0 = q_n[b][i][0];
      n1 = q_n[b][i][1];
      n2 = q_n[b][i][2];
      n3 = q_n[b][i][3];
      itmp[n0] = b+1;
      itmp[n1] = b+1;
      itmp[n2] = b+1;
      itmp[n3] = b+1;
    }
    for (j=b+1; j < nb; j++)
    {
      if (!bft[j])
        continue;
      for (i=0; i < nt[j]; i++)
      {
        n0 = t_n[j][i][0];
        n1 = t_n[j][i][1];
        n2 = t_n[j][i][2];
        if (itmp[n0] == b+1) itmp[n0] = j+1;
        if (itmp[n1] == b+1) itmp[n1] = j+1;
        if (itmp[n2] == b+1) itmp[n2] = j+1;
      }
      for (i=0; i < nq[j]; i++)
      {
        n0 = q_n[j][i][0];
        n1 = q_n[j][i][1];
        n2 = q_n[j][i][2];
        n3 = q_n[j][i][3];
        if (itmp[n0] == b+1) itmp[n0] = j+1;
        if (itmp[n1] == b+1) itmp[n1] = j+1;
        if (itmp[n2] == b+1) itmp[n2] = j+1;
        if (itmp[n3] == b+1) itmp[n3] = j+1;
      }
      for (n=0; n < nn; n++)
      {
        if (tag[n] < 0 || itmp[n] != j+1)
          continue;
        if (ncn==0)
          (*ch_node) = (Chain_node*)malloc(sizeof(Chain_node));
        else
          (*ch_node) = (Chain_node*)realloc((void*)(*ch_node),(ncn+1)*sizeof(Chain_node));

        (*ch_node)[ncn].index = n;
        (*ch_node)[ncn].b1 = b;
        (*ch_node)[ncn].b2 = j;
        ncn++;
      }
    }
  }
  delete[] itmp;

  b = ncn;
#ifdef PARALLEL
  MPI_Allreduce(&ncn,&b,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of floating adjacent boundary chain nodes = %d\n",b);
    fflush(out_f);
  }

  return(ncn);
}

int Smesh_obj::create_floating_nodes(int wsmoo, int ftag[], Bnormal **bf_node)
{
  // create floating boundary node structure array
  int g, i, n, nfn, n0, n1, n2, n3;

  int *itmp = new int[nn];
  nfn=0;
  for (g=0; g < nb; g++)
  {
    for (n=0; n < nn; n++)
      itmp[n] = 0;
    for (i=0; i < nt[g]; i++)
    {
      n0 = t_n[g][i][0];
      n1 = t_n[g][i][1];
      n2 = t_n[g][i][2];
      if (wsmoo < 0 || ftag[n0] == 1) itmp[n0] = 1;
      if (wsmoo < 0 || ftag[n1] == 1) itmp[n1] = 1;
      if (wsmoo < 0 || ftag[n2] == 1) itmp[n2] = 1;
    }
    for (i=0; i < nq[g]; i++)
    {
      n0 = q_n[g][i][0];
      n1 = q_n[g][i][1];
      n2 = q_n[g][i][2];
      n3 = q_n[g][i][3];
      if (wsmoo < 0 || ftag[n0] == 1) itmp[n0] = 1;
      if (wsmoo < 0 || ftag[n1] == 1) itmp[n1] = 1;
      if (wsmoo < 0 || ftag[n2] == 1) itmp[n2] = 1;
      if (wsmoo < 0 || ftag[n3] == 1) itmp[n3] = 1;
    }
    for (n=0; n < nn; n++)
      if (itmp[n] == 1) nfn++;
  }
  g = nfn;
#ifdef PARALLEL
  MPI_Allreduce(&nfn,&g,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of floating adjacent boundary nodes = %d\n",g);
    fflush(out_f);
  }

  if (nfn > 0)
  {
    *bf_node = new Bnormal[nfn];
    nfn = 0;
    for (g=0; g < nb; g++)
    {
      for (n=0; n < nn; n++)
        itmp[n] = 0;
      for (i=0; i < nt[g]; i++)
      {
        n0 = t_n[g][i][0];
        n1 = t_n[g][i][1];
        n2 = t_n[g][i][2];
        if (wsmoo < 0 || ftag[n0] == 1) itmp[n0] = 1;
        if (wsmoo < 0 || ftag[n1] == 1) itmp[n1] = 1;
        if (wsmoo < 0 || ftag[n2] == 1) itmp[n2] = 1;
      }
      for (i=0; i < nq[g]; i++)
      {
        n0 = q_n[g][i][0];
        n1 = q_n[g][i][1];
        n2 = q_n[g][i][2];
        n3 = q_n[g][i][3];
        if (wsmoo < 0 || ftag[n0] == 1) itmp[n0] = 1;
        if (wsmoo < 0 || ftag[n1] == 1) itmp[n1] = 1;
        if (wsmoo < 0 || ftag[n2] == 1) itmp[n2] = 1;
        if (wsmoo < 0 || ftag[n3] == 1) itmp[n3] = 1;
      }
      for (n=0; n < nn; n++)
      {
        if (itmp[n] == 1)
        {
          (*bf_node)[nfn].index = n;
          (*bf_node)[nfn].boundary = g;
          (*bf_node)[nfn].norm = Vector(0.0,0.0,0.0);
          nfn++;
        }
      }
    }
  }

  delete[] itmp;

  return(nfn);
}

int Smesh_obj::create_compress_row_storage(int **ia, int **iau, int **ja)
{
  List **nhash;
  int i, j, n;
  int mdim=0;

  // create node-to-node hash list
  nhash = new List*[nn];
  for (n=0; n < nn; n++)
    nhash[n] = new List();

  create_node_to_node(nhash);

  // create compressed row storage pointers
  *ia = new int[nn+1];
  for (n=0; n < nn; n++)
  {
    nhash[n]->Check_List(n);
    (*ia)[n] = mdim;
    mdim += nhash[n]->max;
  }
  (*ia)[nn] = mdim;
  *iau = new int[nn];
  *ja = new int[mdim];
  for (mdim=n=0; n < nn; n++)
  {
    for (i=0; i < nhash[n]->max; i++)
    {
      j = nhash[n]->list[i];
      if (j == n)
        (*iau)[n]=mdim;
      (*ja)[mdim++] = j;
    }
  }

  for (n=0; n < nn; n++)
    delete nhash[n];
  delete[] nhash;

  return(mdim);
}

int Smesh_obj::mesh_move(geometry *geom, int osmoo, int lsmoo, int vtag,
                          int cnvrg_order, int nsearch,
                          double geo1, double geo2)
{
  int b, c, g, i, ir, jc, j, k, l, m, n, sweep, p, rec;
  int n0, n1, n2, n3, n4, n5, n6, n7, ncp, nfn, ncn, ne;
  int start, end, inc, nsrch, gint, gdouble, source, dest, ptag;
  int *cpn;
  double *u, *v, *w;
  double dx, dy, dz, ds, dsmn, dsmx, mag, Lnorm, res, relax;
  List **ltet, **lpyr, **lpri, **lhex;
  Bnormal *bf_node;
  Chain_node *ch_node;
  //internal vars
  int *tag, *ftag, *bt, *ia, *iau, *ja, mdim, order;
  double (*mat)[3][3];
  double (*matinv)[3][3];
  double *wgt, *Enode;
  double (*sol)[3], (*BB)[3];
  double rhs[3], dot, globalLnorm;
  Point *dpsave;
  Point p0, p1, p2, p3, p4, p5, p6, p7, cp, cg, cf, del, pmn, pmx;
  Vector v1, v2, v3, norm;
  bool point_implicit;
  //use struct to get MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmx_g_in, cmx_g_out; 

  int success = 1;

#ifdef PARALLEL
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (my_rank == 0) fprintf(out_f,"\nBegin smoothing mesh.\n");
  if (my_rank == 0) fflush(out_f);

  relax = 0.75;

  // identify critical points in mesh
  ncp = critical_points(geom,&cpn);
  
  tag = new int[nn];

  // set tag to 0 for interior
  // set to -1 for boundary where layers flag is set to vtag
  //           and for non-adjacent boundaries
  // set ftag to 1 for floating adjacent boundary nodes & where layers tag is negative
  //            -1 for boundaries where layers flag is set to vtag
  //             0 for all others
  bt = new int[nb];
  ftag = new int[nn];

  for (n=0; n < nn; n++)
  {
    tag[n] = 0;
    ftag[n] = 0;
  }
  for (b=0; b < nb; b++)
  {
    bt[b] = 1;
    if (vtag > 0 && geom->layers[b] != vtag)
      continue;
    if (vtag <= 0 && geom->layers[b] <= 0)
      continue;
    bt[b] = 0;
    if (my_rank == 0) fprintf(out_f,"\nBody %d is marked viscous.", b+1);
    if (my_rank == 0) fflush(out_f);
    for (i=0; i < nt[b]; i++)
      for (j=0; j < 3; j++)
      {
        n = t_n[b][i][j];
        tag[n] = -1;
        ftag[n] = -1;
      }
    for (i=0; i < nq[b]; i++)
      for (j=0; j < 4; j++)
      {
        n = q_n[b][i][j];
        tag[n] = -1;
        ftag[n] = -1;
      }
  }
  int gfloat;
  for (b=0; b < nb; b++)
  {
    if (bt[b] && geom->layers[b] >= 0)
    {
      k = 0;
      for (i=0; i < nt[b] && !k; i++)
        for (j=0; j < 3 && !k; j++)
        {
          n = t_n[b][i][j];
          if (tag[n] == -1)
            k = 1;
        }
      for (i=0; i < nq[b] && !k; i++)
        for (j=0; j < 4 && !k; j++)
        {
          n = q_n[b][i][j];
          if (tag[n] == -1)
            k = 1;
        }
      if (!k)
        bt[b] = 0;
    }
    if (bt[b] && geom->layers[b] < 0)
      bt[b] = 0;
    
    gfloat = bt[b];
#ifdef PARALLEL
    MPI_Allreduce(&bt[b],&gfloat,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
    bt[b] = gfloat;
  }
  // set tag for floating boundary nodes
  for (b=0; b < nb; b++)
  {
    if (my_rank == 0) fprintf(out_f,"\nBoundary %d, layers flag = %d, boundary float tag = %d",b+1,geom->layers[b],bt[b]);
    if (my_rank == 0) fflush(out_f);
    if (bt[b])
    {
      for (i=0; i < nt[b]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[b][i][j];
          if (tag[n] == 0)
            ftag[n] = 1;
        }
      }
      for (i=0; i < nq[b]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[b][i][j];
          if (tag[n] == 0)
            ftag[n] = 1;
        }
      }
    }
  }
  // turn off movement of non-adjacent boundary points and ghost nodes (-2)
  for (b=0; b < nb; b++)
  {
    if (!bt[b])
    {
      for (i=0; i < nt[b]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[b][i][j];
          tag[n] = -1;
        }
      }
      for (i=0; i < nq[b]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[b][i][j];
          tag[n] = -1;
        }
      }
    }
  }
  //tag ghost nodes
  for (n = 0; n < nn; n++)
    if (pmap[n][1] != my_rank)
      tag[n] = -2;

  order = 3;
  
  // create floating boundary node structure array
  nfn = create_floating_nodes(0,ftag,&bf_node);
  
  // create floating chain node structure array
  ncn = create_chain_nodes(bt,tag,&ch_node);
  
  delete[] ftag;

  mdim=0;
  
  mdim = create_compress_row_storage(&ia, &iau, &ja);

  if (mdim <= 0)
  {
    fprintf(stderr,"\nError creating compressed row storage!");
    fprintf(stderr,"\nTotal dimension of compressed row matrix = %d",mdim);
    fflush(stderr);
#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(0);
  }

  // begin main smoothing loop
  wgt = new double[nn];
  Enode = new double[nn];
  u   = new double[nn];
  v   = new double[nn];
  w   = new double[nn];
  mat = new double[mdim][3][3];
  sol = new double[nn][3];
  BB  = new double[nn][3];
  matinv = new double[nn][3][3];
    

  // create local node-to-element hash lists & global element lists
  ltet = new List*[nn];
  lpyr = new List*[nn];
  lpri = new List*[nn];
  lhex = new List*[nn];
  for (n=0; n < nn; n++)
  {
    ltet[n] = new List();
    lpyr[n] = new List();
    lpri[n] = new List();
    lhex[n] = new List();
  }
  
  for (c=0; c < ntet; c++)
  {
    n0 = tet_n[c][0];
    n1 = tet_n[c][1];
    n2 = tet_n[c][2];
    n3 = tet_n[c][3];
    ltet[n0]->Add_To_List(c);
    ltet[n1]->Add_To_List(c);
    ltet[n2]->Add_To_List(c);
    ltet[n3]->Add_To_List(c);
  }
  for (c=0; c < npyr; c++)
  {
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];
    lpyr[n0]->Add_To_List(c);
    lpyr[n1]->Add_To_List(c);
    lpyr[n2]->Add_To_List(c);
    lpyr[n3]->Add_To_List(c);
    lpyr[n4]->Add_To_List(c);
  }
  for (c=0; c < npri; c++)
  {
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];
    lpri[n0]->Add_To_List(c);
    lpri[n1]->Add_To_List(c);
    lpri[n2]->Add_To_List(c);
    lpri[n3]->Add_To_List(c);
    lpri[n4]->Add_To_List(c);
    lpri[n5]->Add_To_List(c);
  }
  for (c=0; c < nhex; c++)
  {
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];
    lhex[n0]->Add_To_List(c);
    lhex[n1]->Add_To_List(c);
    lhex[n2]->Add_To_List(c);
    lhex[n3]->Add_To_List(c);
    lhex[n4]->Add_To_List(c);
    lhex[n5]->Add_To_List(c);
    lhex[n6]->Add_To_List(c);
    lhex[n7]->Add_To_List(c);
  }
  
  int e, outer, pflag;
  Vector vec1, vec2, vec3, pert, cross1, cross2;
  double vol, cnvrg, icnvrg, diverge, gdsmx;

  int fflag = 0;
  outer = 0;

  for (outer=1; outer <= abs(osmoo) && success; outer++)
  {
    //
    // move boundary nodes and ghost nodes and use linear-elasticity to move interior
    //
    if (success)
    {
      // boundary nodes moving
      // position of boundary nodes is stored in computational coordinates

      //P_VLI:  note that based on the following, ghost nodes are NOT considered in computing max bd movement, but they are used as boundary values
      //also, we will set up the matrix only once per outer loop (using dirichlet conditions for real nodes and ghost nodes (used as bd val))
      //but u,v,w for all interior and ghost (used as bd) nodes will be updated between each gmres iteration for the rhs
      //thus, we have explict bc for real boundaries and implicit-esque bc for ghost nodes used as boundaries, and normal update for interior nodes

      if (outer == 1)
      {
        dpsave = new Point[nn];
        for (n=0; n < nn; n++)
        {
          if (tag[n] < 0)
          {
            Point dp = (nodec[n] - nodep[n]);
            double fac = 0.0;
            for (i=0; i < abs(osmoo); i++)
              fac += pow(geo1,(double)i);
            dpsave[n] = dp/fac;
            //dpsave[n] = (nodec[n] - nodep[n])/abs(osmoo);
          } else
            dpsave[n] = Point(0.0,0.0,0.0);
        }
      }

      // check cell volumes
      success = check_volumes();
      if (success)
        success = check_Jacobians();

#ifdef PARALLEL
      int success1 = success;
 
      MPI_Barrier(MPI_COMM_WORLD); //assure all ready before reduce

      MPI_Allreduce(&success1,&success,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD); //P_VLI: assure that we don't do this until all have finished trying, and one has no success
#endif

      if (success)
      {
        nsrch = nsearch;
        
        //P_VLI:set u,v,w for ghost nodes too
        dsmx = 0.0;
        for (n=0; n < nn; n++)
        {
          u[n] = 0.0;
          v[n] = 0.0;
          w[n] = 0.0;
          wgt[n] = 0.0;
          if (tag[n] < 0)
          {
            u[n] = dpsave[n][0];
            v[n] = dpsave[n][1];
            w[n] = dpsave[n][2];
            dpsave[n] *= geo1;
            if (tag[n] == -2)
              continue; //P_VLI:  do not include ghost "bd" in bd spacing max
            dsmx = MAX(dsmx,u[n]*u[n]+v[n]*v[n]+w[n]*w[n]);
          }
        }

        gdsmx = dsmx = sqrt(dsmx);
#ifdef PARALLEL
        MPI_Allreduce(&dsmx,&gdsmx,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); //P_VLI: assure that we don't do this until all have finished trying, and one has no success
#endif
        if (my_rank == 0) fprintf(out_f,"\nMaximum boundary movement = %lg",gdsmx);
        if (my_rank == 0) fflush(out_f);

        // initialize 3X3 matrices
        for (n=0; n < mdim; n++)
        {
          mat[n][0][0] = 0.0;
          mat[n][0][1] = 0.0;
          mat[n][0][2] = 0.0;
          mat[n][1][0] = 0.0;
          mat[n][1][1] = 0.0;
          mat[n][1][2] = 0.0;
          mat[n][2][0] = 0.0;
          mat[n][2][1] = 0.0;
          mat[n][2][2] = 0.0;
        }

        // compute nodal values for E
        double E, vol;
        for (n=0; n < nn; n++)
          wgt[n] = Enode[n] = 0.0;
        for (c=0; c < ntet; c++)
        {
          p0 = nodep[n0=tet_n[c][0]];
          p1 = nodep[n1=tet_n[c][1]];
          p2 = nodep[n2=tet_n[c][2]];
          p3 = nodep[n3=tet_n[c][3]];
          switch (E_mode)
          {
            case 1: case 5: E = tetrahedral_aspect_ratio(p0,p1,p2,p3); break;
            case 2: case 6: E = 1.0/tetrahedral_volume(p0,p1,p2,p3); break;
            case 3: case 7: E = tetrahedral_aspect_ratio(p0,p1,p2,p3)/tetrahedral_volume(p0,p1,p2,p3); break;
            default: E = 1.0; break;
          }
          Enode[n0] += E;
          Enode[n1] += E;
          Enode[n2] += E;
          Enode[n3] += E;
          wgt[n0] += 1.0;
          wgt[n1] += 1.0;
          wgt[n2] += 1.0;
          wgt[n3] += 1.0;
        }
        for (c=0; c < npyr; c++)
        {
          p0 = nodep[n0=pyr_n[c][0]];
          p1 = nodep[n1=pyr_n[c][1]];
          p2 = nodep[n2=pyr_n[c][2]];
          p3 = nodep[n3=pyr_n[c][3]];
          p4 = nodep[n4=pyr_n[c][4]];
          switch (E_mode)
          {
            case 1: case 5: E = pyramid_aspect_ratio(p0,p1,p2,p3,p4); break;
            case 2: case 6: E = 1.0/pyramid_volume(p0,p1,p2,p3,p4); break;
            case 3: case 7: E = pyramid_aspect_ratio(p0,p1,p2,p3,p4)/pyramid_volume(p0,p1,p2,p3,p4); break;
            default: E = 1.0; break;
          }
          Enode[n0] += E;
          Enode[n1] += E;
          Enode[n2] += E;
          Enode[n3] += E;
          Enode[n4] += E;
          wgt[n0] += 1.0;
          wgt[n1] += 1.0;
          wgt[n2] += 1.0;
          wgt[n3] += 1.0;
          wgt[n4] += 1.0;
        }
        for (c=0; c < npri; c++)
        {
          p0 = nodep[n0=pri_n[c][0]];
          p1 = nodep[n1=pri_n[c][1]];
          p2 = nodep[n2=pri_n[c][2]];
          p3 = nodep[n3=pri_n[c][3]];
          p4 = nodep[n4=pri_n[c][4]];
          p5 = nodep[n5=pri_n[c][5]];
          switch (E_mode)
          {
            case 1: case 5: E = prism_aspect_ratio(p0,p1,p2,p3,p4,p5); break;
            case 2: case 6: E = 1.0/prism_volume(p0,p1,p2,p3,p4,p5); break;
            case 3: case 7: E = prism_aspect_ratio(p0,p1,p2,p3,p4,p5)/prism_volume(p0,p1,p2,p3,p4,p5); break;
            default: E = 1.0; break;
          }
          Enode[n0] += E;
          Enode[n1] += E;
          Enode[n2] += E;
          Enode[n3] += E;
          Enode[n4] += E;
          Enode[n5] += E;
          wgt[n0] += 1.0;
          wgt[n1] += 1.0;
          wgt[n2] += 1.0;
          wgt[n3] += 1.0;
          wgt[n4] += 1.0;
          wgt[n5] += 1.0;
        }
        for (c=0; c < nhex; c++)
        {
          p0 = nodep[n0=hex_n[c][0]];
          p1 = nodep[n1=hex_n[c][1]];
          p2 = nodep[n2=hex_n[c][2]];
          p3 = nodep[n3=hex_n[c][3]];
          p4 = nodep[n4=hex_n[c][4]];
          p5 = nodep[n5=hex_n[c][5]];
          p6 = nodep[n6=hex_n[c][6]];
          p7 = nodep[n7=hex_n[c][7]];
          switch (E_mode)
          {
            case 1: case 5: E = hexahedral_aspect_ratio(p0,p1,p2,p3,p4,p5,p6,p7); break;
            case 2: case 6: E = 1.0/hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7); break;
            case 3: case 7: E = hexahedral_aspect_ratio(p0,p1,p2,p3,p4,p5,p6,p7)/hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7); break;
            default: E = 1.0; break;
          }
          Enode[n0] += E;
          Enode[n1] += E;
          Enode[n2] += E;
          Enode[n3] += E;
          Enode[n4] += E;
          Enode[n5] += E;
          Enode[n6] += E;
          Enode[n7] += E;
          wgt[n0] += 1.0;
          wgt[n1] += 1.0;
          wgt[n2] += 1.0;
          wgt[n3] += 1.0;
          wgt[n4] += 1.0;
          wgt[n5] += 1.0;
          wgt[n6] += 1.0;
          wgt[n7] += 1.0;
        }
        for (n=0; n < nn; n++)
        {
          Enode[n] /= wgt[n];
          wgt[n] = 0.0;
          if (E_mode > 3)
          {
            ds = 1e20;
            cg = nodep[n];
            for (b=0; b < nb; b++)
            {
              g = b;
              if (geom->layers[g] <= 0) continue;
              if (cg[0]+ds < geom->blo[g][0] || cg[1]+ds < geom->blo[g][1] || cg[2]+ds < geom->blo[g][2] ||
                  cg[0]-ds > geom->bhi[g][0] || cg[1]-ds > geom->bhi[g][1] || cg[2]-ds > geom->bhi[g][2])
                continue;
              norm = Vector(0.0,0.0,0.0);
              cf = geom->closest(cg,g,norm);
              ds = MIN(ds,distance(cg,cf));
            }
            Enode[n] /= MAX(1e-10,ds);
          }
        }
        exchange_double(Enode);
        double Emin = 1.0e20;
        double Emax = 0.0;
        for (n=0; n < nn; n++)
        {
          Emin = MIN(Emin,Enode[n]);
          Emax = MAX(Emax,Enode[n]);
        }
#ifdef PARALLEL
        if (num_procs > 1)
        {
          E = Emin;
          MPI_Allreduce(&E,&Emin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
          E = Emax;
          MPI_Allreduce(&E,&Emax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        }
#endif
        E = MIN(999.0,Emax-1.0);
        if (E_mode > 3 && my_rank == 0)
          fprintf(out_f,"\nComputed Young's Modulus range is %lg  -  %lg",Emin,Emax);
        //if (my_rank == 0)
        //  fprintf(out_f,"\nComputed Young's Modulus range is %lg  -  %lg",Emin,Emax);
        //  fprintf(out_f,"\nComputed Young's Modulus range is %lg  -  %lg, rescaling to 1.0  -  %lg",Emin,Emax,E+1.0);
        //for (n=0; n < nn; n++)
        //  Enode[n] = (Enode[n]-Emin)/(Emax-Emin)*E + 1.0;

#ifdef PARALLEL
        exchange_double(Enode);
#endif

        // compute equation coefficients...P_VLI: note that there is no need to distinguish between types here, since boundary values reset before gmres
        pflag = 1;
        elliptic(pflag,u,v,w,Enode,mat,ia,iau,ja,ltet,lpyr,lpri,lhex);

        // logic to handle tangent perturbation boundary condition
        for (i=0; i < nfn; i++)
        {
          n = bf_node[i].index;
          b = bf_node[i].boundary;
          norm = Vector(0.0,0.0,0.0);
          p1 = geom->closest(nodep[n],b,norm);
          bf_node[i].norm = norm;

          // reset all matrices on row for node n
          //for (k=ia[n]; k < ia[n+1]; k++)
          //  for (ir=0; ir < order; ir++)
          //    for (jc=0; jc < order; jc++)
          //      mat[k][ir][jc] = 0.0;

          //m = iau[n];
          //mat[m][0][0] += 1.0;
          //mat[m][1][1] += 1.0;
          //mat[m][2][2] += 1.0;
        }

        // overwrite matrices for prescribed boundary points
        for (n=0; n < nn; n++)
        {
          if (tag[n] != -1 && tag[n] != -2)
            continue;
          for (k=ia[n]; k < ia[n+1]; k++)
          {
            m = ja[k];
            for (ir=0; ir < order; ir++)
            {
              for (jc=0; jc < order; jc++)
              {
                if (ir==jc && m == n)
                  mat[k][ir][jc] = 1.0;
                else
                  mat[k][ir][jc] = 0.0;
              }
            }
          }
        }

        //output_K(mdim,mat,ia,iau,ja);
        //output_GMRES(mdim,mat,ia,iau,ja,u,v,w);

        //cnvrg = -1.0;
        cnvrg = gdsmx/pow(10.0,cnvrg_order);
        if (my_rank == 0) fprintf(out_f,"\nConvergence Lnorm = %lg",cnvrg);
        //diverge = gdsmx*pow(10.0,cnvrg_order);
        diverge = -1.0;
        icnvrg = -1.0;
        point_implicit = false;
        if (nsrch <= 0) point_implicit = true;
        //if (fabs(kappa) > 1.0e-12) point_implicit = true;
        if (my_rank == 0 && point_implicit) fprintf(out_f,"\nUsing point implicit scheme.");
        if (my_rank == 0 && !point_implicit) fprintf(out_f,"\nUsing GMRES scheme.");

        // NOTE: invert diagonal if using point implicit scheme
        if (point_implicit)
        {
          for (n=0; n < nn; n++)
          {
            m = iau[n];
            // precondition matrix NOT RIGHT!!!!!
            //for (ir=0; ir < order; ir++)
            //{
            //  mag = 1.0/mat[m][ir][ir];
            //  for (k=ia[n]; k < ia[n+1]; k++)
            //    for (jc=0; jc < order; jc++)
            //      mat[k][ir][jc] *= mag;
            //}
            // compute inverse of diagonal
            for (ir=0; ir < 3; ir++)
              for (jc=0; jc < 3; jc++)
                matinv[n][ir][jc] = mat[iau[n]][ir][jc];
            invert3X3(matinv[n]);

            //double tmp[3][3];
            //if (!point_implicit)
            //  for (k=ia[n]; k < ia[n+1]; k++)
            //  {
            //    for (ir=0; ir < 3; ir++)
            //    {
            //      for (jc=0; jc < 3; jc++)
            //      {
            //        mag = 0.0;
            //        for (i=0; i < 3; i++)
            //          mag += matinv[n][ir][i]*mat[k][i][jc];
            //        tmp[ir][jc] = mag;
            //      }
            //    }
            //    for (ir=0; ir < 3; ir++)
            //      for (jc=0; jc < 3; jc++)
            //        mat[k][ir][jc] = tmp[ir][jc];
            //  }
          }
        }

        if (my_rank == 0) fprintf(out_f,"\nLinear-Elasticity");
        if (my_rank == 0) fprintf(out_f,"\nouter  sweep #    Lnorm      Max delta");
        if (my_rank == 0) fprintf(out_f,"\n----- -------- ------------ ------------");
        if (my_rank == 0) fflush(out_f);

        sweep = 0;
        do
        {
          sweep++;
          res=Lnorm = 0.0;
          dsmx= 0.0;
          
          for (n=0; n < nn; n++)
          {
            if (tag[n] >= 0)
            {
              //rhs[0]=rhs[1]=rhs[2]=0.0;
              //for (i=ia[n]; i < ia[n+1]; i++)
              //{
              //  m = ja[i];
              //  rhs[0] += (mat[i][0][0]*u[m]+mat[i][0][1]*v[m]+mat[i][0][2]*w[m]);
              //  rhs[1] += (mat[i][1][0]*u[m]+mat[i][1][1]*v[m]+mat[i][1][2]*w[m]);
              //  rhs[2] += (mat[i][2][0]*u[m]+mat[i][2][1]*v[m]+mat[i][2][2]*w[m]);
              //}
              //res += rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2];
              BB[n][0] = BB[n][1] = BB[n][2] = 0.0;
            } else if (tag[n] == -2)
            {
              BB[n][0] = u[n];
              BB[n][1] = v[n];
              BB[n][2] = w[n];
            } else
            {
              BB[n][0] = u[n];
              BB[n][1] = v[n];
              BB[n][2] = w[n];
            }
            sol[n][0] = u[n];
            sol[n][1] = v[n];
            sol[n][2] = w[n];
          }
          //for (i=0; i < nfn; i++)
          //{
          //  n = bf_node[i].index;
          //  BB[n][0] += bf_node[i].norm[0];
          //  BB[n][1] += bf_node[i].norm[1];
          //  BB[n][2] += bf_node[i].norm[2];
          //}

          relax = 1.0;
          if (point_implicit)
          {
            relax = 0.75;
            relax = MIN(0.5,MAX(0.001/abs(lsmoo),1.0/sqrt(Enode[n])));
            for (n=0; n < nn; n++)
            {
              rhs[0]=rhs[1]=rhs[2]=0.0;
              for (i=ia[n]; i < ia[n+1]; i++)
              {
                m = ja[i];
                if (m == n)
                  continue;
                rhs[0] -= (mat[i][0][0]*u[m]+mat[i][0][1]*v[m]+mat[i][0][2]*w[m]);
                rhs[1] -= (mat[i][1][0]*u[m]+mat[i][1][1]*v[m]+mat[i][1][2]*w[m]);
                rhs[2] -= (mat[i][2][0]*u[m]+mat[i][2][1]*v[m]+mat[i][2][2]*w[m]);
              }
              sol[n][0] = matinv[n][0][0]*rhs[0]+matinv[n][0][1]*rhs[1]+matinv[n][0][2]*rhs[2];
              sol[n][1] = matinv[n][1][0]*rhs[0]+matinv[n][1][1]*rhs[1]+matinv[n][1][2]*rhs[2];
              sol[n][2] = matinv[n][2][0]*rhs[0]+matinv[n][2][1]*rhs[1]+matinv[n][2][2]*rhs[2];
            }
          } else
            gmres(nsrch,(double)cnvrg_order,nn,mdim,ia,ja,iau,mat,sol,BB);

          //
          // impose tangential movement of adjacent boundary points along geometry
          //
          for (i=0; i < nfn; i++)
          {
            n = bf_node[i].index;
            pert = Vector(sol[n][0],sol[n][1],sol[n][2]);
            //ds = pert.magnitude();
            dot = pert * bf_node[i].norm;
            pert -= bf_node[i].norm*dot;
            //pert.normalize();
            //pert *= ds;
            sol[n][0] = pert[0];
            sol[n][1] = pert[1];
            sol[n][2] = pert[2];
          }

          // force boundary nodes of adjacent boundaries to lie on geometry chain
          for (i=0; i < ncn; i++)
          {
            n = ch_node[i].index;
            g = ch_node[i].b1;
            b = ch_node[i].b2;
            pert = Vector(sol[n][0],sol[n][1],sol[n][2]);
            //ds = pert.magnitude();
            p1 = Point(nodep[n][0]+sol[n][0],nodep[n][1]+sol[n][1],nodep[n][2]+sol[n][2]);
            p2 = geom->closest_chain(p1,g,b);
            pert = Vector(nodep[n],p2);
            //pert.normalize();
            //pert *= ds;
            sol[n][0] = pert[0];
            sol[n][1] = pert[1];
            sol[n][2] = pert[2];
          }

          for (n=0; n < nn; n++)
          {
            if (tag[n] >= 0)
            {
              del[0]=sol[n][0];
              del[1]=sol[n][1];
              del[2]=sol[n][2];
              del[0] = (del[0]-u[n]);
              del[1] = (del[1]-v[n]);
              del[2] = (del[2]-w[n]);
              ds=(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
              Lnorm = MAX(Lnorm,sqrt(ds));
              u[n] += del[0]*relax;
              v[n] += del[1]*relax;
              w[n] += del[2]*relax;
              ds = u[n]*u[n]+v[n]*v[n]+w[n]*w[n];
              if (ds > dsmx)
              {
                dsmx = ds;
                pmx = nodep[n];
              }
            }
          }

          dsmx = sqrt(dsmx);

          //MPI_Allreduce(&res,&globalLnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
          //res = globalLnorm;

#ifdef PARALLEL
          MPI_Barrier(MPI_COMM_WORLD);
#endif

          globalLnorm = Lnorm;
#ifdef PARALLEL
          MPI_Allreduce(&Lnorm,&globalLnorm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

          if (diverge < 0.0)
            diverge = globalLnorm*pow(10.0,cnvrg_order);

          if (icnvrg < 0.0)
            icnvrg = globalLnorm/100.0;

#ifdef PARALLEL
          cmx_g_in.val = dsmx;
          cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              dsmx = cmx_g_out.val;
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif

          //if (my_rank == 0) fprintf(out_f,"\n%5d %8d %12.6e %12.6e %12.6e @ %d ( %g, %g, %g )",
          //       outer,sweep,res,globalLnorm,cmx_g_out.val,cmx_g_out.proc,pmx[0],pmx[1],pmx[2]);
          if (my_rank == 0) fprintf(out_f,"\n%5d %8d %12.6e %12.6e @ %d ( %g, %g, %g )",
                 outer,sweep,globalLnorm,cmx_g_out.val,cmx_g_out.proc,pmx[0],pmx[1],pmx[2]);
          if (my_rank == 0) fflush(out_f);

#ifdef PARALLEL
          MPI_Barrier(MPI_COMM_WORLD);
          // exchange ghost node values of displacements
          exchange_double(u);
          exchange_double(v);
          exchange_double(w);
#endif

          if (globalLnorm > diverge)
          {
            if (my_rank == 0) fprintf(out_f,"\nSolution diverging!");
            if (my_rank == 0) fflush(out_f);
            success = 0;
          }
          //if (point_implicit && globalLnorm < icnvrg)
          //{
          //  if (my_rank == 0) fprintf(out_f,"\nSwitching to GMRES");
          //  if (my_rank == 0) fflush(out_f);
          //  point_implicit = false;
          //}

        } while (sweep < abs(lsmoo) && globalLnorm > cnvrg && success);

        //
        // impose tangential movement of adjacent boundary points along geometry
        //
        //for (i=0; i < nfn; i++)
        //{
        //  n = bf_node[i].index;
        //  pert = Vector(u[n],v[n],w[n]);
        //  //ds = pert.magnitude();
        //  dot = pert * bf_node[i].norm;
        //  pert -= bf_node[i].norm*dot;
        //  //pert.normalize();
        //  //pert *= ds;
        //  u[n] = pert[0];
        //  v[n] = pert[1];
        //  w[n] = pert[2];
        //}

        //
        // add projection of adjacent boundary points to geometry
        //
        for (n=0; n < nn; n++)
          wgt[n] = 0.0;

        for (g=0; g < nb; g++)
        {
          if (!bt[g])
            continue;
          for (i=0; i < nt[g]; i++)
          {
            n0 = t_n[g][i][0];
            n1 = t_n[g][i][1];
            n2 = t_n[g][i][2];
            wgt[n0] = (g+1);
            wgt[n1] = (g+1);
            wgt[n2] = (g+1);
          }
          for (i=0; i < nq[g]; i++)
          {
            n0 = q_n[g][i][0];
            n1 = q_n[g][i][1];
            n2 = q_n[g][i][2];
            n3 = q_n[g][i][3];
            wgt[n0] = (g+1);
            wgt[n1] = (g+1);
            wgt[n2] = (g+1);
            wgt[n3] = (g+1);
          }

          for (n=0; n < nn; n++)
          {
            if (tag[n] < 0 || fabs(wgt[n]-(g+1)) > 1.e-10)
              continue;
            b = g;
            p1 = Point(nodep[n][0]+u[n],nodep[n][1]+v[n],nodep[n][2]+w[n]);
            vec1 = Vector(0.0,0.0,0.0);
            p2 = geom->closest(p1,b,vec1);
            u[n] = p2[0] - nodep[n][0];
            v[n] = p2[1] - nodep[n][1];
            w[n] = p2[2] - nodep[n][2];
          }
        }

        // force boundary nodes of adjacent boundaries to lie on geometry chain
        for (i=0; i < ncn; i++)
        {
          n = ch_node[i].index;
          g = ch_node[i].b1;
          b = ch_node[i].b2;
          p1 = Point(nodep[n][0]+u[n],nodep[n][1]+v[n],nodep[n][2]+w[n]);
          p2 = geom->closest_chain(p1,g,b);
          u[n] = p2[0] - nodep[n][0];
          v[n] = p2[1] - nodep[n][1];
          w[n] = p2[2] - nodep[n][2];
        }

        for (n=0; n < nn; n++)
        {
          nodep[n][0] += u[n];
          nodep[n][1] += v[n];
          nodep[n][2] += w[n];
        }

      }

      if (outer == abs(osmoo) || success == 0)
        delete[] dpsave;

    }
	
#ifdef _DEBUG
    char fname[128];
    sprintf(fname,"Debug_%d.cgns",my_rank);
    smooth_io(1,geom,fname);
#endif
  }

  if (success)
    success = check_volumes();
#ifdef PARALLEL
  MPI_Allreduce(&success,&i,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  success = i;
#endif
  if (success)
    success = check_Jacobians();
#ifdef PARALLEL
  MPI_Allreduce(&success,&i,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  success = i;
#endif
  
  //mesh_stats();

  // release list memory
  for (n=0; n < nn; n++)
  {
    delete ltet[n];
    delete lpyr[n];
    delete lpri[n];
    delete lhex[n];
  }

  delete[] ltet;
  delete[] lpyr;
  delete[] lpri;
  delete[] lhex;

  if (nfn > 0)
    delete[] bf_node;
  if (ncn > 0)
    free(ch_node);
  delete[] bt;
  delete[] wgt;
  delete[] Enode;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] tag;
  if (sol > 0)
    delete[] sol;
  if (BB > 0)
  delete[] BB;
  delete[] mat;
  delete[] ia;
  delete[] iau;
  delete[] ja;
  delete[] matinv;
  if (nb == geom->ngb)
    delete[] cpn;

  return(success);
}

int Smesh_obj::check_volumes()
{
  int c, i, k, n, n0, n1, n2, n3, n4, n5, n6, n7, success, source, dest, local, tag;
  double vol, vmin, gdouble;
  Point cp, p0, p1, p2, p3, p4, p5, p6, p7;
  int gntet, gnpyr, gnpri, gnhex;
#ifdef PARALLEL
  MPI_Status status;
#endif
  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 

  success = 1;

  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;
  if (num_procs > 1)
  {
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    gntet=local=0;
    for (n=0; n < ntet; n++)
      local = MAX(local,tet_map[n]+1);
    gntet=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gntet,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnpyr=local=0;
    for (n=0; n < npyr; n++)
      local = MAX(local,pyr_map[n]+1);
    gnpyr=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpyr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnpri=local=0;
    for (n=0; n < npri; n++)
      local = MAX(local,pri_map[n]+1);
    gnpri=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpri,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnhex=local=0;
    for (n=0; n < nhex; n++)
      local = MAX(local,hex_map[n]+1);
    gnhex=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnhex,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  }

  // check cell volumes

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gntet > 0)
  {
    for (k=c=0; c < ntet; c++)
    {
      n0 = tet_n[c][0];
      n1 = tet_n[c][1];
      n2 = tet_n[c][2];
      n3 = tet_n[c][3];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      vol = tetrahedral_volume(p0,p1,p2,p3);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          vmin = vol;
          cp = (p0+p1+p2+p3)/4.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative tetrahedral volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpyr > 0)
  {
    for (k=c=0; c < npyr; c++)
    {
      n0 = pyr_n[c][0];
      n1 = pyr_n[c][1];
      n2 = pyr_n[c][2];
      n3 = pyr_n[c][3];
      n4 = pyr_n[c][4];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      vol = pyramid_volume(p0,p1,p2,p3,p4);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          vmin = vol;
          cp = (p0+p1+p2+p3+p4)/5.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative pyramid volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpri > 0)
  {
    for (k=c=0; c < npri; c++)
    {
      n0 = pri_n[c][0];
      n1 = pri_n[c][1];
      n2 = pri_n[c][2];
      n3 = pri_n[c][3];
      n4 = pri_n[c][4];
      n5 = pri_n[c][5];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      vol = prism_volume(p0,p1,p2,p3,p4,p5);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          vmin = vol;
          cp = (p0+p1+p2+p3+p4+p5)/6.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative prism volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnhex > 0)
  {
    for (k=c=0; c < nhex; c++)
    {
      n0 = hex_n[c][0];
      n1 = hex_n[c][1];
      n2 = hex_n[c][2];
      n3 = hex_n[c][3];
      n4 = hex_n[c][4];
      n5 = hex_n[c][5];
      n6 = hex_n[c][6];
      n7 = hex_n[c][7];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      vol = hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          vmin = vol;
          cp = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative hexahedral volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(&success,&i,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  success = i;
#endif
  return(success);
}

int Smesh_obj::check_Jacobians()
{
  int c, i, k, n, n0, n1, n2, n3, n4, n5, n6, n7, success, source, dest, local, tag;
  double jac, jvg, jmn, vmin, gdouble;
  Point cp, p0, p1, p2, p3, p4, p5, p6, p7;
  Vector v1, v2, v3;
  int gntet, gnpyr, gnpri, gnhex;
#ifdef PARALLEL
  MPI_Status status;
#endif
  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 

  success = 1;

  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;
  if (num_procs > 1)
  {
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    gntet=local=0;
    for (n=0; n < ntet; n++)
      local = MAX(local,tet_map[n]+1);
    gntet=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gntet,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnpyr=local=0;
    for (n=0; n < npyr; n++)
      local = MAX(local,pyr_map[n]+1);
    gnpyr=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpyr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnpri=local=0;
    for (n=0; n < npri; n++)
      local = MAX(local,pri_map[n]+1);
    gnpri=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpri,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

    gnhex=local=0;
    for (n=0; n < nhex; n++)
      local = MAX(local,hex_map[n]+1);
    gnhex=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnhex,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  }

  // check corner Jacobians

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gntet > 0)
  {
    for (k=c=0; c < ntet; c++)
    {
      n0 = tet_n[c][0];
      n1 = tet_n[c][1];
      n2 = tet_n[c][2];
      n3 = tet_n[c][3];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p2);
      v3 = Vector(p0,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      if (jac < 1.0e-40)
      {
        if (jac < vmin)
        {
          vmin = jac;
          cp = (p0+p1+p2+p3)/4.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral Jacobian = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative tetrahedral Jacobians = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpyr > 0)
  {
    for (k=c=0; c < npyr; c++)
    {
      n0 = pyr_n[c][0];
      n1 = pyr_n[c][1];
      n2 = pyr_n[c][2];
      n3 = pyr_n[c][3];
      n4 = pyr_n[c][4];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p3);
      v3 = Vector(p0,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jvg = scalar_triple_product(v1,v2,v3);
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      jvg/=4.0;
      if (jvg < 1.0e-40)
      {
        if (jvg < vmin)
        {
          vmin = jvg;
          cp = (p0+p1+p2+p3+p4)/5.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid Jacobian = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative pyramid Jacobians = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpri > 0)
  {
    for (k=c=0; c < npri; c++)
    {
      n0 = pri_n[c][0];
      n1 = pri_n[c][1];
      n2 = pri_n[c][2];
      n3 = pri_n[c][3];
      n4 = pri_n[c][4];
      n5 = pri_n[c][5];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p2);
      v3 = Vector(p0,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jvg = scalar_triple_product(v1,v2,v3);
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p2,p0);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p3,p5);
      v2 = Vector(p3,p4);
      v3 = Vector(p3,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p4,p3);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p3);
      v3 = Vector(p5,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      jvg/=6.0;
      if (jvg < 1.0e-40)
      {
        if (jvg < vmin)
        {
          vmin = jvg;
          cp = (p0+p1+p2+p3+p4+p5)/6.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism Jacobian = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative prism Jacobians = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnhex > 0)
  {
    for (k=c=0; c < nhex; c++)
    {
      n0 = hex_n[c][0];
      n1 = hex_n[c][1];
      n2 = hex_n[c][2];
      n3 = hex_n[c][3];
      n4 = hex_n[c][4];
      n5 = hex_n[c][5];
      n6 = hex_n[c][6];
      n7 = hex_n[c][7];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p3);
      v3 = Vector(p0,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jvg = scalar_triple_product(v1,v2,v3);
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p6);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p7);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p4,p7);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p6);
      v3 = Vector(p5,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p6,p5);
      v2 = Vector(p6,p7);
      v3 = Vector(p6,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      v1 = Vector(p7,p6);
      v2 = Vector(p7,p4);
      v3 = Vector(p7,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      jac = scalar_triple_product(v1,v2,v3);
      jvg+=jac;
      jvg/=8.0;
      if (jvg < 1.0e-40)
      {
        if (jvg < vmin)
        {
          vmin = jvg;
          cp = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
        }
        success=0;
        k++;
      }
    }
    i=k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral Jacobian = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      fprintf(out_f,"\nTotal number of negative hexahedral Jacobians = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(&success,&i,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  success = i;
#endif
  return(success);
}

void Smesh_obj::elliptic(const int pflag, double u[], double v[], double w[], double Enode[],
                         double mat[][3][3], int ia[], int iau[], int ja[],
                         List **ltet, List **lpyr, List **lpri, List **lhex)
{
  int c, e, i, j, k, n, n1, n2, n3, n4, ncells, type, tflag, source, dest, tag;
  Point p1, p2, p3, p4, ci, cc;
  Point q1, q2, q3, q4;
  double vc, rvc, E, nu;
  double t1u1, t1u2, t1u3, t1u4;
  double t1v1, t1v2, t1v3, t1v4;
  double t1w1, t1w2, t1w3, t1w4;
  double t2u1, t2u2, t2u3, t2u4;
  double t2v1, t2v2, t2v3, t2v4;
  double t2w1, t2w2, t2w3, t2w4;
  double t3u1, t3u2, t3u3, t3u4;
  double t3v1, t3v2, t3v3, t3v4;
  double t3w1, t3w2, t3w3, t3w4;
  double a11, a12, a13;
  double a21, a22, a23;
  double a31, a32, a33;
  double t11, t12, t13;
  double t21, t22, t23;
  double t31, t32, t33;
  double ux, uy, uz, vx, vy, vz, wx, wy, wz, gdouble;
  Vector a1, a2, a3, a4, t, vec1, vec2, vec3;
  double tol = 1.0e-20;
  double Emin = 1.0e20;
  double Emax = 0.0;
  double numin = 1.0e20;
  double numax = 0.0;
  Point *node;
  Point pmin,pmax,qmin,qmax;
#ifdef PARALLEL
  MPI_Status status;
#endif
  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 


  node = nodep;

  for (n=0; n < nn; n++)
  {
    // outer loop switches between 4 different cell types
    for (type=0; type < 4; type++)
    {
      switch(type)
      {
        case 0: ncells = ltet[n]->max; break;
        case 1: ncells = lpyr[n]->max; break;
        case 2: ncells = lpri[n]->max; break;
        case 3: ncells = lhex[n]->max; break;
        default: ncells = 0; break;
      }

      for (c=0; c < ncells; c++)
      {
        // loop over nodes
        for (i=0; i < 8; i++)
        {
          switch (i)
          {
            case 0:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = tet_n[ltet[n]->list[c]][0];
                  n2 = tet_n[ltet[n]->list[c]][1];
                  n3 = tet_n[ltet[n]->list[c]][2];
                  n4 = tet_n[ltet[n]->list[c]][3];
                  tflag = type;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][0];
                  n2 = pyr_n[lpyr[n]->list[c]][1];
                  n3 = pyr_n[lpyr[n]->list[c]][3];
                  n4 = pyr_n[lpyr[n]->list[c]][4];
                  tflag = type;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][0];
                  n2 = pri_n[lpri[n]->list[c]][1];
                  n3 = pri_n[lpri[n]->list[c]][2];
                  n4 = pri_n[lpri[n]->list[c]][3];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][0];
                  n2 = hex_n[lhex[n]->list[c]][1];
                  n3 = hex_n[lhex[n]->list[c]][3];
                  n4 = hex_n[lhex[n]->list[c]][4];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 1:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = tet_n[ltet[n]->list[c]][1];
                  n2 = tet_n[ltet[n]->list[c]][2];
                  n3 = tet_n[ltet[n]->list[c]][0];
                  n4 = tet_n[ltet[n]->list[c]][3];
                  tflag = type;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][1];
                  n2 = pyr_n[lpyr[n]->list[c]][2];
                  n3 = pyr_n[lpyr[n]->list[c]][0];
                  n4 = pyr_n[lpyr[n]->list[c]][4];
                  tflag = type;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][1];
                  n2 = pri_n[lpri[n]->list[c]][2];
                  n3 = pri_n[lpri[n]->list[c]][0];
                  n4 = pri_n[lpri[n]->list[c]][4];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][1];
                  n2 = hex_n[lhex[n]->list[c]][2];
                  n3 = hex_n[lhex[n]->list[c]][0];
                  n4 = hex_n[lhex[n]->list[c]][5];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 2:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = tet_n[ltet[n]->list[c]][2];
                  n2 = tet_n[ltet[n]->list[c]][0];
                  n3 = tet_n[ltet[n]->list[c]][1];
                  n4 = tet_n[ltet[n]->list[c]][3];
                  tflag = type;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][2];
                  n2 = pyr_n[lpyr[n]->list[c]][3];
                  n3 = pyr_n[lpyr[n]->list[c]][1];
                  n4 = pyr_n[lpyr[n]->list[c]][4];
                  tflag = type;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][2];
                  n2 = pri_n[lpri[n]->list[c]][0];
                  n3 = pri_n[lpri[n]->list[c]][1];
                  n4 = pri_n[lpri[n]->list[c]][5];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][2];
                  n2 = hex_n[lhex[n]->list[c]][3];
                  n3 = hex_n[lhex[n]->list[c]][1];
                  n4 = hex_n[lhex[n]->list[c]][6];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 3:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = tet_n[ltet[n]->list[c]][3];
                  n2 = tet_n[ltet[n]->list[c]][0];
                  n3 = tet_n[ltet[n]->list[c]][2];
                  n4 = tet_n[ltet[n]->list[c]][1];
                  tflag = type;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][3];
                  n2 = pyr_n[lpyr[n]->list[c]][0];
                  n3 = pyr_n[lpyr[n]->list[c]][2];
                  n4 = pyr_n[lpyr[n]->list[c]][4];
                  tflag = type;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][3];
                  n2 = pri_n[lpri[n]->list[c]][5];
                  n3 = pri_n[lpri[n]->list[c]][4];
                  n4 = pri_n[lpri[n]->list[c]][0];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][3];
                  n2 = hex_n[lhex[n]->list[c]][0];
                  n3 = hex_n[lhex[n]->list[c]][2];
                  n4 = hex_n[lhex[n]->list[c]][7];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 4:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = -1;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][4];
                  n2 = pyr_n[lpyr[n]->list[c]][0];
                  n3 = pyr_n[lpyr[n]->list[c]][2];
                  n4 = pyr_n[lpyr[n]->list[c]][1];
                  tflag = 0;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][4];
                  n2 = pri_n[lpri[n]->list[c]][3];
                  n3 = pri_n[lpri[n]->list[c]][5];
                  n4 = pri_n[lpri[n]->list[c]][1];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][4];
                  n2 = hex_n[lhex[n]->list[c]][7];
                  n3 = hex_n[lhex[n]->list[c]][5];
                  n4 = hex_n[lhex[n]->list[c]][0];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 5: 
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = -1;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][4];
                  n2 = pyr_n[lpyr[n]->list[c]][1];
                  n3 = pyr_n[lpyr[n]->list[c]][3];
                  n4 = pyr_n[lpyr[n]->list[c]][2];
                  tflag = 0;
                  break;
                case 2: // Prism
                  n1 = pri_n[lpri[n]->list[c]][5];
                  n2 = pri_n[lpri[n]->list[c]][4];
                  n3 = pri_n[lpri[n]->list[c]][3];
                  n4 = pri_n[lpri[n]->list[c]][2];
                  tflag = type;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][5];
                  n2 = hex_n[lhex[n]->list[c]][4];
                  n3 = hex_n[lhex[n]->list[c]][6];
                  n4 = hex_n[lhex[n]->list[c]][1];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 6:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = -1;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][4];
                  n2 = pyr_n[lpyr[n]->list[c]][2];
                  n3 = pyr_n[lpyr[n]->list[c]][0];
                  n4 = pyr_n[lpyr[n]->list[c]][3];
                  tflag = 0;
                  break;
                case 2: // Prism
                  n1 = -1;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][6];
                  n2 = hex_n[lhex[n]->list[c]][5];
                  n3 = hex_n[lhex[n]->list[c]][7];
                  n4 = hex_n[lhex[n]->list[c]][2];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            case 7:
              switch (type)
              {
                case 0: // Tetrahedra
                  n1 = -1;
                  break;
                case 1: // Pyramid
                  n1 = pyr_n[lpyr[n]->list[c]][4];
                  n2 = pyr_n[lpyr[n]->list[c]][3];
                  n3 = pyr_n[lpyr[n]->list[c]][1];
                  n4 = pyr_n[lpyr[n]->list[c]][0];
                  tflag = 0;
                  break;
                case 2: // Prism
                  n1 = -1;
                  break;
                case 3: // Hexahdera
                  n1 = hex_n[lhex[n]->list[c]][7];
                  n2 = hex_n[lhex[n]->list[c]][6];
                  n3 = hex_n[lhex[n]->list[c]][4];
                  n4 = hex_n[lhex[n]->list[c]][3];
                  tflag = type;
                  break;
                default:
                  n1=n2=n3=n4= -1;
                  break;
              }
              break;
            default:
              n1=n2=n3=n4= -1;
              break;
          }
          if (n1 < 0 || n1 != n)
            continue;

          p1 = node[n1];
          p2 = node[n2];
          p3 = node[n3];
          p4 = node[n4];

          vec1 = Vector(p2,p3);
          vec2 = Vector(p2,p4);
          a1 = (vec1 % vec2)*-0.5;
          vec1 = Vector(p1,p4);
          vec2 = Vector(p1,p3);
          a2 = (vec1 % vec2)*-0.5;
          vec1 = Vector(p1,p2);
          vec2 = Vector(p1,p4);
          a3 = (vec1 % vec2)*-0.5;
          vec1 = Vector(p1,p3);
          vec2 = Vector(p1,p2);
          a4 = (vec1 % vec2)*-0.5;
          vec1 = Vector(p1,p2);
          vec2 = Vector(p1,p3);
          vec3 = Vector(p1,p4);
          vc = scalar_triple_product(vec1,vec2,vec3)/6.0;
          if (vc < tol)
            continue;

          rvc = 1.0/vc/3.0;

          if (E_mode > 3)
            E = (Enode[n1]+Enode[n2]+Enode[n3]+Enode[n4])*0.25;
          else
            E = Youngs_Modulus(p1,p2,p3,p4,type);
          nu = Poissons_Ratio(p1,p2,p3,p4);
          //if (pflag)
          //{
          //  E = Youngs_Modulus(p1,p2,p3,p4,type);
          //  nu = Poissons_Ratio(p1,p2,p3,p4);
          //} else
          //{
          //  q1 = p1 + Point(u[n1],v[n1],w[n1]);
          //  q2 = p2 + Point(u[n2],v[n2],w[n2]);
          //  q3 = p3 + Point(u[n3],v[n3],w[n3]);
          //  q4 = p4 + Point(u[n4],v[n4],w[n4]);
          //  E = Youngs_Modulus(q1,q2,q3,q4,type);
          //  nu = Poissons_Ratio(q1,q2,q3,q4);
          //}
          nu = MAX(0.0,MIN(0.5-1.0e-06,nu));
          if (E < Emin)
          {
            Emin = E;
            pmin = (p1+p2+p3+p4)*0.25;
          }
          if (E > Emax)
          {
            Emax = E;
            pmax = (p1+p2+p3+p4)*0.25;
          }
          if (nu < numin)
          {
            numin = nu;
            qmin = (p1+p2+p3+p4)*0.25;
          }
          //numin = MIN(numin,nu);
          if (nu > numax)
          {
            numax = nu;
            qmax = (p1+p2+p3+p4)*0.25;
          }
          //numax = MAX(numax,nu);

          a11 = a22 = a33 = E*(1.0-nu)/(1.0+nu)/(1.0-2*nu);
          a12 = a13 = a21 = a23 = a31 = a32 = 0.5*E/(1.0+nu);
          t11 = t22 = t33 = E*nu/(1.0+nu)/(1.0-2*nu);
          t12 = t13 = t21 = t23 = t31 = t32 = 0.5*E/(1.0+nu);

          vec1= Vector(p2,p3);
          vec2= Vector(p2,p4);
          t = (vec1 % vec2)*0.5;

          t1u1=(a11*a1[0]*t[0]+ a12*a1[1]*t[1]+ a13*a1[2]*t[2])*rvc;
          t1u2=(a11*a2[0]*t[0]+ a12*a2[1]*t[1]+ a13*a2[2]*t[2])*rvc;
          t1u3=(a11*a3[0]*t[0]+ a12*a3[1]*t[1]+ a13*a3[2]*t[2])*rvc;
          t1u4=(a11*a4[0]*t[0]+ a12*a4[1]*t[1]+ a13*a4[2]*t[2])*rvc;
          t1v1=(t11*a1[1]*t[0]+t12*a1[0]*t[1])*rvc;
          t1v2=(t11*a2[1]*t[0]+t12*a2[0]*t[1])*rvc;
          t1v3=(t11*a3[1]*t[0]+t12*a3[0]*t[1])*rvc;
          t1v4=(t11*a4[1]*t[0]+t12*a4[0]*t[1])*rvc;
          t1w1=(t11*a1[2]*t[0]+t13*a1[0]*t[2])*rvc;
          t1w2=(t11*a2[2]*t[0]+t13*a2[0]*t[2])*rvc;
          t1w3=(t11*a3[2]*t[0]+t13*a3[0]*t[2])*rvc;
          t1w4=(t11*a4[2]*t[0]+t13*a4[0]*t[2])*rvc;
          t1u1-= 0.5*kappa*(a1[1]*t[1]+a1[2]*t[2])*rvc;
          t1u2-= 0.5*kappa*(a2[1]*t[1]+a2[2]*t[2])*rvc;
          t1u3-= 0.5*kappa*(a3[1]*t[1]+a3[2]*t[2])*rvc;
          t1u4-= 0.5*kappa*(a4[1]*t[1]+a4[2]*t[2])*rvc;
          t1v1+= 0.5*kappa*a1[0]*t[1]*rvc;
          t1v2+= 0.5*kappa*a2[0]*t[1]*rvc;
          t1v3+= 0.5*kappa*a3[0]*t[1]*rvc;
          t1v4+= 0.5*kappa*a4[0]*t[1]*rvc;
          t1w1+= 0.5*kappa*a1[0]*t[2]*rvc;
          t1w2+= 0.5*kappa*a2[0]*t[2]*rvc;
          t1w3+= 0.5*kappa*a3[0]*t[2]*rvc;
          t1w4+= 0.5*kappa*a4[0]*t[2]*rvc;
          //t1v1+= kappa*a1[0]*t[1]*rvc;
          //t1v2+= kappa*a2[0]*t[1]*rvc;
          //t1v3+= kappa*a3[0]*t[1]*rvc;
          //t1v4+= kappa*a4[0]*t[1]*rvc;
          //t1w1+= kappa*a1[0]*t[2]*rvc;
          //t1w2+= kappa*a2[0]*t[2]*rvc;
          //t1w3+= kappa*a3[0]*t[2]*rvc;
          //t1w4+= kappa*a4[0]*t[2]*rvc;

          t2u1=(t21*a1[1]*t[0]+t22*a1[0]*t[1])*rvc;
          t2u2=(t21*a2[1]*t[0]+t22*a2[0]*t[1])*rvc;
          t2u3=(t21*a3[1]*t[0]+t22*a3[0]*t[1])*rvc;
          t2u4=(t21*a4[1]*t[0]+t22*a4[0]*t[1])*rvc;
          t2v1=(a21*a1[0]*t[0]+ a22*a1[1]*t[1]+ a23*a1[2]*t[2])*rvc;
          t2v2=(a21*a2[0]*t[0]+ a22*a2[1]*t[1]+ a23*a2[2]*t[2])*rvc;
          t2v3=(a21*a3[0]*t[0]+ a22*a3[1]*t[1]+ a23*a3[2]*t[2])*rvc;
          t2v4=(a21*a4[0]*t[0]+ a22*a4[1]*t[1]+ a23*a4[2]*t[2])*rvc;
          t2w1=(t22*a1[2]*t[1]+t23*a1[1]*t[2])*rvc;
          t2w2=(t22*a2[2]*t[1]+t23*a2[1]*t[2])*rvc;
          t2w3=(t22*a3[2]*t[1]+t23*a3[1]*t[2])*rvc;
          t2w4=(t22*a4[2]*t[1]+t23*a4[1]*t[2])*rvc;
          t2u1+= 0.5*kappa*a1[1]*t[0]*rvc;
          t2u2+= 0.5*kappa*a2[1]*t[0]*rvc;
          t2u3+= 0.5*kappa*a3[1]*t[0]*rvc;
          t2u4+= 0.5*kappa*a4[1]*t[0]*rvc;
          t2v1-= 0.5*kappa*(a1[0]*t[0]+a1[2]*t[2])*rvc;
          t2v2-= 0.5*kappa*(a2[0]*t[0]+a2[2]*t[2])*rvc;
          t2v3-= 0.5*kappa*(a3[0]*t[0]+a3[2]*t[2])*rvc;
          t2v4-= 0.5*kappa*(a4[0]*t[0]+a4[2]*t[2])*rvc;
          t2w1+= 0.5*kappa*a1[1]*t[2]*rvc;
          t2w2+= 0.5*kappa*a2[1]*t[2]*rvc;
          t2w3+= 0.5*kappa*a3[1]*t[2]*rvc;
          t2w4+= 0.5*kappa*a4[1]*t[2]*rvc;
          //t2u1+= kappa*a1[1]*t[0]*rvc;
          //t2u2+= kappa*a2[1]*t[0]*rvc;
          //t2u3+= kappa*a3[1]*t[0]*rvc;
          //t2u4+= kappa*a4[1]*t[0]*rvc;
          //t2w1+= kappa*a1[1]*t[2]*rvc;
          //t2w2+= kappa*a2[1]*t[2]*rvc;
          //t2w3+= kappa*a3[1]*t[2]*rvc;
          //t2w4+= kappa*a4[1]*t[2]*rvc;

          t3u1=(t31*a1[2]*t[0]+t33*a1[0]*t[2])*rvc;
          t3u2=(t31*a2[2]*t[0]+t33*a2[0]*t[2])*rvc;
          t3u3=(t31*a3[2]*t[0]+t33*a3[0]*t[2])*rvc;
          t3u4=(t31*a4[2]*t[0]+t33*a4[0]*t[2])*rvc;
          t3v1=(t32*a1[2]*t[1]+t33*a1[1]*t[2])*rvc;
          t3v2=(t32*a2[2]*t[1]+t33*a2[1]*t[2])*rvc;
          t3v3=(t32*a3[2]*t[1]+t33*a3[1]*t[2])*rvc;
          t3v4=(t32*a4[2]*t[1]+t33*a4[1]*t[2])*rvc;
          t3w1=(a31*a1[0]*t[0]+ a32*a1[1]*t[1]+ a33*a1[2]*t[2])*rvc;
          t3w2=(a31*a2[0]*t[0]+ a32*a2[1]*t[1]+ a33*a2[2]*t[2])*rvc;
          t3w3=(a31*a3[0]*t[0]+ a32*a3[1]*t[1]+ a33*a3[2]*t[2])*rvc;
          t3w4=(a31*a4[0]*t[0]+ a32*a4[1]*t[1]+ a33*a4[2]*t[2])*rvc;
          t3u1+= 0.5*kappa*a1[2]*t[0]*rvc;
          t3u2+= 0.5*kappa*a2[2]*t[0]*rvc;
          t3u3+= 0.5*kappa*a3[2]*t[0]*rvc;
          t3u4+= 0.5*kappa*a4[2]*t[0]*rvc;
          t3v1+= 0.5*kappa*a1[2]*t[1]*rvc;
          t3v2+= 0.5*kappa*a2[2]*t[1]*rvc;
          t3v3+= 0.5*kappa*a3[2]*t[1]*rvc;
          t3v4+= 0.5*kappa*a4[2]*t[1]*rvc;
          t3w1-= 0.5*kappa*(a1[0]*t[0]+a1[1]*t[1])*rvc;
          t3w2-= 0.5*kappa*(a2[0]*t[0]+a2[1]*t[1])*rvc;
          t3w3-= 0.5*kappa*(a3[0]*t[0]+a3[1]*t[1])*rvc;
          t3w4-= 0.5*kappa*(a4[0]*t[0]+a4[1]*t[1])*rvc;
          //t3u1+= kappa*a1[2]*t[0]*rvc;
          //t3u2+= kappa*a2[2]*t[0]*rvc;
          //t3u3+= kappa*a3[2]*t[0]*rvc;
          //t3u4+= kappa*a4[2]*t[0]*rvc;
          //t3v1+= kappa*a1[2]*t[1]*rvc;
          //t3v2+= kappa*a2[2]*t[1]*rvc;
          //t3v3+= kappa*a3[2]*t[1]*rvc;
          //t3v4+= kappa*a4[2]*t[1]*rvc;

          mat[iau[n1]][0][0] += t1u1;
          mat[iau[n1]][0][1] += t1v1;
          mat[iau[n1]][0][2] += t1w1;
          mat[iau[n1]][1][0] += t2u1;
          mat[iau[n1]][1][1] += t2v1;
          mat[iau[n1]][1][2] += t2w1;
          mat[iau[n1]][2][0] += t3u1;
          mat[iau[n1]][2][1] += t3v1;
          mat[iau[n1]][2][2] += t3w1;
          for (j=ia[n1]; j < ia[n1+1]; j++)
          {
            if (ja[j] == n2)
            {
              mat[j][0][0] += t1u2;
              mat[j][0][1] += t1v2;
              mat[j][0][2] += t1w2;
              mat[j][1][0] += t2u2;
              mat[j][1][1] += t2v2;
              mat[j][1][2] += t2w2;
              mat[j][2][0] += t3u2;
              mat[j][2][1] += t3v2;
              mat[j][2][2] += t3w2;
              break;
            }
          }
          for (j=ia[n1]; j < ia[n1+1]; j++)
          {
            if (ja[j] == n3)
            {
              mat[j][0][0] += t1u3;
              mat[j][0][1] += t1v3;
              mat[j][0][2] += t1w3;
              mat[j][1][0] += t2u3;
              mat[j][1][1] += t2v3;
              mat[j][1][2] += t2w3;
              mat[j][2][0] += t3u3;
              mat[j][2][1] += t3v3;
              mat[j][2][2] += t3w3;
              break;
            }
          }
          for (j=ia[n1]; j < ia[n1+1]; j++)
          {
            if (ja[j] == n4)
            {
              mat[j][0][0] += t1u4;
              mat[j][0][1] += t1v4;
              mat[j][0][2] += t1w4;
              mat[j][1][0] += t2u4;
              mat[j][1][1] += t2v4;
              mat[j][1][2] += t2w4;
              mat[j][2][0] += t3u4;
              mat[j][2][1] += t3v4;
              mat[j][2][2] += t3w4;
              break;
            }
          }
        }
      }
    }
  }

  if (pflag)
  {
#ifdef PARALLEL
    MPI_Allreduce(&Emin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = Emin;
    cmn_g_in.proc = my_rank;
    Emin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmin[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmin[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmin[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmin[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmin[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmin[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }

    MPI_Allreduce(&Emax,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = Emax;
    cmx_g_in.proc = my_rank;
    Emax = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmax[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmax[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmax[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmax[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmax[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmax[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }

    MPI_Allreduce(&numin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = numin;
    cmn_g_in.proc = my_rank;
    numin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&qmin[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmin[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmin[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&qmin[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmin[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmin[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }

    MPI_Allreduce(&numax,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = numax;
    cmx_g_in.proc = my_rank;
    numax = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&qmax[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmax[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmax[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&qmax[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmax[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmax[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum Young's Modulus = %g, @ (%lf, %lf, %lf)",Emin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nMaximum Young's Modulus = %g, @ (%lf, %lf, %lf)",Emax,pmax[0],pmax[1],pmax[2]);
      fprintf(out_f,"\nMinimum Poisson's Ratio = %g, @ (%lf, %lf, %lf)",numin,qmin[0],qmin[1],qmin[2]);
      fprintf(out_f,"\nMaximum Poisson's Ratio = %g, @ (%lf, %lf, %lf)",numax,qmax[0],qmax[1],qmax[2]);
      fflush(out_f);
    }
  }

  return;
}

void Smesh_obj::gg_gradients(double f[], double &fx, double &fy, double &fz,
                         List *ltet, List *lpyr, List *lpri, List *lhex)
{
  int c, i, n0, n1, n2, n3, n4, n5, n6, n7;
  Point p0, p1, p2, p3, p4, p5, p6, p7;
  double sum, vol;
  Vector cgrad;

  // Uses Green-Gauss to compute gradient at nodes

  sum = fx = fy = fz = 0.0;

  for (i=0; i < ltet->max; i++)
  {
    c = ltet->list[i];
    n0 = tet_n[c][0];
    n1 = tet_n[c][1];
    n2 = tet_n[c][2];
    n3 = tet_n[c][3];
    p0 = nodep[n0];
    p1 = nodep[n1];
    p2 = nodep[n2];
    p3 = nodep[n3];

    cgrad = tetrahedral_gradient(p0,p1,p2,p3,f[n0],f[n1],f[n2],f[n3]);
    vol = tetrahedral_volume(p0,p1,p2,p3);

    cgrad *= vol;
    fx += cgrad[0];
    fy += cgrad[1];
    fz += cgrad[2];
    sum += vol;
  }

  for (i=0; i < lpyr->max; i++)
  {
    c = lpyr->list[i];
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];
    p0 = nodep[n0];
    p1 = nodep[n1];
    p2 = nodep[n2];
    p3 = nodep[n3];
    p4 = nodep[n4];

    cgrad = pyramid_gradient(p0,p1,p2,p3,p4,f[n0],f[n1],f[n2],f[n3],f[n4]);
    vol = pyramid_volume(p0,p1,p2,p3,p4);

    cgrad *= vol;
    fx += cgrad[0];
    fy += cgrad[1];
    fz += cgrad[2];
    sum += vol;
  }

  for (i=0; i < lpri->max; i++)
  {
    c = lpri->list[i];
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];
    p0 = nodep[n0];
    p1 = nodep[n1];
    p2 = nodep[n2];
    p3 = nodep[n3];
    p4 = nodep[n4];
    p5 = nodep[n5];

    cgrad = prism_gradient(p0,p1,p2,p3,p4,p5,f[n0],f[n1],f[n2],f[n3],f[n4],f[n5]);
    vol = prism_volume(p0,p1,p2,p3,p4,p5);

    cgrad *= vol;
    fx += cgrad[0];
    fy += cgrad[1];
    fz += cgrad[2];
    sum += vol;
  }

  for (i=0; i < lhex->max; i++)
  {
    c = lhex->list[i];
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];
    p0 = nodep[n0];
    p1 = nodep[n1];
    p2 = nodep[n2];
    p3 = nodep[n3];
    p4 = nodep[n4];
    p5 = nodep[n5];
    p6 = nodep[n6];
    p7 = nodep[n7];

    cgrad = hexahedral_gradient(p0,p1,p2,p3,p4,p5,p6,p7,f[n0],f[n1],f[n2],f[n3],f[n4],f[n5],f[n6],f[n7]);
    vol = hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7);

    cgrad *= vol;
    fx += cgrad[0];
    fy += cgrad[1];
    fz += cgrad[2];
    sum += vol;
  }

  fx /= sum;
  fy /= sum;
  fz /= sum;
}
