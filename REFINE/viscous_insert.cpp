#ifdef PARALLEL
#include "mpi.h"
#endif
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
extern int stack_growth;

double stack_height(int nl, double g1, double g2, double g_min, double g_max, double &tfac)
{
  int i;
  double geo, dtotal;

  tfac = dtotal = 1.0;
  geo = g1;
  for (i=1; i < nl; i++)
  {
    tfac *= geo;
    dtotal += tfac;
    geo *= g2;
    geo = MAX(g_min,MIN(geo,g_max));
  }

  return(dtotal);
}

int viscous_insert(int v_layers, geometry *geom, char sname[],
                   int osmoo1, int osmoo2, int lsmoo, int nsmoo, int nsearch,
                   double geo1, double geo2, double geo_min, double geo_max, double aao_fac, int cflag, int restart)
{
  int b, c, e, i, j, k, l, lay, m, n, n0, n1, n2, n3, n4, n5, n6, n7, s, p, q, z, r, flag;
  int b0, b1, b2, b3, ne, nfn, gint, source, dest, local, ptag;
  int orig_nn, old_nn, nn, ntet, npyr, npri, nhex, maxtag, tgn, m0, m1, m2, m3;
  int *tag, *bft, *itmp;
  int fflag = 0;
  char fname[80];
  double mag, maxds, gdouble;
  double fac, tfac, geo;
  double tol=1.0e-12;
  double relax = 0.25;
  Point cg, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, pi, pj, pa, pb, pt;
  Vector norm, vec, v0, v1, v2, v3;
  SW_edge *edge;
  BN_FLIST *bn_face;
  Bnormal *bf_node;
  const int bdim = 132; //buffer dim
  char buff[bdim]; //buffer
#ifdef PARALLEL
  //P_VLI: declare MPI vars
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim2, *bdim21;
  char **sbuff, **rbuff;
  MPI_Status status;
#endif
  //P_VLI: gathering vars
  int *glmax = 0;
  int **new_elem = 0;
  int *nelem_owned = 0;
  int *own = 0;
  int temp_node, temp_node2, test;
  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 

  Smesh_obj *mesh1 = new Smesh_obj();
  
  sprintf(fname,"failed.cgns");
  if (num_procs > 1)
    sprintf(fname,"failed_%d.cgns",my_rank);

  if (v_layers >= 0) aao_fac = 1.0;

  //have to make file name explicit for given processor, since input deck only takes generic file name
  if (num_procs > 1)
  {
    sprintf(buff,"%s",sname);
    char *ptr = strstr(buff,".cgns");
    if (ptr == NULL)
    {
      fprintf(stderr,"\nCGNS suffix <.cgns> not found in file name!");
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    } else
    {
      //reset cursor to overwrite extension with new extension and proc num
      *ptr = '\0';
    }
    sprintf(sname,"%s_%d.cgns",buff,my_rank);
  }

  mesh1->smooth_io(-1,geom,sname);
  
  maxtag = 0;
  for (b=0; b < mesh1->nb; b++)
    maxtag = MAX(maxtag,geom->layers[b]);

  bft = new int[mesh1->nb];

  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of non-zero boundary tags = %d",maxtag);
    fflush(out_f);
  }

  // store max node for tet elements
  orig_nn = 0;
  for (c=0; c < mesh1->ntet; c++)
    for (i=0; i < 4; i++)
      orig_nn=MAX(orig_nn,mesh1->tet_n[c][i]);

  int lay_inc = 1;
  if (v_layers < 0)
    lay_inc = abs(v_layers);

  // P_VLI:  In this routine, insertion can proceed normally until the end of
  // each lock-step layer when "ghost elements" which are created on both the owning and
  // non-owning procs must be mapped to one global element.  Likewise, new nodes must get
  // the same treatment.  They will be created identically (and this will be checked for)
  // and needed in the parallel files in case they are never re-comped, but they must have
  // consistent numbering.  The beauty here is all intersections are checked by either
  // the geometry (which all procs have full access to) or neighboring elements, which are
  // also held by all interested parties.

  for (lay = abs(v_layers); lay >= 1 && fflag==0; lay -= lay_inc)
  {

    int mxlay = 0;
    int mxnew = 0;
    for (tgn=1; tgn <= maxtag; tgn++)
    {

      if (my_rank == 0)
      {
        if (v_layers > 0)
          fprintf(out_f,"\n\nInserting layer %d for boundary tag %d\n",lay,tgn);
        else
          fprintf(out_f,"\n\nInserting layers for boundary tag %d\n",tgn);
        fflush(out_f);
      }

#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      //
      // mark viscous boundary nodes 0, other boundary nodes -1, interior -2
      //
	  
      tag = (int*)malloc(mesh1->nn*sizeof(int));
      for (n=0; n < mesh1->nn; n++)
        tag[n] = -2;
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] == tgn)
          continue;
        for (i=0; i < mesh1->nt[b]; i++)
        {
          n0 = mesh1->t_n[b][i][0];
          n1 = mesh1->t_n[b][i][1];
          n2 = mesh1->t_n[b][i][2];
          tag[n0] = -1;
          tag[n1] = -1;
          tag[n2] = -1;
        }
        for (i=0; i < mesh1->nq[b]; i++)
        {
          n0 = mesh1->q_n[b][i][0];
          n1 = mesh1->q_n[b][i][1];
          n2 = mesh1->q_n[b][i][2];
          n3 = mesh1->q_n[b][i][3];
          tag[n0] = -1;
          tag[n1] = -1;
          tag[n2] = -1;
          tag[n3] = -1;
        }
      }
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] != tgn)
          continue;
        for (i=0; i < mesh1->nt[b]; i++)
        {
          n0 = mesh1->t_n[b][i][0];
          n1 = mesh1->t_n[b][i][1];
          n2 = mesh1->t_n[b][i][2];
          tag[n0] = 0;
          tag[n1] = 0;
          tag[n2] = 0;
        }
        for (i=0; i < mesh1->nq[b]; i++)
        {
          n0 = mesh1->q_n[b][i][0];
          n1 = mesh1->q_n[b][i][1];
          n2 = mesh1->q_n[b][i][2];
          n3 = mesh1->q_n[b][i][3];
          tag[n0] = 0;
          tag[n1] = 0;
          tag[n2] = 0;
          tag[n3] = 0;
        }
      }
      // identify adjacent floating boundaries
      for (b=0; b < mesh1->nb; b++)
      {
        bft[b] = 0;
        if (geom->layers[b] == tgn)
          continue;
        j = 0;
        for (i=0; i < mesh1->nt[b] && !j; i++)
        {
          n0 = mesh1->t_n[b][i][0];
          n1 = mesh1->t_n[b][i][1];
          n2 = mesh1->t_n[b][i][2];
          if (tag[n0]==0 || tag[n1]==0 || tag[n2]==0) j=1;
        }
        for (i=0; i < mesh1->nq[b] && !j; i++)
        {
          n0 = mesh1->q_n[b][i][0];
          n1 = mesh1->q_n[b][i][1];
          n2 = mesh1->q_n[b][i][2];
          n3 = mesh1->q_n[b][i][3];
          if (tag[n0]==0 || tag[n1]==0 || tag[n2]==0 || tag[n3]==0) j=1;
        }
        bft[b]=j;
      }

      // create floating boundary node structure array
      itmp = new int[mesh1->nn];
      nfn = 0;
      for (b=0; b < mesh1->nb; b++)
      {
        if (!bft[b])
          continue;
        for (n=0; n < mesh1->nn; n++)
          itmp[n] = 0;
        for (i=0; i < mesh1->nt[b]; i++)
        {
          for (j=0; j < 3; j++)
          {
            n = mesh1->t_n[b][i][j];
            if (tag[n]==0) itmp[n]=1;
          }
        }
        for (i=0; i < mesh1->nq[b]; i++)
        {
          for (j=0; j < 4; j++)
          {
            n = mesh1->q_n[b][i][j];
            if (tag[n]==0) itmp[n]=1;
          }
        }
        for (n=0; n < mesh1->nn; n++)
          nfn += itmp[n];
      }

      gint = nfn;
#ifdef PARALLEL
      MPI_Allreduce(&nfn,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0)
      {
        fprintf(out_f,"\nLayer %d, # of floating boundary nodes = %d",lay,gint);
        fflush(out_f);
      }
      if (nfn > 0)
      {
        bf_node = new Bnormal[nfn];
        nfn = 0;
        for (b=0; b < mesh1->nb; b++)
        {
          if (!bft[b])
            continue;
          for (n=0; n < mesh1->nn; n++)
            itmp[n] = 0;
          for (i=0; i < mesh1->nt[b]; i++)
          {
            for (j=0; j < 3; j++)
            {
              n = mesh1->t_n[b][i][j];
              if (tag[n]==0) itmp[n]=1;
            }
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            for (j=0; j < 4; j++)
            {
              n = mesh1->q_n[b][i][j];
              if (tag[n]==0) itmp[n]=1;
            }
          }
          for (n=0; n < mesh1->nn; n++)
          {
            if (itmp[n] == 1)
            {
              bf_node[nfn].index = n;
              bf_node[nfn].boundary = b;
              i = b;
              norm = Vector(0.0,0.0,0.0);
              p1 = geom->closest(mesh1->nodep[n],i,norm);
              bf_node[nfn].norm = norm;
              nfn++;
            }
          }
        }
      }
      delete[] itmp;

      // create list of faces per boundary node
      bn_face = new BN_FLIST[mesh1->nn];
      for (n=0; n < mesh1->nn; n++)
      {
        bn_face[n].nf = 0;
        bn_face[n].face = 0;
      }
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] != tgn)
          continue;
        for (i=0; i < mesh1->nt[b]; i++)
        {
          n0 = mesh1->t_n[b][i][0];
          n1 = mesh1->t_n[b][i][1];
          n2 = mesh1->t_n[b][i][2];
          bn_face[n0].face = (FACE*)realloc((void*)bn_face[n0].face,(bn_face[n0].nf+1)*sizeof(FACE));
          bn_face[n0].face[bn_face[n0].nf].node[0] = n0;
          bn_face[n0].face[bn_face[n0].nf].node[1] = n1;
          bn_face[n0].face[bn_face[n0].nf].node[2] = n2;
          bn_face[n0].face[bn_face[n0].nf].node[3] = -1;
          bn_face[n0].nf++;
          bn_face[n1].face = (FACE*)realloc((void*)bn_face[n1].face,(bn_face[n1].nf+1)*sizeof(FACE));
          bn_face[n1].face[bn_face[n1].nf].node[0] = n1;
          bn_face[n1].face[bn_face[n1].nf].node[1] = n2;
          bn_face[n1].face[bn_face[n1].nf].node[2] = n0;
          bn_face[n1].face[bn_face[n1].nf].node[3] = -1;
          bn_face[n1].nf++;
          bn_face[n2].face = (FACE*)realloc((void*)bn_face[n2].face,(bn_face[n2].nf+1)*sizeof(FACE));
          bn_face[n2].face[bn_face[n2].nf].node[0] = n2;
          bn_face[n2].face[bn_face[n2].nf].node[1] = n0;
          bn_face[n2].face[bn_face[n2].nf].node[2] = n1;
          bn_face[n2].face[bn_face[n2].nf].node[3] = -1;
          bn_face[n2].nf++;
        }
        for (i=0; i < mesh1->nq[b]; i++)
        {
          n0 = mesh1->q_n[b][i][0];
          n1 = mesh1->q_n[b][i][1];
          n2 = mesh1->q_n[b][i][2];
          n3 = mesh1->q_n[b][i][3];
          bn_face[n0].face = (FACE*)realloc((void*)bn_face[n0].face,(bn_face[n0].nf+1)*sizeof(FACE));
          bn_face[n0].face[bn_face[n0].nf].node[0] = n0;
          bn_face[n0].face[bn_face[n0].nf].node[1] = n1;
          bn_face[n0].face[bn_face[n0].nf].node[2] = n2;
          bn_face[n0].face[bn_face[n0].nf].node[3] = n3;
          bn_face[n0].nf++;
          bn_face[n1].face = (FACE*)realloc((void*)bn_face[n1].face,(bn_face[n1].nf+1)*sizeof(FACE));
          bn_face[n1].face[bn_face[n1].nf].node[0] = n1;
          bn_face[n1].face[bn_face[n1].nf].node[1] = n2;
          bn_face[n1].face[bn_face[n1].nf].node[2] = n3;
          bn_face[n1].face[bn_face[n1].nf].node[3] = n0;
          bn_face[n1].nf++;
          bn_face[n2].face = (FACE*)realloc((void*)bn_face[n2].face,(bn_face[n2].nf+1)*sizeof(FACE));
          bn_face[n2].face[bn_face[n2].nf].node[0] = n2;
          bn_face[n2].face[bn_face[n2].nf].node[1] = n3;
          bn_face[n2].face[bn_face[n2].nf].node[2] = n0;
          bn_face[n2].face[bn_face[n2].nf].node[3] = n1;
          bn_face[n2].nf++;
          bn_face[n3].face = (FACE*)realloc((void*)bn_face[n3].face,(bn_face[n3].nf+1)*sizeof(FACE));
          bn_face[n3].face[bn_face[n3].nf].node[0] = n3;
          bn_face[n3].face[bn_face[n3].nf].node[1] = n0;
          bn_face[n3].face[bn_face[n3].nf].node[2] = n1;
          bn_face[n3].face[bn_face[n3].nf].node[3] = n2;
          bn_face[n3].nf++;
        }
      }

      // for boundary nodes, compute bisector plane normal for 2 facets with minimum angle between
      double dot, dotmax;
      Vector *bnorm, vn, vt, normi, normj;
      bnorm = new Vector[mesh1->nn];
      int ii, jj;
      b=c=0;
      for (n=0; n < mesh1->nn; n++)
      {
        bnorm[n]=Vector(0.0,0.0,0.0);
        if (tag[n] != 0)
          continue;
        b++;
        dotmax=0.50;
        ii=jj= -1;
        for (i=0; i < bn_face[n].nf; i++)
        {
          m0 = bn_face[n].face[i].node[0];
          m1 = bn_face[n].face[i].node[1];
          m2 = bn_face[n].face[i].node[2];
          m3 = bn_face[n].face[i].node[3];
          p0 = mesh1->nodep[m0];
          p1 = mesh1->nodep[m1];
          p2 = mesh1->nodep[m2];
          if (m3 >= 0)
          {
            p3 = mesh1->nodep[m3];
            v1 = Vector(p0,p2);
            v2 = Vector(p1,p3);
            m = m3;
            pi = (p0+p1+p2+p3)*0.25;
          } else
          {
            v1 = Vector(p0,p1);
            v2 = Vector(p0,p2);
            m = m2;
            pi = (p0+p1+p2)/3.0;
          }
          normi = v1 % v2;
          normi.normalize();
          for (j=0; j < bn_face[n].nf; j++)
          {
            if (i==j) continue;
            n0 = bn_face[n].face[j].node[0];
            n1 = bn_face[n].face[j].node[1];
            n2 = bn_face[n].face[j].node[2];
            n3 = bn_face[n].face[j].node[3];
            if (n1 != m)
              continue;
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            if (n3 >= 0)
            {
              p3 = mesh1->nodep[n3];
              v1 = Vector(p0,p2);
              v2 = Vector(p1,p3);
              pj = (p0+p1+p2+p3)*0.25;
            } else
            {
              v1 = Vector(p0,p1);
              v2 = Vector(p0,p2);
              pj = (p0+p1+p2)/3.0;
            }
            normj = v1 % v2;
            normj.normalize();
            dot = normi * normj;
            v3 = Vector(mesh1->nodep[n],(pi+pj)*0.5);
            // save facets with minimum concave angle
            //if (dot < dotmax && (normi*v3 > 0.0 && normj*v3 > 0.0))
            // save facets with minimum concave or convex angle
            if (dot < dotmax && ((normi*v3 > 0.0 && normj*v3 > 0.0) ||
                                 (normi*v3 < 0.0 && normj*v3 < 0.0)))
            {
              dotmax = dot;
              ii = i;
              jj = j;
            }
          }
        }
        if (ii < 0 || jj < 0)
          continue;
        c++;
        m0 = bn_face[n].face[ii].node[0];
        m1 = bn_face[n].face[ii].node[1];
        m2 = bn_face[n].face[ii].node[2];
        m3 = bn_face[n].face[ii].node[3];
        p0 = mesh1->nodep[m0];
        p1 = mesh1->nodep[m1];
        p2 = mesh1->nodep[m2];
        if (m3 >= 0)
        {
          p3 = mesh1->nodep[m3];
          v1 = Vector(p0,p2);
          v2 = Vector(p1,p3);
        } else
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p0,p2);
        }
        normi = v1 % v2;
        normi.normalize();
        n0 = bn_face[n].face[jj].node[0];
        n1 = bn_face[n].face[jj].node[1];
        n2 = bn_face[n].face[jj].node[2];
        n3 = bn_face[n].face[jj].node[3];
        p0 = mesh1->nodep[n0];
        p1 = mesh1->nodep[n1];
        p2 = mesh1->nodep[n2];
        if (n3 >= 0)
        {
          p3 = mesh1->nodep[n3];
          v1 = Vector(p0,p2);
          v2 = Vector(p1,p3);
        } else
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p0,p2);
        }
        normj = v1 % v2;
        normj.normalize();
        v1 = (normi + normj)*0.5;
        v2 = normi % normj;
        bnorm[n] = v1 % v2;
        bnorm[n].normalize();
      }
      for (n=0; n < mesh1->nn; n++)
        free(bn_face[n].face);
      delete[] bn_face;

      gint = b;
#ifdef PARALLEL
      MPI_Allreduce(&b,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nLayer %d, # of marching boundary nodes = %d",lay,gint);
      gint = c;
#ifdef PARALLEL
      MPI_Allreduce(&c,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nLayer %d, # of bisector planes computed= %d",lay,gint);
      if (my_rank == 0) fflush(out_f);

      old_nn = mesh1->nn;
      // compute desired normals & spacings..P_VLI:  will compute incorrectly for ghost nodes and have them recomm'd before smoothing
      double *desired, ds, ang, ang0, ang1, ang2, ang3;
      Vector *nrm, *tng, dv;
      double *wgt;
      nrm = new Vector[mesh1->nn];
      tng = new Vector[mesh1->nn];
      desired = new double[mesh1->nn];
      wgt = new double[mesh1->nn];
      itmp = new int[mesh1->nn];
      
      //before determining normals, sync tag with owning proc
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
      n=0;
      mesh1->exchange_tag(tag,n);
#endif

      for (n=0; n < mesh1->nn; n++)
      {
        wgt[n] = 0.0;
        nrm[n] = Vector(0.0,0.0,0.0);
        tng[n] = Vector(0.0,0.0,0.0);
      }
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] == tgn)
        {
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            v1 = Vector(p0,p1);
            v1.normalize();
            v2 = Vector(p0,p2);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n0];
            norm -= bnorm[n0]*dot;
            norm.normalize();
            nrm[n0] += norm*ang;
            wgt[n0] += ang;
            v1 = Vector(p1,p2);
            v1.normalize();
            v2 = Vector(p1,p0);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n1];
            norm -= bnorm[n1]*dot;
            norm.normalize();
            nrm[n1] += norm*ang;
            wgt[n1] += ang;
            v1 = Vector(p2,p0);
            v1.normalize();
            v2 = Vector(p2,p1);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n2];
            norm -= bnorm[n2]*dot;
            norm.normalize();
            nrm[n2] += norm*ang;
            wgt[n2] += ang;
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            v1 = Vector(p0,p1);
            v1.normalize();
            v2 = Vector(p0,p3);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n0];
            norm -= bnorm[n0]*dot;
            norm.normalize();
            nrm[n0] += norm*ang;
            wgt[n0] += ang;
            v1 = Vector(p1,p2);
            v1.normalize();
            v2 = Vector(p1,p0);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n1];
            norm -= bnorm[n1]*dot;
            norm.normalize();
            nrm[n1] += norm*ang;
            wgt[n1] += ang;
            v1 = Vector(p2,p3);
            v1.normalize();
            v2 = Vector(p2,p1);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n2];
            norm -= bnorm[n2]*dot;
            norm.normalize();
            nrm[n2] += norm*ang;
            wgt[n2] += ang;
            v1 = Vector(p3,p0);
            v1.normalize();
            v2 = Vector(p3,p2);
            v2.normalize();
            ang = acos(v1*v2);
            norm = v1 % v2;
            dot = norm * bnorm[n3];
            norm -= bnorm[n3]*dot;
            norm.normalize();
            nrm[n3] += norm*ang;
            wgt[n3] += ang;
          }
        }
      }

      // Reset tangent nodes
      for (i=0; i < nfn; i++)
      {
        n = bf_node[i].index;
        nrm[n]=Vector(0.0,0.0,0.0);
        wgt[n]=0.0;
      }
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] != tgn)
        {
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            if (tag[n0]!=0 && tag[n1]!=0 && tag[n2]!=0)
              continue;
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            v1 = Vector(p0,p1);
            v2 = Vector(p0,p2);
            norm = v1 % v2;
            norm.normalize();
            if (tag[n0]==0 && tag[n1]==0)
            {
              v1 = Vector(p1,p0);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n0] += v2;
              wgt[n0] += 1.0;
              nrm[n1] += v2;
              wgt[n1] += 1.0;
            }
            if (tag[n1]==0 && tag[n2]==0)
            {
              v1 = Vector(p2,p1);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n1] += v2;
              wgt[n1] += 1.0;
              nrm[n2] += v2;
              wgt[n2] += 1.0;
            }
            if (tag[n2]==0 && tag[n0]==0)
            {
              v1 = Vector(p0,p2);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n2] += v2;
              wgt[n2] += 1.0;
              nrm[n0] += v2;
              wgt[n0] += 1.0;
            }
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            if (tag[n0]!=0 && tag[n1]!=0 && tag[n2]!=0 && tag[n3]!=0)
              continue;
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            v1 = Vector(p0,p2);
            v2 = Vector(p1,p3);
            norm = v1 % v2;
            norm.normalize();
            if (tag[n0]==0 && tag[n1]==0)
            {
              v1 = Vector(p1,p0);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n0] += v2;
              wgt[n0] += 1.0;
              nrm[n1] += v2;
              wgt[n1] += 1.0;
            }
            if (tag[n1]==0 && tag[n2]==0)
            {
              v1 = Vector(p2,p1);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n1] += v2;
              wgt[n1] += 1.0;
              nrm[n2] += v2;
              wgt[n2] += 1.0;
            }
            if (tag[n2]==0 && tag[n3]==0)
            {
              v1 = Vector(p3,p2);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n2] += v2;
              wgt[n2] += 1.0;
              nrm[n3] += v2;
              wgt[n3] += 1.0;
            }
            if (tag[n3]==0 && tag[n0]==0)
            {
              v1 = Vector(p0,p3);
              v2 = v1 % norm;
              v2.normalize();
              nrm[n3] += v2;
              wgt[n3] += 1.0;
              nrm[n0] += v2;
              wgt[n0] += 1.0;
            }
          }
        }
      }

      for (n=0; n < mesh1->nn; n++)
      {
        if (tag[n] > -2 && wgt[n] > 1.0e-20)
        {
          nrm[n] /= wgt[n];
          nrm[n].normalize();
        }
      }

      // detect sidewalls
      int sidewall = 0;
      ne = 0;
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] != tgn)
        {
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            if (tag[n0] == 0 || tag[n1] == 0 || tag[n2] == 0)
              sidewall = 1;
            if (tag[n0] == 0 && tag[n1] == 0)
              ne++;
            if (tag[n1] == 0 && tag[n2] == 0)
              ne++;
            if (tag[n2] == 0 && tag[n0] == 0)
              ne++;
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            if (tag[n0] == 0 || tag[n1] == 0 || tag[n2] == 0 || tag[n3] == 0)
              sidewall = 1;
            if (tag[n0] == 0 && tag[n1] == 0)
              ne++;
            if (tag[n1] == 0 && tag[n2] == 0)
              ne++;
            if (tag[n2] == 0 && tag[n3] == 0)
              ne++;
            if (tag[n3] == 0 && tag[n0] == 0)
              ne++;
          }
        }
      }
      if (ne > 0)
      {
        edge = new SW_edge[ne];
        ne = 0;
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] != tgn)
          {
            for (i=0; i < mesh1->nt[b]; i++)
            {
              for (s=0; s < 3; s++)
              {
                switch (s)
                {
                  case 0: n0 = mesh1->t_n[b][i][0]; n1 = mesh1->t_n[b][i][1]; break;
                  case 1: n0 = mesh1->t_n[b][i][1]; n1 = mesh1->t_n[b][i][2]; break;
                  case 2: n0 = mesh1->t_n[b][i][2]; n1 = mesh1->t_n[b][i][0]; break;
                  default: n0=n1= -1;
                }
                if (n0 >= 0 && n1 >= 0 && tag[n0] == 0 && tag[n1] == 0)
                {
                  k=0;
                  for (j=0; j < ne && !k; j++)
                  {
                    if (edge[j].nodes[1] == n0 && edge[j].nodes[0] == n1)
                    {
                      k=1;
                      edge[j].nodes[0] = edge[ne-1].nodes[0];
                      edge[j].nodes[1] = edge[ne-1].nodes[1];
                      edge[j].boundary = edge[ne-1].boundary;
                      ne--;
                    }
                  }
                  if (!k)
                  {
                    edge[ne].nodes[0] = n0;
                    edge[ne].nodes[1] = n1;
                    edge[ne].boundary = b;
                    ne++;
                  }
                }
              }
            }
            for (i=0; i < mesh1->nq[b]; i++)
            {
              for (s=0; s < 4; s++)
              {
                switch (s)
                {
                  case 0: n0 = mesh1->q_n[b][i][0]; n1 = mesh1->q_n[b][i][1]; break;
                  case 1: n0 = mesh1->q_n[b][i][1]; n1 = mesh1->q_n[b][i][2]; break;
                  case 2: n0 = mesh1->q_n[b][i][2]; n1 = mesh1->q_n[b][i][3]; break;
                  case 3: n0 = mesh1->q_n[b][i][3]; n1 = mesh1->q_n[b][i][0]; break;
                  default: n0=n1= -1;
                }
                if (n0 >= 0 && n1 >= 0 && tag[n0] == 0 && tag[n1] == 0)
                {
                  k=0;
                  for (j=0; j < ne && !k; j++)
                  {
                    if (edge[j].nodes[1] == n0 && edge[j].nodes[0] == n1)
                    {
                      k=1;
                      edge[j].nodes[0] = edge[ne-1].nodes[0];
                      edge[j].nodes[1] = edge[ne-1].nodes[1];
                      edge[j].boundary = edge[ne-1].boundary;
                      ne--;
                    }
                  }
                  if (!k)
                  {
                    edge[ne].nodes[0] = n0;
                    edge[ne].nodes[1] = n1;
                    edge[ne].boundary = b;
                    ne++;
                  }
                }
              }
            }
          }
        }
      }
      gint = ne;
#ifdef PARALLEL
      MPI_Allreduce(&ne,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nNumber of sidewall edges = %d",gint);

      // impose sidewall tangency
      for (i=0; i < nfn; i++)
      {
        n = bf_node[i].index;
        norm = bf_node[i].norm;
        dot = nrm[n] * bnorm[n];
        nrm[n] -= bnorm[n]*dot;
        dot = norm * nrm[n];
        nrm[n] -= norm*dot;
        nrm[n].normalize();
      }
      // create node-to-node hash
      List **nhash;
      nhash = new List*[mesh1->nn];
      for (n=0; n < mesh1->nn; n++)
        nhash[n] = new List();

      mesh1->create_node_to_node(nhash);

      // set critical floating node to line up with chain
      for (b=0; b < mesh1->nb; b++)
      {
        if (!bft[b])
          continue;

        for (n=0; n < mesh1->nn; n++)
          itmp[n] = 0;
        for (i=0; i < mesh1->nt[b]; i++)
        {
          n0 = mesh1->t_n[b][i][0];
          n1 = mesh1->t_n[b][i][1];
          n2 = mesh1->t_n[b][i][2];
          itmp[n0] = 1;
          itmp[n1] = 1;
          itmp[n2] = 1;
        }
        for (i=0; i < mesh1->nq[b]; i++)
        {
          n0 = mesh1->q_n[b][i][0];
          n1 = mesh1->q_n[b][i][1];
          n2 = mesh1->q_n[b][i][2];
          n3 = mesh1->q_n[b][i][3];
          itmp[n0] = 1;
          itmp[n1] = 1;
          itmp[n2] = 1;
          itmp[n3] = 1;
        }
        for (j=0; j < mesh1->nb; j++)
        {
          if (!bft[j] || j==b)
            continue;
          for (i=0; i < mesh1->nt[j]; i++)
          {
            n0 = mesh1->t_n[j][i][0];
            n1 = mesh1->t_n[j][i][1];
            n2 = mesh1->t_n[j][i][2];
            if (itmp[n0] == 1) itmp[n0] = 2;
            if (itmp[n1] == 1) itmp[n1] = 2;
            if (itmp[n2] == 1) itmp[n2] = 2;
          }
          for (i=0; i < mesh1->nq[j]; i++)
          {
            n0 = mesh1->q_n[j][i][0];
            n1 = mesh1->q_n[j][i][1];
            n2 = mesh1->q_n[j][i][2];
            n3 = mesh1->q_n[j][i][3];
            if (itmp[n0] == 1) itmp[n0] = 2;
            if (itmp[n1] == 1) itmp[n1] = 2;
            if (itmp[n2] == 1) itmp[n2] = 2;
            if (itmp[n3] == 1) itmp[n3] = 2;
          }
          for (n=0; n < mesh1->nn; n++)
          {
            if (tag[n] != 0 || itmp[n] != 2)
              continue;
            for (i=0; i < nhash[n]->max; i++)
            {
              m = nhash[n]->list[i];
              if (tag[m] < 0 && itmp[m] == 2)
              {
                nrm[n] = Vector(mesh1->nodep[n],mesh1->nodep[m]);
                nrm[n].normalize();
                break;
              }
            }
          }
        }
      }

      // optimize & smooth normals 
      int nopt = 0;
      int optmax = 100;
      int knt;
      double rms;
      int *ntag = new int[mesh1->nn];
      for (n=0; n < mesh1->nn; n++)
      {
        if (bnorm[n].magnitude() > 0.01)
          ntag[n] = 1;
        else
          ntag[n] = 0;
      }
      rms = 0.0;
      knt = 0;
      
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
      mesh1->exchange_Vector(nrm);
#endif

      List *punt = new List();

      double globalrms, negmax;
      int globalknt, globalnrms, nmax;
      globalrms = rms;
      globalknt = knt;
      for (nopt=1; nopt <= MAX(1,nsmoo) || (nopt <= optmax && globalknt > 0); nopt++)
      {
        negmax = 0.0;
        nmax = -1;
        //begin normal smoothing
            
        for (n=0; n < mesh1->nn; n++)
        {
          tng[n] = Vector(0.0,0.0,0.0);
          wgt[n] = 0.0;
        }
        j=0;
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] == tgn)
          {
            for (i=0; i < mesh1->nt[b]; i++)
            {
              n0 = mesh1->t_n[b][i][0];
              n1 = mesh1->t_n[b][i][1];
              n2 = mesh1->t_n[b][i][2];
              p0 = mesh1->nodep[n0];
              p1 = mesh1->nodep[n1];
              p2 = mesh1->nodep[n2];

              v1 = Vector(p0,p1);
              v2 = Vector(p0,p2);
              norm = v1 % v2;
              norm.normalize();

              mag = 1.0 - norm * nrm[n0];
              if (mag > 1.0)
              {
                j++;
                ntag[n0] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n0;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p0[0],p0[1],p0[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n0);
              }
              if (ntag[n0])
              {
                //v1 = Vector(p0,p1);
                //v1.normalize();
                //v2 = Vector(p0,p2);
                //v2.normalize();
                //ang0 = acos(v1*v2);
                dot = norm * bnorm[n0];
                v1 = norm - bnorm[n0]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n0];
                tng[n0] += v1*mag;
                wgt[n0] += mag;
              }

              mag = 1.0 - norm * nrm[n1];
              if (mag > 1.0)
              {
                j++;
                ntag[n1] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n1;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p1[0],p1[1],p1[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n1);
              }
              if (ntag[n1])
              {
                //v1 = Vector(p1,p2);
                //v1.normalize();
                //v2 = Vector(p1,p0);
                //v2.normalize();
                //ang1 = acos(v1*v2);
                dot = norm * bnorm[n1];
                v1 = norm - bnorm[n1]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n1];
                tng[n1] += v1*mag;
                wgt[n1] += mag;
              }

              mag = 1.0 - norm * nrm[n2];
              if (mag > 1.0)
              {
                j++;
                ntag[n2] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n2;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p2[0],p2[1],p2[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n2);
              }
              if (ntag[n2])
              {
                //v1 = Vector(p2,p0);
                //v1.normalize();
                //v2 = Vector(p2,p1);
                //v2.normalize();
                //ang2 = acos(v1*v2);
                dot = norm * bnorm[n2];
                v1 = norm - bnorm[n2]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n2];
                tng[n2] += v1*mag;
                wgt[n2] += mag;
              }

              ds = distance(p0,p1);
              mag = (1.0-ntag[n0])*(1.0-nrm[n0]*nrm[n1])/ds;
              dot = nrm[n1] * bnorm[n0];
              v1 = nrm[n1] - bnorm[n0]*dot;
              v1.normalize();
              tng[n0] += v1*mag;
              wgt[n0] += mag;
              mag = (1.0-ntag[n1])*(1.0-nrm[n0]*nrm[n1])/ds;
              dot = nrm[n0] * bnorm[n1];
              v1 = nrm[n0] - bnorm[n1]*dot;
              v1.normalize();
              tng[n1] += v1*mag;
              wgt[n1] += mag;

              ds = distance(p1,p2);
              mag = (1.0-ntag[n1])*(1.0-nrm[n1]*nrm[n2])/ds;
              dot = nrm[n2] * bnorm[n1];
              v1 = nrm[n2] - bnorm[n1]*dot;
              v1.normalize();
              tng[n1] += v1*mag;
              wgt[n1] += mag;
              mag = (1.0-ntag[n2])*(1.0-nrm[n1]*nrm[n2])/ds;
              dot = nrm[n1] * bnorm[n2];
              v1 = nrm[n1] - bnorm[n2]*dot;
              v1.normalize();
              tng[n2] += v1*mag;
              wgt[n2] += mag;

              ds = distance(p0,p2);
              mag = (1.0-ntag[n2])*(1.0-nrm[n0]*nrm[n2])/ds;
              dot = nrm[n0] * bnorm[n2];
              v1 = nrm[n0] - bnorm[n2]*dot;
              v1.normalize();
              tng[n2] += v1*mag;
              wgt[n2] += mag;
              mag = (1.0-ntag[n0])*(1.0-nrm[n0]*nrm[n2])/ds;
              dot = nrm[n2] * bnorm[n0];
              v1 = nrm[n2] - bnorm[n0]*dot;
              v1.normalize();
              tng[n0] += v1*mag;
              wgt[n0] += mag;
            }
            for (i=0; i < mesh1->nq[b]; i++)
            {
              n0 = mesh1->q_n[b][i][0];
              n1 = mesh1->q_n[b][i][1];
              n2 = mesh1->q_n[b][i][2];
              n3 = mesh1->q_n[b][i][3];
              p0 = mesh1->nodep[n0];
              p1 = mesh1->nodep[n1];
              p2 = mesh1->nodep[n2];
              p3 = mesh1->nodep[n3];

              //v1 = Vector(p0,p1);
              //v2 = Vector(p0,p3);
              v1 = Vector(p0,p2);
              v2 = Vector(p1,p3);
              norm = v1 % v2;
              norm.normalize();

              mag = 1.0 - norm * nrm[n0];
              if (mag > 1.0)
              {
                j++;
                ntag[n0] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n0;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p0[0],p0[1],p0[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n0);
              }
              if (ntag[n0])
              {
                //v1 = Vector(p0,p1);
                //v1.normalize();
                //v2 = Vector(p0,p3);
                //v2.normalize();
                //ang0 = acos(v1*v2);
                dot = norm * bnorm[n0];
                v1 = norm - bnorm[n0]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n0];
                tng[n0] += v1*mag;
                wgt[n0] += mag;
              }

              mag = 1.0 - norm * nrm[n1];
              if (mag > 1.0)
              {
                j++;
                ntag[n1] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n1;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p1[0],p1[1],p1[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n1);
              }
              if (ntag[n1])
              {
                //v1 = Vector(p1,p2);
                //v1.normalize();
                //v2 = Vector(p1,p0);
                //v2.normalize();
                //ang1 = acos(v1*v2);
                dot = norm * bnorm[n1];
                v1 = norm - bnorm[n1]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n1];
                tng[n1] += v1*mag;
                wgt[n1] += mag;
              }

              mag = 1.0 - norm * nrm[n2];
              if (mag > 1.0)
              {
                j++;
                ntag[n2] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n2;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p2[0],p2[1],p2[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n2);
              }
              if (ntag[n2])
              {
                //v1 = Vector(p2,p3);
                //v1.normalize();
                //v2 = Vector(p2,p1);
                //v2.normalize();
                //ang2 = acos(v1*v2);
                dot = norm * bnorm[n2];
                v1 = norm - bnorm[n2]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n2];
                tng[n2] += v1*mag;
                wgt[n2] += mag;
              }

              mag = 1.0 - norm * nrm[n3];
              if (mag > 1.0)
              {
                j++;
                ntag[n3] = 1;
                if (mag > negmax)
                {
                  negmax = mag;
                  nmax = n3;
                }
                if (num_procs == 1)
                {
                  fprintf(out_f,"\nNegative Dot Product occurs @ (%lg, %lg, %lg)",p3[0],p3[1],p3[2]);
                  fflush(out_f);
                }
                if (nopt == optmax)
                  punt->Check_List(n3);
              }
              if (ntag[n3])
              {
                //v1 = Vector(p3,p0);
                //v1.normalize();
                //v2 = Vector(p3,p2);
                //v2.normalize();
                //ang3 = acos(v1*v2);
                dot = norm * bnorm[n3];
                v1 = norm - bnorm[n3]*dot;
                v1.normalize();
                mag = 1.0 - v1 * nrm[n3];
                tng[n3] += v1*mag;
                wgt[n3] += mag;
              }

              ds = distance(p0,p1);
              mag = (1.0-ntag[n0])*(1.0-nrm[n0]*nrm[n1])/ds;
              dot = nrm[n1] * bnorm[n0];
              v1 = nrm[n1] - bnorm[n0]*dot;
              v1.normalize();
              tng[n0] += v1*mag;
              wgt[n0] += mag;
              mag = (1.0-ntag[n1])*(1.0-nrm[n0]*nrm[n1])/ds;
              dot = nrm[n0] * bnorm[n1];
              v1 = nrm[n0] - bnorm[n1]*dot;
              v1.normalize();
              tng[n1] += v1*mag;
              wgt[n1] += mag;

              ds = distance(p1,p2);
              mag = (1.0-ntag[n1])*(1.0-nrm[n1]*nrm[n2])/ds;
              dot = nrm[n2] * bnorm[n1];
              v1 = nrm[n2] - bnorm[n1]*dot;
              v1.normalize();
              tng[n1] += v1*mag;
              wgt[n1] += mag;
              mag = (1.0-ntag[n2])*(1.0-nrm[n1]*nrm[n2])/ds;
              dot = nrm[n1] * bnorm[n2];
              v1 = nrm[n1] - bnorm[n2]*dot;
              v1.normalize();
              tng[n2] += v1*mag;
              wgt[n2] += mag;

              ds = distance(p2,p3);
              mag = (1.0-ntag[n2])*(1.0-nrm[n2]*nrm[n3])/ds;
              dot = nrm[n3] * bnorm[n2];
              v1 = nrm[n3] - bnorm[n2]*dot;
              v1.normalize();
              tng[n2] += v1*mag;
              wgt[n2] += mag;
              mag = (1.0-ntag[n3])*(1.0-nrm[n2]*nrm[n3])/ds;
              dot = nrm[n2] * bnorm[n3];
              v1 = nrm[n2] - bnorm[n3]*dot;
              v1.normalize();
              tng[n3] += v1*mag;
              wgt[n3] += mag;

              ds = distance(p3,p0);
              mag = (1.0-ntag[n3])*(1.0-nrm[n0]*nrm[n3])/ds;
              dot = nrm[n0] * bnorm[n3];
              v1 = nrm[n0] - bnorm[n3]*dot;
              v1.normalize();
              tng[n3] += v1*mag;
              wgt[n3] += mag;
              mag = (1.0-ntag[n0])*(1.0-nrm[n0]*nrm[n3])/ds;
              dot = nrm[n3] * bnorm[n0];
              v1 = nrm[n3] - bnorm[n0]*dot;
              v1.normalize();
              tng[n0] += v1*mag;
              wgt[n0] += mag;
            }
          }
        }
        gint = j;
#ifdef PARALLEL
        MPI_Allreduce(&j,&gint,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
        if (gint > 0)
        {
          knt = 3;
          cg = Point(0.0,0.0,0.0);
          if (nmax >= 0)
            cg = mesh1->nodep[nmax];
#ifdef PARALLEL
          cmx_g_out.val = cmx_g_in.val = negmax;
          cmx_g_out.proc = cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&(cg[0]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[1]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[2]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&cg[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif
          if (my_rank == 0)
          {
            fprintf(out_f,"\nMAX Negative Dot Product occurs @ (%lg, %lg, %lg)",cg[0],cg[1],cg[2]);
            fflush(out_f);
          }
        } else
          knt = MAX(1,knt)-1;

        // Reset tangent nodes
        for (i=0; i < nfn; i++)
        {
          n = bf_node[i].index;
          tng[n]=Vector(0.0,0.0,0.0);
          wgt[n]=0.0;
        }
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] != tgn)
          {
            for (i=0; i < mesh1->nt[b]; i++)
            {
              n0 = mesh1->t_n[b][i][0];
              n1 = mesh1->t_n[b][i][1];
              n2 = mesh1->t_n[b][i][2];
              if (tag[n0]!=0 && tag[n1]!=0 && tag[n2]!=0)
                continue;
              p0 = mesh1->nodep[n0];
              p1 = mesh1->nodep[n1];
              p2 = mesh1->nodep[n2];
              v1 = Vector(p0,p1);
              v2 = Vector(p0,p2);
              norm = v1 % v2;
              norm.normalize();
              if (tag[n0]==0 && tag[n1]==0)
              {
                v1 = Vector(p1,p0);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n0]);
                tng[n0] += v2*mag;
                wgt[n0] += mag;
                mag = 1.0-nrm[n0]*nrm[n1];
                tng[n0] += nrm[n1]*mag;
                wgt[n0] += mag;
                mag = 1.0 - fabs(v2 * nrm[n1]);
                tng[n1] += v2*mag;
                wgt[n1] += mag;
                mag = 1.0-nrm[n0]*nrm[n1];
                tng[n1] += nrm[n0]*mag;
                wgt[n1] += mag;
              }
              if (tag[n1]==0 && tag[n2]==0)
              {
                v1 = Vector(p2,p1);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n1]);
                tng[n1] += v2*mag;
                wgt[n1] += mag;
                mag = 1.0-nrm[n1]*nrm[n2];
                tng[n1] += nrm[n2]*mag;
                wgt[n1] += mag;
                mag = 1.0 - fabs(v2 * nrm[n2]);
                tng[n2] += v2*mag;
                wgt[n2] += mag;
                mag = 1.0-nrm[n1]*nrm[n2];
                tng[n2] += nrm[n1]*mag;
                wgt[n2] += mag;
              }
              if (tag[n2]==0 && tag[n0]==0)
              {
                v1 = Vector(p0,p2);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n2]);
                tng[n2] += v2*mag;
                wgt[n2] += mag;
                mag = 1.0-nrm[n0]*nrm[n2];
                tng[n2] += nrm[n0]*mag;
                wgt[n2] += mag;
                mag = 1.0 - fabs(v2 * nrm[n0]);
                tng[n0] += v2*mag;
                wgt[n0] += mag;
                mag = 1.0-nrm[n0]*nrm[n2];
                tng[n0] += nrm[n2]*mag;
                wgt[n0] += mag;
              }
            }
            for (i=0; i < mesh1->nq[b]; i++)
            {
              n0 = mesh1->q_n[b][i][0];
              n1 = mesh1->q_n[b][i][1];
              n2 = mesh1->q_n[b][i][2];
              n3 = mesh1->q_n[b][i][3];
              if (tag[n0]!=0 && tag[n1]!=0 && tag[n2]!=0 && tag[n3]!=0)
                continue;
              p0 = mesh1->nodep[n0];
              p1 = mesh1->nodep[n1];
              p2 = mesh1->nodep[n2];
              p3 = mesh1->nodep[n3];
              v1 = Vector(p0,p2);
              v2 = Vector(p1,p3);
              norm = v1 % v2;
              norm.normalize();
              if (tag[n0]==0 && tag[n1]==0)
              {
                v1 = Vector(p1,p0);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n0]);
                tng[n0] += v2*mag;
                wgt[n0] += mag;
                mag = 1.0-nrm[n0]*nrm[n1];
                tng[n0] += nrm[n1]*mag;
                wgt[n0] += mag;
                mag = 1.0 - fabs(v2 * nrm[n1]);
                tng[n1] += v2*mag;
                wgt[n1] += mag;
                mag = 1.0-nrm[n0]*nrm[n1];
                tng[n1] += nrm[n0]*mag;
                wgt[n1] += mag;
              }
              if (tag[n1]==0 && tag[n2]==0)
              {
                v1 = Vector(p2,p1);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n1]);
                tng[n1] += v2*mag;
                wgt[n1] += mag;
                mag = 1.0-nrm[n1]*nrm[n2];
                tng[n1] += nrm[n2]*mag;
                wgt[n1] += mag;
                mag = 1.0 - fabs(v2 * nrm[n2]);
                tng[n2] += v2*mag;
                wgt[n2] += mag;
                mag = 1.0-nrm[n1]*nrm[n2];
                tng[n2] += nrm[n1]*mag;
                wgt[n2] += mag;
              }
              if (tag[n2]==0 && tag[n3]==0)
              {
                v1 = Vector(p3,p2);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n2]);
                tng[n2] += v2*mag;
                wgt[n2] += mag;
                mag = 1.0-nrm[n3]*nrm[n2];
                tng[n2] += nrm[n3]*mag;
                wgt[n2] += mag;
                mag = 1.0 - fabs(v2 * nrm[n3]);
                tng[n3] += v2*mag;
                wgt[n3] += mag;
                mag = 1.0-nrm[n3]*nrm[n2];
                tng[n3] += nrm[n2]*mag;
                wgt[n3] += mag;
              }
              if (tag[n3]==0 && tag[n0]==0)
              {
                v1 = Vector(p0,p3);
                v2 = v1 % norm;
                v2.normalize();
                mag = 1.0 - fabs(v2 * nrm[n3]);
                tng[n3] += v2*mag;
                wgt[n3] += mag;
                mag = 1.0-nrm[n0]*nrm[n3];
                tng[n3] += nrm[n0]*mag;
                wgt[n3] += mag;
                mag = 1.0 - fabs(v2 * nrm[n0]);
                tng[n0] += v2*mag;
                wgt[n0] += mag;
                mag = 1.0-nrm[n0]*nrm[n3];
                tng[n0] += nrm[n3]*mag;
                wgt[n0] += mag;
              }
            }
          }
        }

        //don't include ghost normals in rms, wait for comm to reset

        rms = 0.0;
        int nrms = 0;
        if (knt > 0 || nsmoo > 0)
        {
          for (n=0; n < mesh1->nn; n++)
          {
            if (mesh1->pmap[n][1] == my_rank && tag[n] == 0 && wgt[n] > 1.0e-20)
            {
              if ((nopt <= nsmoo && ntag[n] == 0) || ntag[n] == 1)
              {
                tng[n].normalize();
                tng[n] -= nrm[n];
                nrm[n] += tng[n]*relax;
                dot = nrm[n] * bnorm[n];
                nrm[n] -= bnorm[n]*dot;
                nrm[n].normalize();
                //nrm[n] += tng[n]*relax;
                rms += tng[n]*tng[n];
                nrms++;
              }
            }
          }

          // impose sidewall tangency
          for (i=0; i < nfn; i++)
          {
            n = bf_node[i].index;
            if (mesh1->pmap[n][1] != my_rank)
              continue;
            norm = bf_node[i].norm;
            dot = norm * nrm[n];
            nrm[n] -= norm*dot;
            nrm[n].normalize();
          }
        }
        // set critical floating node to line up with chain
        for (b=0; b < mesh1->nb; b++)
        {
          if (!bft[b])
            continue;

          for (n=0; n < mesh1->nn; n++)
            itmp[n] = 0;
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            itmp[n0] = 1;
            itmp[n1] = 1;
            itmp[n2] = 1;
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            itmp[n0] = 1;
            itmp[n1] = 1;
            itmp[n2] = 1;
            itmp[n3] = 1;
          }
          for (j=0; j < mesh1->nb; j++)
          {
            if (!bft[j] || j==b)
              continue;
            for (i=0; i < mesh1->nt[j]; i++)
            {
              n0 = mesh1->t_n[j][i][0];
              n1 = mesh1->t_n[j][i][1];
              n2 = mesh1->t_n[j][i][2];
              if (itmp[n0] == 1) itmp[n0] = 2;
              if (itmp[n1] == 1) itmp[n1] = 2;
              if (itmp[n2] == 1) itmp[n2] = 2;
            }
            for (i=0; i < mesh1->nq[j]; i++)
            {
              n0 = mesh1->q_n[j][i][0];
              n1 = mesh1->q_n[j][i][1];
              n2 = mesh1->q_n[j][i][2];
              n3 = mesh1->q_n[j][i][3];
              if (itmp[n0] == 1) itmp[n0] = 2;
              if (itmp[n1] == 1) itmp[n1] = 2;
              if (itmp[n2] == 1) itmp[n2] = 2;
              if (itmp[n3] == 1) itmp[n3] = 2;
            }
            for (n=0; n < mesh1->nn; n++)
            {
              if (tag[n] != 0 || itmp[n] != 2)
                continue;
              for (i=0; i < nhash[n]->max; i++)
              {
                m = nhash[n]->list[i];
                if (tag[m] < 0 && itmp[m] == 2)
                {
                  nrm[n] = Vector(mesh1->nodep[n],mesh1->nodep[m]);
                  nrm[n].normalize();
                  break;
                }
              }
            }
          }
        }
        
        //exchange nrm after each pass here 
        mesh1->exchange_Vector(nrm);

        // compute global rms
        globalrms = rms;
#ifdef PARALLEL
        MPI_Allreduce(&rms,&globalrms,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
#endif
        rms = sqrt(rms/MAX(1,nrms));
        globalnrms = nrms;
#ifdef PARALLEL
        MPI_Allreduce(&nrms,&globalnrms,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
        globalrms = sqrt(globalrms/MAX(1,globalnrms));
        globalknt = knt;
#ifdef PARALLEL
        MPI_Allreduce(&knt,&globalknt,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);        
#endif

        if (my_rank == 0)
        {
          fprintf(out_f,"\nLayer %d, normal optimizing loop %d, residual = %g, # of negative dot products = %d",lay,nopt,globalrms,gint);
          fflush(out_f);
        }
      }
      if (gint > 0)
      {
        if (my_rank == 0)
        {
          //fprintf(out_f,"\nNormals contain negative dot products. Cannot continue!");
          fprintf(out_f,"\nNormals contain negative dot products. Restricting layers for those points!");
          fflush(out_f);
        }
        //MPI_Abort(MPI_COMM_WORLD,0);
        //exit(0);
      }

      delete[] bnorm;

      // compute desired wall spacing
      // adjust for convex/concave situations
      for (n=0; n < mesh1->nn; n++)
      {
        desired[n] = 0.0;
        wgt[n] = 0.0;
      }
      //dotmax = cos(atan(1.0));
      //dotmax = 0.5;
      for (b=0; b < mesh1->nb; b++)
      {
        if (geom->layers[b] == tgn)
        {
          ds = geom->n_space[b];
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];

            desired[n0] += ds;
            wgt[n0] += 1.0;
            desired[n1] += ds;
            wgt[n1] += 1.0;
            desired[n2] += ds;
            wgt[n2] += 1.0;

            //p0 = mesh1->nodep[n0];
            //p1 = mesh1->nodep[n1];
            //p2 = mesh1->nodep[n2];
            //v1 = Vector(p0,p1);
            //v2 = Vector(p0,p2);
            //norm = v1 % v2;
            //norm.normalize();

            //dot = MAX(dotmax,nrm[n0]*norm);
            //desired[n0] += ds/dot;
            //wgt[n0] += 1.0;
            //dot = MAX(dotmax,nrm[n1]*norm);
            //desired[n1] += ds/dot;
            //wgt[n1] += 1.0;
            //dot = MAX(dotmax,nrm[n2]*norm);
            //desired[n2] += ds/dot;
            //wgt[n2] += 1.0;

            //v1 = Vector(p0,p1);
            //v1.normalize();
            //v2 = Vector(p0,p2);
            //v2.normalize();
            //ang = acos(v1*v2);
            //v1 = Vector(p0,p1);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n0]));
            //desired[n0] += ds*dot*ang;
            //wgt[n0] += ang;
            //v1 = Vector(p0,p2);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n0]));
            //desired[n0] += ds*dot*ang;
            //wgt[n0] += ang;

            //v1 = Vector(p1,p0);
            //v1.normalize();
            //v2 = Vector(p1,p2);
            //v2.normalize();
            //ang = acos(v1*v2);
            //v1 = Vector(p1,p0);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n1]));
            //desired[n1] += ds*dot*ang;
            //wgt[n1] += ang;
            //v1 = Vector(p1,p2);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n1]));
            //desired[n1] += ds*dot*ang;
            //wgt[n1] += ang;

            //v1 = Vector(p2,p0);
            //v1.normalize();
            //v2 = Vector(p2,p1);
            //v2.normalize();
            //ang = acos(v1*v2);
            //v1 = Vector(p2,p0);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n2]));
            //desired[n2] += ds*dot*ang;
            //wgt[n2] += ang;
            //v1 = Vector(p2,p1);
            //v1.normalize();
            //dot = (1.0 + (v1 * nrm[n2]));
            //desired[n2] += ds*dot*ang;
            //wgt[n2] += ang;
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];

            desired[n0] += ds;
            wgt[n0] += 1.0;
            desired[n1] += ds;
            wgt[n1] += 1.0;
            desired[n2] += ds;
            wgt[n2] += 1.0;
            desired[n3] += ds;
            wgt[n3] += 1.0;

            //p0 = mesh1->nodep[n0];
            //p1 = mesh1->nodep[n1];
            //p2 = mesh1->nodep[n2];
            //p3 = mesh1->nodep[n3];
            //v1 = Vector(p0,p2);
            //v2 = Vector(p1,p3);
            //norm = v1 % v2;
            //norm.normalize();

            //dot = MAX(dotmax,nrm[n0]*norm);
            //desired[n0] += ds/dot;
            //wgt[n0] += 1.0;
            //dot = MAX(dotmax,nrm[n1]*norm);
            //desired[n1] += ds/dot;
            //wgt[n1] += 1.0;
            //dot = MAX(dotmax,nrm[n2]*norm);
            //desired[n2] += ds/dot;
            //wgt[n2] += 1.0;
            //dot = MAX(dotmax,nrm[n3]*norm);
            //desired[n3] += ds/dot;
            //wgt[n3] += 1.0;

            //v1 = Vector(p0,p1);
            //v1.normalize();
            //v2 = Vector(p0,p3);
            //v2.normalize();
            //ang = acos(v1*v2);
            //dot = (1.0 + (v1 * nrm[n0]));
            //desired[n0] += ds*dot*ang;
            //wgt[n0] += ang;
            //dot = (1.0 + (v2 * nrm[n0]));
            //desired[n0] += ds*dot*ang;
            //wgt[n0] += ang;

            //v1 = Vector(p1,p0);
            //v1.normalize();
            //v2 = Vector(p1,p2);
            //v2.normalize();
            //ang = acos(v1*v2);
            //dot = (1.0 + (v1 * nrm[n1]));
            //desired[n1] += ds*dot*ang;
            //wgt[n1] += ang;
            //dot = (1.0 + (v2 * nrm[n1]));
            //desired[n1] += ds*dot*ang;
            //wgt[n1] += ang;

            //v1 = Vector(p2,p1);
            //v1.normalize();
            //v2 = Vector(p2,p3);
            //v2.normalize();
            //ang = acos(v1*v2);
            //dot = (1.0 + (v1 * nrm[n2]));
            //desired[n2] += ds*dot*ang;
            //wgt[n2] += ang;
            //dot = (1.0 + (v2 * nrm[n2]));
            //desired[n2] += ds*dot*ang;
            //wgt[n2] += ang;

            //v1 = Vector(p3,p0);
            //v1.normalize();
            //v2 = Vector(p3,p2);
            //v2.normalize();
            //ang = acos(v1*v2);
            //dot = (1.0 + (v1 * nrm[n3]));
            //desired[n3] += ds*dot*ang;
            //wgt[n3] += ang;
            //dot = (1.0 + (v2 * nrm[n3]));
            //desired[n3] += ds*dot*ang;
            //wgt[n3] += ang;
          }
        }
      }

      for (i=0; i < punt->max; i++)
      {
        n = punt->list[i];
        desired[n] = 1.0e20;
      }
      delete punt;

      double dsmin, dsmax;
      dsmin = 1.0e20;
      dsmax = 0.0;
      for (n=0; n < mesh1->nn; n++)
      {
        if (mesh1->pmap[n][1] == my_rank && wgt[n] > 1.0e-15)
        {
          desired[n] /= wgt[n];
          dsmin = MIN(dsmin,desired[n]);
          dsmax = MAX(dsmax,desired[n]);
        }
      }
      gdouble = dsmin;
#ifdef PARALLEL
      MPI_Allreduce(&dsmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);        
#endif
      if (my_rank == 0) fprintf(out_f,"\nMinimum initial marching distance = %lf",gdouble);
      gdouble = dsmax;
#ifdef PARALLEL
      MPI_Allreduce(&dsmax,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);        
#endif
      if (my_rank == 0) fprintf(out_f,"\nMaximum initial marching distance = %lf",gdouble);
      if (my_rank == 0) fflush(out_f);

#ifdef PARALLEL
      //exchange desired distances
      mesh1->exchange_double(desired);
#endif

      for (n=0; n < mesh1->nn; n++)
        if (tag[n] == 0) tag[n] = lay;

#ifdef PARALLEL
      m=0;
      mesh1->exchange_tag(tag,m);
#endif

      //
      // create hash tables for elements of marching boundary nodes only
      //
      List **tethash, **pyrhash, **prihash, **hexhash;
      tethash = new List*[mesh1->nn];
      pyrhash = new List*[mesh1->nn];
      prihash = new List*[mesh1->nn];
      hexhash = new List*[mesh1->nn];
      for (n=0; n < mesh1->nn; n++)
      {
        tethash[n] = new List();
        pyrhash[n] = new List();
        prihash[n] = new List();
        hexhash[n] = new List();
      }

      for (c=0; c < mesh1->ntet; c++)
        for (i=0; i < 4; i++)
        {
          n=mesh1->tet_n[c][i];
          if (tag[n] >= 0)
            tethash[n]->Add_To_List(c);
        }
      for (c=0; c < mesh1->npyr; c++)
        for (i=0; i < 5; i++)
        {
          n=mesh1->pyr_n[c][i];
          if (tag[n] >= 0)
            pyrhash[n]->Add_To_List(c);
        }
      for (c=0; c < mesh1->npri; c++)
        for (i=0; i < 6; i++)
        {
          n=mesh1->pri_n[c][i];
          if (tag[n] >= 0)
            prihash[n]->Add_To_List(c);
        }
      for (c=0; c < mesh1->nhex; c++)
        for (i=0; i < 8; i++)
        {
          n=mesh1->hex_n[c][i];
          if (tag[n] >= 0)
            hexhash[n]->Add_To_List(c);
        }

      for (n=0; n < mesh1->nn; n++)
      {
        if (tag[n] <= 0)
          continue;

        // compute minimum distance to interior nodes
        double dtotal, tfac;
        dtotal = stack_height(lay,geo1,geo2,geo_min,geo_max,tfac);

        if (v_layers < 0)
          mag = dtotal*aao_fac;
        else
          mag = tfac;

        ds = mag*desired[n];
        for (i=0; i < tethash[n]->max; i++)
        {
          c = tethash[n]->list[i];
          for (j=0; j < 4; j++)
            if ((m = mesh1->tet_n[c][j]) != n && (tag[m] == -2 || tag[m] == -1))
              ds = MIN(ds,distance(mesh1->nodep[n],mesh1->nodep[m]));
          pa = mesh1->nodep[n];
          pb = pa + Point(nrm[n][0],nrm[n][1],nrm[n][2])*ds;
          p0 = mesh1->nodep[n0=mesh1->tet_n[c][0]];
          p1 = mesh1->nodep[n1=mesh1->tet_n[c][1]];
          p2 = mesh1->nodep[n2=mesh1->tet_n[c][2]];
          p3 = mesh1->nodep[n3=mesh1->tet_n[c][3]];

          if (n == n0 && line_facet_intersect(pa,pb,p1,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n1 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n2 && line_facet_intersect(pa,pb,p0,p1,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n3 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
        }
        for (i=0; i < pyrhash[n]->max; i++)
        {
          c = pyrhash[n]->list[i];
          for (j=0; j < 5; j++)
            if ((m = mesh1->pyr_n[c][j]) != n && (tag[m] == -2 || tag[m] == -1))
              ds = MIN(ds,distance(mesh1->nodep[n],mesh1->nodep[m]));
          pa = mesh1->nodep[n];
          pb = pa + Point(nrm[n][0],nrm[n][1],nrm[n][2])*ds;
          p0 = mesh1->nodep[n0=mesh1->pyr_n[c][0]];
          p1 = mesh1->nodep[n1=mesh1->pyr_n[c][1]];
          p2 = mesh1->nodep[n2=mesh1->pyr_n[c][2]];
          p3 = mesh1->nodep[n3=mesh1->pyr_n[c][3]];
          p4 = mesh1->nodep[n4=mesh1->pyr_n[c][4]];

          if (n == n0 && line_facet_intersect(pa,pb,p1,p2,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p2,p3,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n1 && line_facet_intersect(pa,pb,p2,p3,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n2 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p0,p1,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n3 && line_facet_intersect(pa,pb,p0,p1,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p1,p2,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n4 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
        }
        for (i=0; i < prihash[n]->max; i++)
        {
          c = prihash[n]->list[i];
          for (j=0; j < 6; j++)
            if ((m = mesh1->pri_n[c][j]) != n && (tag[m] == -2 || tag[m] == -1))
              ds = MIN(ds,distance(mesh1->nodep[n],mesh1->nodep[m]));
          pa = mesh1->nodep[n];
          pb = pa + Point(nrm[n][0],nrm[n][1],nrm[n][2])*ds;
          p0 = mesh1->nodep[n0=mesh1->pri_n[c][0]];
          p1 = mesh1->nodep[n1=mesh1->pri_n[c][1]];
          p2 = mesh1->nodep[n2=mesh1->pri_n[c][2]];
          p3 = mesh1->nodep[n3=mesh1->pri_n[c][3]];
          p4 = mesh1->nodep[n4=mesh1->pri_n[c][4]];
          p5 = mesh1->nodep[n5=mesh1->pri_n[c][5]];

          if (n == n0 && line_facet_intersect(pa,pb,p3,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p1,p2,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p2,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n1 && line_facet_intersect(pa,pb,p3,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p2,p0,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p2,p3,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n2 && line_facet_intersect(pa,pb,p3,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p0,p1,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p1,p3,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n3 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p1,p2,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p2,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n4 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p2,p0,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p2,p3,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n5 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p0,p1,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p1,p3,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
        }
        for (i=0; i < hexhash[n]->max; i++)
        {
          c = hexhash[n]->list[i];
          for (j=0; j < 8; j++)
            if ((m = mesh1->hex_n[c][j]) != n && (tag[m] == -2 || tag[m] == -1))
              ds = MIN(ds,distance(mesh1->nodep[n],mesh1->nodep[m]));
          pa = mesh1->nodep[n];
          pb = pa + Point(nrm[n][0],nrm[n][1],nrm[n][2])*ds;
          p0 = mesh1->nodep[n0=mesh1->hex_n[c][0]];
          p1 = mesh1->nodep[n1=mesh1->hex_n[c][1]];
          p2 = mesh1->nodep[n2=mesh1->hex_n[c][2]];
          p3 = mesh1->nodep[n3=mesh1->hex_n[c][3]];
          p4 = mesh1->nodep[n4=mesh1->hex_n[c][4]];
          p5 = mesh1->nodep[n5=mesh1->hex_n[c][5]];
          p6 = mesh1->nodep[n6=mesh1->hex_n[c][6]];
          p7 = mesh1->nodep[n7=mesh1->hex_n[c][7]];

          if (n == n0 && line_facet_intersect(pa,pb,p1,p2,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p2,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p2,p3,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p3,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p4,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n0 && line_facet_intersect(pa,pb,p4,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n1 && line_facet_intersect(pa,pb,p2,p3,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p3,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p3,p7,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p4,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n1 && line_facet_intersect(pa,pb,p4,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n2 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p3,p7,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p0,p1,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p0,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p4,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n2 && line_facet_intersect(pa,pb,p4,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n3 && line_facet_intersect(pa,pb,p0,p1,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p0,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p1,p2,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p2,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p4,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n3 && line_facet_intersect(pa,pb,p4,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n4 && line_facet_intersect(pa,pb,p1,p2,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p2,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p2,p3,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p3,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n4 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n5 && line_facet_intersect(pa,pb,p2,p3,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p3,p6,p7,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p3,p7,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n5 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n6 && line_facet_intersect(pa,pb,p3,p0,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n6 && line_facet_intersect(pa,pb,p3,p7,p4,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n6 && line_facet_intersect(pa,pb,p0,p1,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n6 && line_facet_intersect(pa,pb,p0,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n6 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n6 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));

          if (n == n7 && line_facet_intersect(pa,pb,p0,p1,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n7 && line_facet_intersect(pa,pb,p0,p4,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n7 && line_facet_intersect(pa,pb,p1,p2,p5,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n7 && line_facet_intersect(pa,pb,p2,p5,p6,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n7 && line_facet_intersect(pa,pb,p0,p1,p2,pt,tol))
            ds = MIN(ds,distance(pa,pt));
          if (n == n7 && line_facet_intersect(pa,pb,p0,p2,p3,pt,tol))
            ds = MIN(ds,distance(pa,pt));
        }

        dtotal = tfac = 1.0;
        geo = geo1;
        for (i=1; i < lay; i++)
        {
          tfac *= geo;
          geo *= geo2;
          geo = MAX(geo_min,MIN(geo,geo_max));
          dtotal += tfac;

          if (v_layers < 0)
            mag = dtotal/aao_fac;
          else
            mag = tfac;

          if (mag*desired[n] > ds)
          {
            tag[n] = i-1;
            break;
          }
        }
      }
      m=1;
      mesh1->exchange_tag(tag,m);

      for (n=0; n < mesh1->nn; n++)
      {
        delete tethash[n];
        delete pyrhash[n];
        delete prihash[n];
        delete hexhash[n];
      }
      delete[] tethash;
      delete[] pyrhash;
      delete[] prihash;
      delete[] hexhash;
      
      k=abs(v_layers);
      m=0;
      for (n=0; n < mesh1->nn; n++)
        if (tag[n] >= 0)
        {
          k = MIN(k,tag[n]);
          m = MAX(m,tag[n]);
        }

      gint = k;
#ifdef PARALLEL
      MPI_Allreduce(&k,&gint,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nPreliminary minimum number of layers = %d",gint);
      gint = m;
#ifdef PARALLEL
      MPI_Allreduce(&m,&gint,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nPreliminary maximum number of layers = %d",gint);
      if (my_rank == 0) fflush(out_f);

      // perform geometry intersection test
      for (n=0; n < mesh1->nn; n++)
        ntag[n] = 0;
      for (i=0; i < nfn; i++)
        ntag[bf_node[i].index] = 1;
      List ilist;
      ilist.Redimension(0);
      for (n=0; n < mesh1->nn; n++)
      {
        if (tag[n] > 0 && ntag[0] == 0)
          ilist.Add_To_List(n);
      }
      gint = m = ilist.max;
#ifdef PARALLEL
      MPI_Allreduce(&m,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0)
      {
        fprintf(out_f,"\nLayer %d, performing intersection test for %d nodes.",lay,gint);
        fflush(out_f);
      }
      int tenth = (int)((double)ilist.max/10.0+1);
      int tens[11];
      for (n=0; n < 11; n++)
        tens[n] = MIN(MAX(0,ilist.max-1),n*tenth);
      int kprint = 0;
      for (k=l=0; l < 10; l++)
      {
        for (m=tens[l]; m < tens[l+1]; m++)
        {
          n = ilist.list[m];
          mag = 0.0;
          for (i=0; i < nhash[n]->max; i++)
            mag += distance(mesh1->nodep[n],mesh1->nodep[nhash[n]->list[i]]);
          mag /= nhash[n]->max;
          mag *= 0.1;

          int iflag = 1;
          double dtotal = 3.0*stack_height(tag[n],geo1,geo2,geo_min,geo_max,tfac);
          p0 = mesh1->nodep[n] + Point(nrm[n][0],nrm[n][1],nrm[n][2])*mag;
          p1 = p0 + Point(nrm[n][0],nrm[n][1],nrm[n][2])*dtotal*desired[n];
          if ((iflag = geom->Intersect(p0,p1,p2,0)) == 1 && distance(mesh1->nodep[n],p2) > 1e-12)
          {
            k++;
            ds = distance(mesh1->nodep[n],p2)/3.0;
            while (iflag && tag[n] > 0)
            {
              tag[n]--;
              dtotal = stack_height(tag[n],geo1,geo2,geo_min,geo_max,tfac);
              if (dtotal*desired[n] > ds)
                iflag=1;
              else
                iflag=0;
            }
          
            if (num_procs == 1 && kprint < 10)
            {
              fprintf(out_f,"\nIntersection detected from (%lg,%lg,%lg) to (%lg,%lg,%lg)",
                      p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
              fprintf(out_f,"\n  mag = %lg, dtotal = %lg, desired = %lg",mag,dtotal,desired[n]);
              fprintf(out_f,"\n  Normal vector = (%lg,%lg,%lg)",nrm[n][0],nrm[n][1],nrm[n][2]);
              fprintf(out_f,"\n  Intersection point = (%lg,%lg,%lg)",p2[0],p2[1],p2[2]);
              for (int bb=0; bb < geom->ngb; bb++)
              {
                if (geom->Intersect(p0,p1,p2,bb) == 1)
                  fprintf(out_f,"\n  Intersection detected with boundary %d",bb+1);
              }
              fflush(out_f);
              kprint++;
            }
          }
        }
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (my_rank == 0) fprintf(out_f,"\n%d%%",10*(l+1));
        if (my_rank == 0) fflush(out_f);
      }
      ilist.Redimension(0);
      gint = k;
#ifdef PARALLEL
      MPI_Allreduce(&k,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0 && gint > 0)
        fprintf(out_f,"\nNumber of intersections detected = %d",gint);
      if (my_rank == 0) fflush(out_f);

      j = 0;
      gint = nfn;
#ifdef PARALLEL
      MPI_Allreduce(&nfn,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0)
      {
        fprintf(out_f,"\nLayer %d, performing intersection test for %d sidewall nodes.\n",lay,gint);
        fflush(out_f);
      }
      for (i=0; i < nfn; i++)
      {
        n = bf_node[i].index;
        if (tag[n] >= 0)
          tng[n] = Vector(0.0,0.0,0.0);
      }
      for (i=0; i < nfn; i++)
      {
        n = bf_node[i].index;
        if (tag[n] >= 0)
          tng[n] += bf_node[i].norm;
      }
      kprint = 0;
      for (k=m=0; m < nfn; m++)
      {
        n = bf_node[m].index;
        if (tag[n] >= 0)
        {
          //norm = bf_node[i].norm;
          vec = tng[n];
          vec.normalize();
          vec *= 0.025;
          mag = 0.0;
          for (i=0; i < nhash[n]->max; i++)
            mag += distance(mesh1->nodep[n],mesh1->nodep[nhash[n]->list[i]]);
          mag /= nhash[n]->max;
          mag *= 0.1;

          int iflag = 1;
          double dtotal = 3.0*stack_height(tag[n],geo1,geo2,geo_min,geo_max,tfac);
          v1 = nrm[n] + vec;
          v1.normalize();
          p0 = mesh1->nodep[n] + Point(v1[0],v1[1],v1[2])*mag;
          p1 = p0 + Point(v1[0],v1[1],v1[2])*dtotal*desired[n];
          if ((iflag = geom->Intersect(p0,p1,p2,0)) == 1 && distance(mesh1->nodep[n],p2) > 1e-12)
          {
            k++;
            ds = distance(mesh1->nodep[n],p2)/3.0;
            while (iflag && tag[n] > 0)
            {
              tag[n]--;
              dtotal = stack_height(tag[n],geo1,geo2,geo_min,geo_max,tfac);
              if (dtotal*desired[n] > ds)
                iflag=1;
              else
                iflag=0;
            }
        
            if (num_procs == 1 && kprint < 10)
            {
              fprintf(out_f,"\nIntersection detected from (%lg,%lg,%lg) to (%lg,%lg,%lg)",
                      p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
              fprintf(out_f,"\n  mag = %lg, dtotal = %lg, desired = %lg",mag,dtotal,desired[n]);
              fprintf(out_f,"\n  Normal vector = (%lg,%lg,%lg)",nrm[n][0],nrm[n][1],nrm[n][2]);
              fprintf(out_f,"\n  Intersection point = (%lg,%lg,%lg)",p2[0],p2[1],p2[2]);
              for (int bb=0; bb < geom->ngb; bb++)
              {
                if (geom->Intersect(p0,p1,p2,bb) == 1)
                  fprintf(out_f,"\n  Intersection detected with boundary %d",bb+1);
              }
              fflush(out_f);
              kprint++;
            }
          }
        }
      }
      gint = k;
#ifdef PARALLEL
      MPI_Allreduce(&k,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
      if (my_rank == 0 && gint > 0)
        fprintf(out_f,"\nNumber of sidewall node intersections detected = %d",gint);
      if (my_rank == 0) fflush(out_f);

      if (nfn > 0)
        delete[] bf_node;
      delete[] ntag;
      delete[] tng;
      delete[] wgt;
      delete[] itmp;
      for (n=0; n < mesh1->nn; n++)
        delete nhash[n];
      delete[] nhash;
     
      // force one layer to be inserted
      if (v_layers == 1)
        for (n=0; n < mesh1->nn; n++)
          if (tag[n] == 0) tag[n] = 1;
      
#ifdef PARALLEL
      //exchange for proper count
      m=0;
      mesh1->exchange_tag(tag,m);
#endif

      m = lay;
      n = 0;
      for (i=0; i < mesh1->nn; i++)
        if (tag[i] >= 0)
        {
          m = MIN(m,tag[i]);
          n = MAX(n,tag[i]);
        }

      gint = m;
#ifdef PARALLEL
      MPI_Allreduce(&m,&gint,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nOriginal minimum number of layers = %d",gint);
      gint = n;
#ifdef PARALLEL
      MPI_Allreduce(&n,&gint,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
      if (my_rank == 0) fprintf(out_f,"\nOriginal maximum number of layers = %d",gint);
      if (my_rank == 0) fprintf(out_f,"\n");
      if (my_rank == 0) fflush(out_f);

      // check element quality
      double vol;

      tol = 1.0e-15;

      int sweep = 0;
      int globali = 0;
      int globalj = 0;
      int globalk = 0;
      int globalm = 0;
      int globaln = 0;
      int itotal= 0;
      int jtotal= 0;
      int ktotal= 0;
      int vtest = 1;
      int jtest = 1;

      int ldiff = 1;
      if (v_layers < 0)
        ldiff = 2*stack_growth;
      else
        ldiff = stack_growth;
      do
      {
        k=1;
        mesh1->exchange_tag(tag,k);
        sweep++;

        if (v_layers == 1 && lay == 1)
          break;

        if (v_layers < 0)
          vtest = jtest = 1;

        i=j=k=n=0;

        l=0;
        do
        {
          k=0;
          for (b=0; b < mesh1->nb; b++)
          {
            if (geom->layers[b] == tgn)
            {
              for (c=0; c < mesh1->nt[b]; c++)
              {
                n0 = mesh1->t_n[b][c][0];
                n1 = mesh1->t_n[b][c][1];
                n2 = mesh1->t_n[b][c][2];
                // only allow ldiff layers difference between nodes
                int mn = MIN(tag[n0],MIN(tag[n1],tag[n2]));
                int mx = MAX(tag[n0],MAX(tag[n1],tag[n2]));
                if (mx-mn <= ldiff)
                  continue;
                // knock down high node
                if (tag[n0] > mn+ldiff)
                {
                  tag[n0]--;
                  k++;
                }
                if (tag[n1] > mn+ldiff)
                {
                  tag[n1]--;
                  k++;
                }
                if (tag[n2] > mn+ldiff)
                {
                  tag[n2]--;
                  k++;
                }
              }
              for (c=0; c < mesh1->nq[b]; c++)
              {
                n0 = mesh1->q_n[b][c][0];
                n1 = mesh1->q_n[b][c][1];
                n2 = mesh1->q_n[b][c][2];
                n3 = mesh1->q_n[b][c][3];

                // allow all with same number of layers
                if (tag[n0]==tag[n1] && tag[n2]==tag[n3] && tag[n0]==tag[n2])
                  continue;

                // allow with same number of layers on opposite edges
                if ((tag[n0]==tag[n1] && tag[n2]==tag[n3]) ||
                    (tag[n0]==tag[n3] && tag[n1]==tag[n2]))
                  continue;

                int count, hi_side, hi_count, m01, m23, m03, m12;
                hi_side = 0;
                hi_count = (tag[n0]+tag[n1]);
                if ((count = tag[n1]+tag[n2]) > hi_count)
                {
                  hi_count = count;
                  hi_side = 1;
                }
                if ((count = tag[n2]+tag[n3]) > hi_count)
                {
                  hi_count = count;
                  hi_side = 2;
                }
                if ((count = tag[n3]+tag[n0]) > hi_count)
                {
                  hi_count = count;
                  hi_side = 3;
                }
                // knock down high node on high and opposite sides
                switch(hi_side)
                {
                  case 0:
                  case 2:
                    //m01 = (tag[n0]+tag[n1])/2;
                    //m23 = (tag[n2]+tag[n3])/2;
                    m01 = MIN(tag[n0],tag[n1]);
                    m23 = MIN(tag[n2],tag[n3]);
                    //m01 = MAX(tag[n0],tag[n1]);
                    //m23 = MAX(tag[n2],tag[n3]);
                    if (tag[n0] != m01)
                    {
                      tag[n0] = (tag[n0] > m01) ? tag[n0]-1 : tag[n0]+1;
                      k++;
                    }
                    if (tag[n1] != m01)
                    {
                      tag[n1] = (tag[n1] > m01) ? tag[n1]-1 : tag[n1]+1;
                      k++;
                    }
                    if (tag[n2] != m23)
                    {
                      tag[n2] = (tag[n2] > m23) ? tag[n2]-1 : tag[n2]+1;
                      k++;
                    }
                    if (tag[n3] != m23)
                    {
                      tag[n3] = (tag[n3] > m23) ? tag[n3]-1 : tag[n3]+1;
                      k++;
                    }
                    break;
                  case 1:
                  case 3:
                    //m03 = (tag[n0]+tag[n3])/2;
                    //m12 = (tag[n1]+tag[n2])/2;
                    m03 = MIN(tag[n0],tag[n3]);
                    m12 = MIN(tag[n1],tag[n2]);
                    //m03 = MAX(tag[n0],tag[n3]);
                    //m12 = MAX(tag[n1],tag[n2]);
                    if (tag[n0] != m03)
                    {
                      tag[n0] = (tag[n0] > m03) ? tag[n0]-1 : tag[n0]+1;
                      k++;
                    }
                    if (tag[n3] != m03)
                    {
                      tag[n3] = (tag[n3] > m03) ? tag[n3]-1 : tag[n3]+1;
                      k++;
                    }
                    if (tag[n1] != m12)
                    {
                      tag[n1] = (tag[n1] > m12) ? tag[n1]-1 : tag[n1]+1;
                      k++;
                    }
                    if (tag[n2] != m12)
                    {
                      tag[n2] = (tag[n2] > m12) ? tag[n2]-1 : tag[n2]+1;
                      k++;
                    }
                    break;
                }
              }
            }
          }

          l += k;
        
          //continue in loop until all have finished element quality checking      
          globalk = k;
#ifdef PARALLEL
          MPI_Allreduce(&k,&globalk,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);         
#endif
        } while(globalk > 0);
        k = l;

        negmax = 1.0e20;
        for (c=0; c < mesh1->ntet && vtest; c++)
        {
          n0 = mesh1->tet_n[c][0];
          n1 = mesh1->tet_n[c][1];
          n2 = mesh1->tet_n[c][2];
          n3 = mesh1->tet_n[c][3];
          // elements must have positive volume
          if ((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0)) ||
              (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0))
          {
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            if ((vol=tetrahedral_volume(p0,p1,p2,p3)) < 0.0)
            {
              n++;
              if (vol < negmax)
              {
                negmax = vol;
                nmax = c;
              }
            } else
            {
              double vol_tol = vol*0.05;
              vol_tol = 0.0;

              p4 = p0 = mesh1->nodep[n0];
              p5 = p1 = mesh1->nodep[n1];
              p6 = p2 = mesh1->nodep[n2];
              p7 = p3 = mesh1->nodep[n3];
              tfac = 1.0;
              geo = geo1;
              for (l=0; l < lay; l++)
              {
                if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2] && l >= tag[n3])
                  break;

                if (v_layers < 0)
                {
                  p0 = p4;
                  p1 = p5;
                  p2 = p6;
                  p3 = p7;
                }

                if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                {
                  vec = nrm[n0];
                  p4 = p0 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n0];
                }
                if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                {
                  vec = nrm[n1];
                  p5 = p1 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n1];
                }
                if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                {
                  vec = nrm[n2];
                  p6 = p2 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n2];
                }
                if ((l < tag[n3] && v_layers < 0) || (v_layers > 0 && l == tag[n3]-1 && l == lay-1))
                {
                  vec = nrm[n3];
                  p7 = p3 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n3];
                }
                if ((v_layers < 0 || l == lay-1) && (vol=tetrahedral_volume(p4,p5,p6,p7)) < vol_tol)
                {
                  if (tag[n0] > l)
                  {
                    tag[n0] = l;
                    i++;
                  }
                  if (tag[n1] > l)
                  {
                    tag[n1] = l;
                    i++;
                  }
                  if (tag[n2] > l)
                  {
                    tag[n2] = l;
                    i++;
                  }
                  if (tag[n3] > l)
                  {
                    tag[n3] = l;
                    i++;
                  }
                }
                tfac *= geo;
                geo *= geo2;
                geo = MAX(geo_min,MIN(geo,geo_max));
              }
            }
          }
        }
        gdouble = negmax;
#ifdef PARALLEL
        MPI_Allreduce(&negmax,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);         
#endif
        if (gdouble < 0.0)
        {
          cg = Point(0.0,0.0,0.0);
          if (nmax >= 0)
          {
            n0 = mesh1->tet_n[nmax][0];
            n1 = mesh1->tet_n[nmax][1];
            n2 = mesh1->tet_n[nmax][2];
            n3 = mesh1->tet_n[nmax][3];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            cg = (p0+p1+p2+p3)*0.25;
          }
#ifdef PARALLEL
          cmx_g_in.val = negmax;
          cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&(cg[0]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[1]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[2]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&cg[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif
          if (my_rank == 0)
          {
            fprintf(out_f,"\n  Minimum tetrahedral volume = %lg @ (%lg, %lg, %lg)",gdouble,cg[0],cg[1],cg[2]);
            fflush(out_f);
          }
        }
        
        negmax = 1.0e20;
        for (c=0; c < mesh1->npyr && vtest; c++)
        {
          n0 = mesh1->pyr_n[c][0];
          n1 = mesh1->pyr_n[c][1];
          n2 = mesh1->pyr_n[c][2];
          n3 = mesh1->pyr_n[c][3];
          n4 = mesh1->pyr_n[c][4];
          if (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n4]>=0) // Test corner volume if 4 nodes marching
          {
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            if ((vol=pyramid_volume(p0,p1,p2,p3,p4)) < 0.0)
            {
              n++;
              if (vol < negmax)
              {
                negmax = vol;
                nmax = c;
              }
            } else
            {
              double vol_tol = vol*0.05;
              vol_tol = 0.0;

              p5 = p0 = mesh1->nodep[n0];
              p6 = p1 = mesh1->nodep[n1];
              p7 = p2 = mesh1->nodep[n2];
              p8 = p3 = mesh1->nodep[n3];
              p9 = p4 = mesh1->nodep[n4];
              tfac = 1.0;
              geo = geo1;
              for (l=0; l < lay; l++)
              {
                if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2] && l >= tag[n3] && l >= tag[n4])
                  break;

                if (v_layers < 0)
                {
                  p0 = p5;
                  p1 = p6;
                  p2 = p7;
                  p3 = p8;
                  p4 = p9;
                }

                if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                {
                  vec = nrm[n0];
                  p5 = p0 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n0];
                }
                if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                {
                  vec = nrm[n1];
                  p6 = p1 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n1];
                }
                if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                {
                  vec = nrm[n2];
                  p7 = p2 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n2];
                }
                if ((l < tag[n3] && v_layers < 0) || (v_layers > 0 && l == tag[n3]-1 && l == lay-1))
                {
                  vec = nrm[n3];
                  p8 = p3 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n3];
                }
                if ((l < tag[n4] && v_layers < 0) || (v_layers > 0 && l == tag[n4]-1 && l == lay-1))
                {
                  vec = nrm[n4];
                  p9 = p4 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n4];
                }
                if (v_layers < 0 || l == lay-1)
                {
                  if (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0)
                  {
                    if ((vol=pyramid_volume(p5,p6,p7,p8,p9)) < vol_tol)
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                  } else
                  {
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n3]>=0 || tag[n4]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n3]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p5,p6,p8,p9)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n4]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p6,p7,p5,p9)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n4]>=0)) ||
                         (tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p7,p8,p6,p9)) < 0.0))
                    {
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n4]>=0)) ||
                         (tag[n0]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p8,p5,p7,p9)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                  }
                }
                tfac *= geo;
                geo *= geo2;
                geo = MAX(geo_min,MIN(geo,geo_max));
              }
            }
          }
        }
        gdouble = negmax;
#ifdef PARALLEL
        MPI_Allreduce(&negmax,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);         
#endif
        if (gdouble < 0.0)
        {
          cg = Point(0.0,0.0,0.0);
          if (nmax >= 0)
          {
            n0 = mesh1->pyr_n[c][0];
            n1 = mesh1->pyr_n[c][1];
            n2 = mesh1->pyr_n[c][2];
            n3 = mesh1->pyr_n[c][3];
            n4 = mesh1->pyr_n[c][4];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            cg = (p0+p1+p2+p3+p4)/5.0;
          }
#ifdef PARALLEL
          cmx_g_in.val = negmax;
          cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&(cg[0]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[1]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[2]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&cg[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif
          if (my_rank == 0)
          {
            fprintf(out_f,"\n  Minimum pyramid volume = %lg @ (%lg, %lg, %lg)",gdouble,cg[0],cg[1],cg[2]);
            fflush(out_f);
          }
        }
        
        negmax = 1.0e20;
        for (c=0; c < mesh1->npri && vtest; c++)
        {
          n0 = mesh1->pri_n[c][0];
          n1 = mesh1->pri_n[c][1];
          n2 = mesh1->pri_n[c][2];
          n3 = mesh1->pri_n[c][3];
          n4 = mesh1->pri_n[c][4];
          n5 = mesh1->pri_n[c][5];
          if (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n4]>=0 || tag[n5]>=0) // Test corner volume if 4 nodes marching
          {
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            p5 = mesh1->nodep[n5];
            if ((vol=prism_volume(p0,p1,p2,p3,p4,p5)) < 0.0)
            {
              n++;
              if (vol < negmax)
              {
                negmax = vol;
                nmax = c;
              }
            } else
            {
              double vol_tol = vol*0.05;
              vol_tol = 0.0;

              p6 = p0 = mesh1->nodep[n0];
              p7 = p1 = mesh1->nodep[n1];
              p8 = p2 = mesh1->nodep[n2];
              p9 = p3 = mesh1->nodep[n3];
              p10= p4 = mesh1->nodep[n4];
              p11= p5 = mesh1->nodep[n5];
              tfac = 1.0;
              geo = geo1;
              for (l=0; l < lay; l++)
              {
                if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2] && l >= tag[n3] && l >= tag[n4] && l >= tag[n5])
                  break;

                if (v_layers < 0)
                {
                  p0 = p6;
                  p1 = p7;
                  p2 = p8;
                  p3 = p9;
                  p4 = p10;
                  p5 = p11;
                }

                if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                {
                  vec = nrm[n0];
                  p6 = p0 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n0];
                }
                if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                {
                  vec = nrm[n1];
                  p7 = p1 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n1];
                }
                if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                {
                  vec = nrm[n2];
                  p8 = p2 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n2];
                }
                if ((l < tag[n3] && v_layers < 0) || (v_layers > 0 && l == tag[n3]-1 && l == lay-1))
                {
                  vec = nrm[n3];
                  p9 = p3 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n3];
                }
                if ((l < tag[n4] && v_layers < 0) || (v_layers > 0 && l == tag[n4]-1 && l == lay-1))
                {
                  vec = nrm[n4];
                  p10 = p4 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n4];
                }
                if ((l < tag[n5] && v_layers < 0) || (v_layers > 0 && l == tag[n5]-1 && l == lay-1))
                {
                  vec = nrm[n5];
                  p11 = p5 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n5];
                }
                if (v_layers < 0 || l == lay-1)
                {
                  if (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0 && tag[n5]>=0)
                  {
                    if ((vol=prism_volume(p6,p7,p8,p9,p10,p11)) < vol_tol)
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                  } else
                  {
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0)) &&
                        ((vol=tetrahedral_volume(p6,p7,p8,p9)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n4]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p7,p8,p6,p10)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n5]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n5]>=0)) &&
                        ((vol=tetrahedral_volume(p8,p6,p7,p11)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n3]>=0 || tag[n4]>=0 || tag[n5]>=0)) ||
                         (tag[n0]>=0 && tag[n3]>=0 && tag[n4]>=0 && tag[n5]>=0)) &&
                        ((vol=tetrahedral_volume(p9,p11,p10,p6)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n1]>=0 || tag[n3]>=0 || tag[n4]>=0 || tag[n5]>=0)) ||
                         (tag[n1]>=0 && tag[n3]>=0 && tag[n4]>=0 && tag[n5]>=0)) &&
                        ((vol=tetrahedral_volume(p10,p9,p11,p7)) < 0.0))
                    {
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n2]>=0 || tag[n3]>=0 || tag[n4]>=0 || tag[n5]>=0)) ||
                         (tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0 && tag[n5]>=0)) &&
                        ((vol=tetrahedral_volume(p11,p10,p9,p8)) < 0.0))
                    {
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                  }
                }
                tfac *= geo;
                geo *= geo2;
                geo = MAX(geo_min,MIN(geo,geo_max));
              }
            }
          }
        }
        gdouble = negmax;
#ifdef PARALLEL
        MPI_Allreduce(&negmax,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);         
#endif
        if (gdouble < 0.0)
        {
          cg = Point(0.0,0.0,0.0);
          if (nmax >= 0)
          {
            n0 = mesh1->pri_n[nmax][0];
            n1 = mesh1->pri_n[nmax][1];
            n2 = mesh1->pri_n[nmax][2];
            n3 = mesh1->pri_n[nmax][3];
            n4 = mesh1->pri_n[nmax][4];
            n5 = mesh1->pri_n[nmax][5];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            p5 = mesh1->nodep[n5];
            cg = (p0+p1+p2+p3+p4+p5)/6.0;
          }
#ifdef PARALLEL
          cmx_g_in.val = negmax;
          cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&(cg[0]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[1]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[2]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&cg[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif
          if (my_rank == 0)
          {
            fprintf(out_f,"\n  Minimum prism volume = %lg @ (%lg, %lg, %lg)",gdouble,cg[0],cg[1],cg[2]);
            fflush(out_f);
          }
        }
        
        negmax = 1.0e20;
        for (c=0; c < mesh1->nhex && vtest; c++)
        {
          n0 = mesh1->hex_n[c][0];
          n1 = mesh1->hex_n[c][1];
          n2 = mesh1->hex_n[c][2];
          n3 = mesh1->hex_n[c][3];
          n4 = mesh1->hex_n[c][4];
          n5 = mesh1->hex_n[c][5];
          n6 = mesh1->hex_n[c][6];
          n7 = mesh1->hex_n[c][7];
          if (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0 ||
              tag[n4]>=0 || tag[n5]>=0 || tag[n6]>=0 || tag[n7]>=0) // Test corner volume if 4 nodes marching
          {
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            p5 = mesh1->nodep[n5];
            p6 = mesh1->nodep[n6];
            p7 = mesh1->nodep[n7];
            if ((vol=hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7)) < 0.0)
            {
              n++;
              if (vol < negmax)
              {
                negmax = vol;
                nmax = c;
              }
            } else
            {
              double vol_tol = vol*0.05;
              vol_tol = 0.0;

              p8 = p0 = mesh1->nodep[n0];
              p9 = p1 = mesh1->nodep[n1];
              p10= p2 = mesh1->nodep[n2];
              p11= p3 = mesh1->nodep[n3];
              p12= p4 = mesh1->nodep[n4];
              p13= p5 = mesh1->nodep[n5];
              p14= p6 = mesh1->nodep[n6];
              p15= p7 = mesh1->nodep[n7];
              tfac = 1.0;
              geo = geo1;
              for (l=0; l < lay; l++)
              {
                if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2] && l >= tag[n3] &&
                    l >= tag[n4] && l >= tag[n5] && l >= tag[n6] && l >= tag[n7])
                  break;

                if (v_layers < 0)
                {
                  p0 = p8;
                  p1 = p9;
                  p2 = p10;
                  p3 = p11;
                  p4 = p12;
                  p5 = p13;
                  p6 = p14;
                  p7 = p15;
                }

                if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                {
                  vec = nrm[n0];
                  p8 = p0 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n0];
                }
                if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                {
                  vec = nrm[n1];
                  p9 = p1 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n1];
                }
                if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                {
                  vec = nrm[n2];
                  p10 = p2 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n2];
                }
                if ((l < tag[n3] && v_layers < 0) || (v_layers > 0 && l == tag[n3]-1 && l == lay-1))
                {
                  vec = nrm[n3];
                  p11 = p3 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n3];
                }
                if ((l < tag[n4] && v_layers < 0) || (v_layers > 0 && l == tag[n4]-1 && l == lay-1))
                {
                  vec = nrm[n4];
                  p12 = p4 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n4];
                }
                if ((l < tag[n5] && v_layers < 0) || (v_layers > 0 && l == tag[n5]-1 && l == lay-1))
                {
                  vec = nrm[n5];
                  p13 = p5 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n5];
                }
                if ((l < tag[n6] && v_layers < 0) || (v_layers > 0 && l == tag[n6]-1 && l == lay-1))
                {
                  vec = nrm[n6];
                  p14 = p6 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n6];
                }
                if ((l < tag[n7] && v_layers < 0) || (v_layers > 0 && l == tag[n7]-1 && l == lay-1))
                {
                  vec = nrm[n7];
                  p15 = p7 + Point(vec[0], vec[1], vec[2])*tfac/aao_fac*desired[n7];
                }
                if (v_layers < 0 || l == lay-1)
                {
                  if (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n4]>=0 && tag[n5]>=0 && tag[n6]>=0 && tag[n7]>=0)
                  {
                    if ((vol=hexahedral_volume(p8,p9,p10,p11,p12,p13,p14,p15)) < vol_tol)
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                      if (tag[n6] > l)
                      {
                        tag[n6] = l;
                        i++;
                      }
                      if (tag[n7] > l)
                      {
                        tag[n7] = l;
                        i++;
                      }
                    }
                  } else
                  {
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n3]>=0 || tag[n4]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n3]>=0 && tag[n4]>=0)) &&
                        ((vol=tetrahedral_volume(p8,p9,p11,p12)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n5]>=0)) ||
                         (tag[n0]>=0 && tag[n1]>=0 && tag[n2]>=0 && tag[n5]>=0)) &&
                        ((vol=tetrahedral_volume(p9,p10,p8,p13)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n6]>=0)) ||
                         (tag[n1]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n6]>=0)) &&
                        ((vol=tetrahedral_volume(p10,p11,p9,p14)) < 0.0))
                    {
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n6] > l)
                      {
                        tag[n6] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n2]>=0 || tag[n3]>=0 || tag[n7]>=0)) ||
                         (tag[n0]>=0 && tag[n2]>=0 && tag[n3]>=0 && tag[n7]>=0)) &&
                        ((vol=tetrahedral_volume(p11,p8,p10,p15)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n7] > l)
                      {
                        tag[n7] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n0]>=0 || tag[n4]>=0 || tag[n5]>=0 || tag[n7]>=0)) ||
                         (tag[n0]>=0 && tag[n4]>=0 && tag[n5]>=0 && tag[n7]>=0)) &&
                        ((vol=tetrahedral_volume(p12,p15,p13,p8)) < 0.0))
                    {
                      if (tag[n0] > l)
                      {
                        tag[n0] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                      if (tag[n7] > l)
                      {
                        tag[n7] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n1]>=0 || tag[n4]>=0 || tag[n5]>=0 || tag[n6]>=0)) ||
                         (tag[n1]>=0 && tag[n4]>=0 && tag[n5]>=0 && tag[n6]>=0)) &&
                        ((vol=tetrahedral_volume(p13,p12,p14,p9)) < 0.0))
                    {
                      if (tag[n1] > l)
                      {
                        tag[n1] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                      if (tag[n6] > l)
                      {
                        tag[n6] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n2]>=0 || tag[n5]>=0 || tag[n6]>=0 || tag[n7]>=0)) ||
                         (tag[n2]>=0 && tag[n5]>=0 && tag[n6]>=0 && tag[n7]>=0)) &&
                        ((vol=tetrahedral_volume(p14,p13,p15,p10)) < 0.0))
                    {
                      if (tag[n2] > l)
                      {
                        tag[n2] = l;
                        i++;
                      }
                      if (tag[n5] > l)
                      {
                        tag[n5] = l;
                        i++;
                      }
                      if (tag[n6] > l)
                      {
                        tag[n6] = l;
                        i++;
                      }
                      if (tag[n7] > l)
                      {
                        tag[n7] = l;
                        i++;
                      }
                    }
                    if (((v_layers < 0 && (tag[n3]>=0 || tag[n4]>=0 || tag[n6]>=0 || tag[n7]>=0)) ||
                         (tag[n3]>=0 && tag[n4]>=0 && tag[n6]>=0 && tag[n7]>=0)) &&
                        ((vol=tetrahedral_volume(p15,p14,p12,p11)) < 0.0))
                    {
                      if (tag[n3] > l)
                      {
                        tag[n3] = l;
                        i++;
                      }
                      if (tag[n4] > l)
                      {
                        tag[n4] = l;
                        i++;
                      }
                      if (tag[n6] > l)
                      {
                        tag[n6] = l;
                        i++;
                      }
                      if (tag[n7] > l)
                      {
                        tag[n7] = l;
                        i++;
                      }
                    }
                  }
                }
                tfac *= geo;
                geo *= geo2;
                geo = MAX(geo_min,MIN(geo,geo_max));
              }
            }
          }
        }
        gdouble = negmax;
#ifdef PARALLEL
        MPI_Allreduce(&negmax,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);         
#endif
        if (gdouble < 0.0)
        {
          cg = Point(0.0,0.0,0.0);
          if (nmax >= 0)
          {
            n0 = mesh1->hex_n[nmax][0];
            n1 = mesh1->hex_n[nmax][1];
            n2 = mesh1->hex_n[nmax][2];
            n3 = mesh1->hex_n[nmax][3];
            n4 = mesh1->hex_n[nmax][4];
            n5 = mesh1->hex_n[nmax][5];
            n6 = mesh1->hex_n[nmax][6];
            n7 = mesh1->hex_n[nmax][7];
            p0 = mesh1->nodep[n0];
            p1 = mesh1->nodep[n1];
            p2 = mesh1->nodep[n2];
            p3 = mesh1->nodep[n3];
            p4 = mesh1->nodep[n4];
            p5 = mesh1->nodep[n5];
            p6 = mesh1->nodep[n6];
            p7 = mesh1->nodep[n7];
            cg = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
          }
#ifdef PARALLEL
          cmx_g_in.val = negmax;
          cmx_g_in.proc = my_rank;
          MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
          if (cmx_g_out.proc != 0)
          {
            if (my_rank == cmx_g_out.proc)
            {
              ptag = my_rank;
              dest = 0;
              MPI_Send(&(cg[0]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[1]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
              MPI_Send(&(cg[2]), 1, MPI_DOUBLE, dest, ptag, MPI_COMM_WORLD);
            }
            if (my_rank == 0)
            {
              ptag = source = cmx_g_out.proc;
              MPI_Recv(&cg[0], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[1], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
              MPI_Recv(&cg[2], 1, MPI_DOUBLE, source, ptag, MPI_COMM_WORLD, &status);
            }
          }
#endif
          if (my_rank == 0)
          {
            fprintf(out_f,"\n  Minimum hexahedral volume = %lg @ (%lg, %lg, %lg)",gdouble,cg[0],cg[1],cg[2]);
            fflush(out_f);
          }
        }
        
        vtest = MIN(1,i);

        gint = n;
#ifdef PARALLEL
        MPI_Allreduce(&n,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
        if (gint > 0 && my_rank == 0)
        {
          fprintf(out_f,"\nLayer %d, # of initial negative volumes = %d",lay,gint);
          fprintf(out_f,"\n\nUNABLE TO ENSURE VALID NEW ELEMENTS");
#ifdef PARALLEL
          MPI_Abort(MPI_COMM_WORLD,0);
#endif
          exit(0);
        }

        for (b=0; b < mesh1->nb && jtest; b++)
        {
          if (geom->layers[b] == tgn)
          {
            for (c=0; c < mesh1->nt[b]; c++)
            {
              n0 = mesh1->t_n[b][c][0];
              n1 = mesh1->t_n[b][c][1];
              n2 = mesh1->t_n[b][c][2];
              if (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0)
              {
                p3 = p0 = mesh1->nodep[n0];
                p4 = p1 = mesh1->nodep[n1];
                p5 = p2 = mesh1->nodep[n2];
                tfac = 1.0;
                geo = geo1;
                for (l=0; l < lay; l++)
                {
                  if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2])
                    break;

                  if (v_layers < 0)
                  {
                    p0 = p3;
                    p1 = p4;
                    p2 = p5;
                  }

                  if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                  {
                    vec = nrm[n0];
                    p3 = p0 + Point(vec[0], vec[1], vec[2])*tfac*desired[n0];
                  }
                  if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                  {
                    vec = nrm[n1];
                    p4 = p1 + Point(vec[0], vec[1], vec[2])*tfac*desired[n1];
                  }
                  if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                  {
                    vec = nrm[n2];
                    p5 = p2 + Point(vec[0], vec[1], vec[2])*tfac*desired[n2];
                  }

                  if (v_layers < 0 || l == lay-1)
                  {
                    double vol0=tetrahedral_volume(p0,p1,p2,p3);
                    double vol1=tetrahedral_volume(p1,p2,p0,p4);
                    double vol2=tetrahedral_volume(p2,p0,p1,p5);
                    double vol3=tetrahedral_volume(p3,p5,p4,p0);
                    double vol4=tetrahedral_volume(p4,p3,p5,p1);
                    double vol5=tetrahedral_volume(p5,p4,p3,p2);
                    double vol_tol = 0.0;

                    if (distance(p0,p3) > tol && (vol0 <= vol_tol || vol3 <= vol_tol)  && tag[n0] > l)
                    {
                      tag[n0] = l;
                      j++;
                    }
                    if (distance(p1,p4) > tol && (vol1 <= vol_tol || vol4 <= vol_tol)  && tag[n1] > l)
                    {
                      tag[n1] = l;
                      j++;
                    }
                    if (distance(p2,p5) > tol && (vol2 <= vol_tol || vol5 <= vol_tol)  && tag[n2] > l)
                    {
                      tag[n2] = l;
                      j++;
                    }
                  }

                  tfac *= geo;
                  geo *= geo2;
                  geo = MAX(geo_min,MIN(geo,geo_max));
                }
              }
            }
            
            for (c=0; c < mesh1->nq[b]; c++)
            {
              n0 = mesh1->q_n[b][c][0];
              n1 = mesh1->q_n[b][c][1];
              n2 = mesh1->q_n[b][c][2];
              n3 = mesh1->q_n[b][c][3];
              if (tag[n0]>=0 || tag[n1]>=0 || tag[n2]>=0 || tag[n3]>=0)
              {
                p4 = p0 = mesh1->nodep[n0];
                p5 = p1 = mesh1->nodep[n1];
                p6 = p2 = mesh1->nodep[n2];
                p7 = p3 = mesh1->nodep[n3];
                tfac = 1.0;
                geo = geo1;
                for (l=0; l < lay; l++)
                {
                  if (l >= tag[n0] && l >= tag[n1] && l >= tag[n2] && l >= tag[n3])
                    break;

                  if (v_layers < 0)
                  {
                    p0 = p4;
                    p1 = p5;
                    p2 = p6;
                    p3 = p7;
                  }

                  if ((l < tag[n0] && v_layers < 0) || (v_layers > 0 && l == tag[n0]-1 && l == lay-1))
                  {
                    vec = nrm[n0];
                    p4 = p0 + Point(vec[0], vec[1], vec[2])*tfac*desired[n0];
                  }
                  if ((l < tag[n1] && v_layers < 0) || (v_layers > 0 && l == tag[n1]-1 && l == lay-1))
                  {
                    vec = nrm[n1];
                    p5 = p1 + Point(vec[0], vec[1], vec[2])*tfac*desired[n1];
                  }
                  if ((l < tag[n2] && v_layers < 0) || (v_layers > 0 && l == tag[n2]-1 && l == lay-1))
                  {
                    vec = nrm[n2];
                    p6 = p2 + Point(vec[0], vec[1], vec[2])*tfac*desired[n2];
                  }
                  if ((l < tag[n3] && v_layers < 0) || (v_layers > 0 && l == tag[n3]-1 && l == lay-1))
                  {
                    vec = nrm[n3];
                    p7 = p3 + Point(vec[0], vec[1], vec[2])*tfac*desired[n3];
                  }

                  if (v_layers < 0 || l == lay-1)
                  {
                    double vol0=tetrahedral_volume(p0,p1,p3,p4);
                    double vol1=tetrahedral_volume(p1,p2,p0,p5);
                    double vol2=tetrahedral_volume(p2,p3,p1,p6);
                    double vol3=tetrahedral_volume(p3,p0,p2,p7);
                    double vol4=tetrahedral_volume(p4,p7,p5,p0);
                    double vol5=tetrahedral_volume(p5,p4,p6,p1);
                    double vol6=tetrahedral_volume(p6,p5,p7,p2);
                    double vol7=tetrahedral_volume(p7,p6,p4,p3);
                    double vol_tol = 0.0;

                    if (distance(p0,p4) > tol && (vol0 <= vol_tol || vol4 <= vol_tol)  && tag[n0] > l)
                    {
                      tag[n0] = l;
                      j++;
                    }
                    if (distance(p1,p5) > tol && (vol1 <= vol_tol || vol5 <= vol_tol)  && tag[n1] > l)
                    {
                      tag[n1] = l;
                      j++;
                    }
                    if (distance(p2,p6) > tol && (vol2 <= vol_tol || vol6 <= vol_tol)  && tag[n2] > l)
                    {
                      tag[n2] = l;
                      j++;
                    }
                    if (distance(p3,p7) > tol && (vol3 <= vol_tol || vol7 <= vol_tol)  && tag[n3] > l)
                    {
                      tag[n3] = l;
                      j++;
                    }
                  }

                  tfac *= geo;
                  geo *= geo2;
                  geo = MAX(geo_min,MIN(geo,geo_max));
                }
              }
            }
          }
        }
        
        jtest = MIN(1,j);

        m = lay;
        n = 0;
        for (l=0; l < mesh1->nn; l++)
          if (tag[l] >= 0)
          {
            m = MIN(m,tag[l]);
            n = MAX(n,tag[l]);
          }
        
        globali = i;
        globalj = j;
        globalk = k;
#ifdef PARALLEL
        //continue in loop until all have finished element quality checking
        MPI_Allreduce(&i,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
        MPI_Allreduce(&j,&globalj,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
        MPI_Allreduce(&k,&globalk,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
#endif

        itotal += globali;
        jtotal += globalj;
        ktotal += globalk;
        if (my_rank == 0)
        {
          fprintf(out_f,"\nLayer check sweep %d:",sweep);
          fprintf(out_f,"\n   Number of changes due to volume check = %d",globali);
          fprintf(out_f,"\n   Number of changes due to Jacobian check = %d",globalj);
          fprintf(out_f,"\n   Number of changes due to stack height check = %d",globalk);
          fflush(out_f);
        }
      } while ((globali > 0 || globalj > 0 || globalk > 0) && sweep < MAXSWEEP);

      globali = itotal;
      globalj = jtotal;
      globalk = ktotal;
#ifdef PARALLEL
      MPI_Allreduce(&itotal,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
      MPI_Allreduce(&jtotal,&globalj,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
      MPI_Allreduce(&ktotal,&globalk,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
#endif
      if (my_rank == 0 && globali > 0)
        fprintf(out_f,"\nLayer %d, total number of nodes modified based on volume test = %d",lay,itotal);
      if (my_rank == 0 && globalj > 0)
        fprintf(out_f,"\nLayer %d, total number of nodes modified based on Jacobian test = %d",lay,jtotal);
      if (my_rank == 0 && globalk > 0)
        fprintf(out_f,"\nLayer %d, total number of nodes modified based on stack height test = %d",lay,ktotal);
      if (my_rank == 0) fflush(out_f);

      if (sweep >= MAXSWEEP && my_rank == 0)
      {
        fprintf(out_f,"\nTotal number of sweeps reached the limit of %d. Something is wrong!",MAXSWEEP);
        fflush(out_f);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }

#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
     
      int globalb = n;
#ifdef PARALLEL
      MPI_Allreduce(&n,&globalb,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
      
      mxlay = MAX(mxlay,globalb);

      globalm = m;
      globaln = n;
#ifdef PARALLEL
      MPI_Allreduce(&m,&globalm,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD); 
      MPI_Allreduce(&n,&globaln,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD); 
#endif
      if (my_rank == 0)
      {
        fprintf(out_f,"\nLayer %d, Min # layers = %d, Max # layers = %d",lay,globalm,globaln);
        fflush(out_f);
      }

      //
      // determine number of new points
      //
      for (b=m=n=0; n < mesh1->nn; n++)
      {
        if (tag[n] >= 0)
        {
          if (v_layers < 0)
            b+=tag[n];
          else if (tag[n] == lay)
            b++;

          m++;
        }
      }

      globalm = m;
      globalb = b;
#ifdef PARALLEL
      MPI_Allreduce(&m,&globalm,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
      MPI_Allreduce(&b,&globalb,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
#endif
      if (my_rank == 0)
      {
        fprintf(out_f,"\nLayer %d, number of nodes tagged for insertion = %d",lay,globalm);
        fprintf(out_f,"\nLayer %d, number of new nodes to be created    = %d",lay,globalb);
        fflush(out_f);
      }
	  
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
     
      //P_VLI:  Do global comm to assure that if anyone has inserted points, we go through this loop!
      
      mxnew = MAX(mxnew,globalb);

      if (globalb > 0)
      {

        //P_VLI:  It is important to note here that all expansions are local, and after creating nodes and elems, 
		//we will check for proc boundary dups

        if (b > 0)
        {
          //
          // expand node tag
          //
          tag = (int*)realloc((void*)tag,(mesh1->nn+b)*sizeof(int));
          for (n=mesh1->nn; n < mesh1->nn+b; n++)
            tag[n] = -2; //set for viscous expansion

          //P_VLI:expand pmap here
          mesh1->pmap = (int**)realloc((void*)mesh1->pmap,(mesh1->nn+b)*sizeof(int*));
          for (n=mesh1->nn; n < mesh1->nn+b; n++)
            mesh1->pmap[n] = (int*)malloc(3*sizeof(int));
	  	  
          //P_VLI:note that all new nodes on this stack have same owning proc
          //still init all fields to -1
          for (n=mesh1->nn; n < mesh1->nn+b; n++)
          {
            mesh1->pmap[n][0] = -1;
            mesh1->pmap[n][1] = -1;
            mesh1->pmap[n][2] = -1;
          }
          
          //
          // expand the node arrays
          //
          mesh1->nodec = (Point*)realloc((void*)mesh1->nodec,(mesh1->nn+b)*sizeof(Point));
          mesh1->nodep = (Point*)realloc((void*)mesh1->nodep,(mesh1->nn+b)*sizeof(Point));
        }
        
        //
        // create new nodes
        //
        nn = mesh1->nn;
        List **stack;

        stack = new List*[mesh1->nn+b];

        for (n=0; n < mesh1->nn+b; n++)
        {
          stack[n] = 0;
          if (tag[n] > 0)
            stack[n] = new List();
        }

        for (n=0; n < mesh1->nn; n++)
        {
          mesh1->nodec[n] = mesh1->nodep[n];
          if (tag[n] > 0)
          {
            // create new nodes in the stack
            tfac = 1.0;
            geo = geo1;
            p0 = p1 = mesh1->nodep[n];
            for (i=0; i < tag[n]; i++)
            {
              vec = nrm[n];
              
              p1 = p0 + Point(vec[0], vec[1], vec[2])*tfac*desired[n];

              if (v_layers < 0 || i == lay-1)
              {
                mesh1->nodec[nn] = p0;
                mesh1->nodep[nn] = p0;
                stack[n]->Add_To_List(nn);
                //P_VLI:  need to add to pmap here
                //will leave global number set to -1 until we finish all new nodes
                mesh1->pmap[nn][0] = -1;
                mesh1->pmap[nn][1] = mesh1->pmap[n][1]; //proc number
                if (mesh1->pmap[nn][1] == my_rank)
                  mesh1->pmap[nn][2] = nn; //local index on owning proc
                else
                  mesh1->pmap[nn][2] = -1;
                nn++;
                p0 = p1;
              }

              tfac *= geo;
              geo *= geo2;
              geo = MAX(geo_min,MIN(geo,geo_max));

            }

            // set new boundary locations in computational mesh
            if (v_layers < 0 || tag[n] == lay)
              mesh1->nodec[n] = p1;
          }
        }
		
        //P_VLI: now that we have added all nodes, communicate to find highest global node number on each proc,
        //determine how many each will add, give them a range, and let them add the all together.  Then, communicate global
        //node numbers to ghost nodes on other procs
		
        //find highest global nn on current proc
        int localmax = -2; //set low for debug
        //need a check to see if no nodes added, and then cycle through list again
        for (n = 0; n < nn; n++)
          localmax = MAX(localmax,mesh1->pmap[n][0]);
          
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD); //to assure all vars ready before reduce
#endif
        
        int globalmax = localmax; //set low to assure actually working
#ifdef PARALLEL
        //now, comm to find globalmax on all procs
        MPI_Allreduce(&localmax,&globalmax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
        
        int *npp = 0;
        npp = (int*)calloc(num_procs,sizeof(int));
        
        int nnowned = 0;
        
        //we cannot just send b, since some nodes were created that are not unique (on other procs too)
        for (n = mesh1->nn; n < nn; n++)
          if (mesh1->pmap[n][1] == my_rank)
            nnowned++;
        	
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD); //to assure all arrays ready before gather  
            
        //now, we do an all gather so procs can determine range for new node numbers
        MPI_Allgather(&nnowned,1,MPI_INT,npp,1,MPI_INT,MPI_COMM_WORLD);
#endif
        
        /*for (n = 0; n < num_procs; n++)
          {
          fprintf(out_f,"\nnpp[%d] = %d\n",n,npp[n]);
          fflush(out_f);
          }*/
        
        //need to increment globalmax by one to begin at next node, which will ripple through all procs
        //this gives us true next node, or boundary for nodes if no new created
        globalmax++;
        
        //now, we can easily create global node numbers for owned nodes
        //get proc starting number
        for (n = 0; n < my_rank; n++)
          globalmax += npp[n];
          
        //will free npp after debug check
          
        for (n = mesh1->nn; n < nn; n++)
        {
          if (mesh1->pmap[n][1] == my_rank)
          { 
            mesh1->pmap[n][0] = globalmax;
        	globalmax++; //will actually leave us with one greater than max node
          }
        } 
        
        z = globalmax;
#ifdef PARALLEL
        MPI_Allreduce(&globalmax,&z,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
        
        globalmax = z; //for debug
        
        if (npp > 0 && num_procs > 0)    
          free(npp);
		
#ifdef PARALLEL
        //finally, get the global/local number for nodes owned by other procs
        sendcnt = new int [num_procs];
        recvcnt = new int [num_procs];
        for (n = 0; n < num_procs; n++)
        {
          sendcnt[n] = 0;
          recvcnt[n] = 0;
        }
        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];

        //cycle thru nodes	
        for (n = 0; n < mesh1->nn; n++)
        {
          if (tag[n] > 0 && stack[n]->max > 0 && mesh1->pmap[n][1] != my_rank)
          {
            sendcnt[mesh1->pmap[n][1]] += stack[n]->max; 
          }
        }

        //now, send and recv count for all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
		  
          MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
		  
          MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }

        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);
		
        delete [] srequest;
        delete [] rrequest;
        delete [] statuses;
	  
        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];
		
        //allocate buffers to send/recv node indices
        bdim2 = new int[num_procs];
        for (n = 0; n < num_procs; n++)
          bdim2[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //indices
        sbuff = new char*[num_procs];
        rbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
        {
          sbuff[n] = new char[bdim2[n]];
          rbuff[n] = new char[bdim2[n]];
        }
		
        //now, pack node indices as per procs that will be sending them
        for (p = 0; p < num_procs; p++)
        {
          if (sendcnt[p] == 0 || p == my_rank)
            continue;

          //set position to 0
          sposition=0;

          //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
          for (n = 0; n < mesh1->nn; n++)
          {
            if (mesh1->pmap[n][1] == p && tag[n] > 0 && stack[n]->max > 0)
            {
              for (q = 0; q < stack[n]->max; q++)
              {
                MPI_Pack(&(mesh1->pmap[n][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(stack[n]->list[q]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
        }

        //now, send and recv packets from all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (sendcnt[p] > 0 && p != my_rank) 
          {
            MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }

          if (recvcnt[p] > 0 && p != my_rank) 
          {
            MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }

        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        delete [] srequest;
        delete [] rrequest;
        delete [] statuses;

        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];
	  
        //now, delete sbuff and resize
        for (n = 0; n < num_procs; n++)
          if (sbuff[n] != 0) delete [] sbuff[n];
        delete [] sbuff;
	  
        //resize
        bdim21 = new int[num_procs];
        for (n = 0; n < num_procs; n++)
          bdim21[n] = MAX(sendcnt[n],recvcnt[n])*3*sizeof(int); //nodes and indices
        sbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
          sbuff[n] = new char[bdim21[n]];
	  
        temp_node = 0; //to hold recv proc local index
        temp_node2 = 0; //to hold send proc local index
	  
        //now, unpack local indices and repack nodes with other index
        for (p = 0; p < num_procs; p++)
        {
          sposition = rposition = 0; //reset pos

          if (recvcnt[p] == 0 || p == my_rank)
            continue; //do not want it to go into do while loop unnecessarily

          n = 0; //reset

          do
          {
            //unpack index initially
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

            if (stack[temp_node] == 0)
            {
              //unpack index to keep in line
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              fprintf(stderr,"\n WARNING:  Node L = %d/ L on proc %d = %d has been expanded on ghost but not on real proc. \n",temp_node,p,temp_node2);
              fflush(stderr);
              continue;
            }

            //begin pack, unpack
            for (q = 0; q < stack[temp_node]->max; q++)
            {
              if (q > 0)
              {
                //unpack index uselessly
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
              }
              //unpack index to send back
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //pack node and index
              MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(mesh1->pmap[stack[temp_node]->list[q]][0]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(mesh1->pmap[stack[temp_node]->list[q]][2]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              n++; //help recvcnt along
            }
          } while (n < recvcnt[p]);
        }

        //now, delete rbuff and bdim2, resize
        delete [] bdim2;
        for (n = 0; n < num_procs; n++)
          if (rbuff[n] != 0) delete [] rbuff[n];
        delete [] rbuff;
 
        //resize
        rbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
          rbuff[n] = new char[bdim21[n]];
	  	
        //now, send and recv packets from all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (recvcnt[p] > 0 && p != my_rank) 
          {
            MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }

          if (sendcnt[p] > 0 && p != my_rank) 
          {
            MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
 
        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        //finally, unpack nodes in proper position
        for (p = 0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
        
          rposition = 0; //reset pos
          for (n = 0; n < sendcnt[p]; n++)
          {
            //unpack index and node
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->pmap[temp_node][0]),1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->pmap[temp_node][2]),1,MPI_INT,MPI_COMM_WORLD);
          }
        }	
 
        //debug
        for (n = 0; n < nn; n++)
        {
          if (mesh1->pmap[n][0] < 0 || mesh1->pmap[n][0] >= globalmax)
          {
            fprintf(stderr,"\nCATASTROPHIC ERROR:  GN for %d = %d, which is out of bounds!\n",n,mesh1->pmap[n][0]);
            fflush(stderr);
          }
          if (mesh1->pmap[n][1] < 0 || mesh1->pmap[n][1] >= num_procs)
          {
            fprintf(stderr,"\nCATASTROPHIC ERROR:  Proc number for %d = %d, which is out of bounds!\n",n,mesh1->pmap[n][1]);
            fflush(stderr);
          }
          //no need for proc local check...only used in comm and other indices would belie it being off, plus no knowable, useful cut-off
        }

        //finally, free MPI mem
        delete[] sendcnt;
        delete[] recvcnt;
        delete[] srequest;
        delete[] rrequest;
        delete[] statuses;
        for (n = 0; n < num_procs; n++)
        {
          if (sbuff[n] != 0) delete [] sbuff[n];
          if (rbuff[n] != 0) delete [] rbuff[n];
        }
        delete[] sbuff;
        delete[] rbuff;
        delete[] bdim21;

        //clean other mem
        delete[] nrm;
        delete[] desired;
#endif
        
        /*//debug
        for (n = 0; n < nn; n++)
          {
          i = mesh1->pmap[n][0];
          for (j = 0; j < nn; j++)
            {
            if (j == n)
              continue;
            if (mesh1->pmap[j][0] == i)
              {
              fprintf(out_f,"\nYou have duplicate global indices @ local %d & %d global %d\n",n,j,i);
              fflush(out_f);
              }
            }
          }*/

        //
        // perform Linear-Elastic smoothing
        //
        int success=0;
        int mxsm[2];
        mxsm[0]=osmoo1;
        mxsm[1]=osmoo2;
        int istart = 1;
        if (v_layers > 0)
          istart = 0;

        for (i=istart; i < 2 && !success; i++)
        {
          if (my_rank == 0) fprintf(out_f,"\nNumber of outer smoothing loops = %d",mxsm[i]);
          if (my_rank == 0) fflush(out_f);
 
          //P_VLI:  this run of smoothing only includes last layer, so all new nodes and elements can be ignored in
          //ghost calculations for nodes to ignore and then communicate new u,v,w at end

          int success1 = 0;
		  
#ifdef PARALLEL
          MPI_Barrier(MPI_COMM_WORLD); //sync
#endif
		  
          success1=mesh1->mesh_move(geom,mxsm[i],lsmoo,tgn,cflag,nsearch,geo1,geo2);

          success = success1;
#ifdef PARALLEL
          MPI_Barrier(MPI_COMM_WORLD); //sync
          MPI_Allreduce(&success1,&success,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD); //P_VLI: assure that we don't do this until all have finished trying, and one has no success
#endif

          if (!success) // reset mesh
          {
            if (osmoo2 <= osmoo1)
              break;
            if (my_rank == 0) fprintf(out_f,"\nResetting mesh for another pass.");
            if (my_rank == 0) fflush(out_f);

            for (n=0; n < mesh1->nn; n++)
            {
              mesh1->nodep[n] = mesh1->nodec[n];
              if (tag[n] > 0 && stack[n]->max > 0)
                mesh1->nodep[n] = mesh1->nodep[stack[n]->list[0]];
            }
          }
        }
        if (!success)
        {
          if (my_rank == 0) fprintf(out_f,"\nMesh smoothing failed!");
          if (my_rank == 0) fflush(out_f);
          //if (abs(v_layers) == 1)
          //{
          //  // SPECIAL CASE: replace the coordinate (P_OPT can smooth the mesh)
          //  for (n=0; n < mesh1->nn; n++)
          //    mesh1->nodep[n] = mesh1->nodec[n];
          //} else
          //{
          //  //MPI_Abort(MPI_COMM_WORLD,0);
          //  //exit(0);
          //}
          fflag = 1;
        }
#ifdef PARALLEL
        i = fflag;
        MPI_Allreduce(&i,&fflag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
    
        // redistribute the stack in case of tangent node movement due to reprojection
        if (v_layers < 0)
        {
          for (n=0; n < mesh1->nn; n++)
          {
            if (tag[n] > 0 && mesh1->pmap[n][1] == my_rank)
            {
              m = stack[n]->list[0];
              fac = distance(mesh1->nodep[m],mesh1->nodep[n])/distance(mesh1->nodec[m],mesh1->nodec[n]);
              p0 = mesh1->nodep[m];
              for (i=1; i < stack[n]->max-1; i++)
              {
                int jm1 = stack[n]->list[i-1];
                j = stack[n]->list[i];
                mesh1->nodep[j] = p0 + (mesh1->nodec[j]-mesh1->nodec[jm1])*fac;
                p0 = mesh1->nodep[j];
              }
            }
          }
        }

        mesh1->nn = nn;

        //
        // count the number of new elements
        //
        if (my_rank == 0) fprintf(out_f,"\nCounting number of new elements.");
        if (my_rank == 0) fflush(out_f);
        ntet = mesh1->ntet;
        npyr = mesh1->npyr;
        npri = mesh1->npri;
        nhex = mesh1->nhex;
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] != tgn)
            continue;
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            for (k=0; k < lay; k++)
            {
              j = 0;
              if (tag[n0] > 0 && k < stack[n0]->max) j++;
              if (tag[n1] > 0 && k < stack[n1]->max) j++;
              if (tag[n2] > 0 && k < stack[n2]->max) j++;
              if (j==3)
                npri++;
              if (j==2)
                npyr++;
              if (j==1)
                ntet++;
              if (j==0)
                break;
            }
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            for (k=0; k < lay; k++)
            {
              j = 0;
              if (tag[n0] > 0 && k < stack[n0]->max) j++;
              if (tag[n1] > 0 && k < stack[n1]->max) j++;
              if (tag[n2] > 0 && k < stack[n2]->max) j++;
              if (tag[n3] > 0 && k < stack[n3]->max) j++;
              if (j==4)
                nhex++;
              if (j == 3)
              {
                fprintf(stderr,"\nQuad with three nodes marked! Not allowed!");
                fflush(stderr);
#ifdef PARALLEL
                MPI_Abort(MPI_COMM_WORLD,0);
#endif
                exit(0);
              }
              if (j == 2)
              {
                if ((tag[n0] >= k+1 && tag[n1] >= k+1) ||
                    (tag[n1] >= k+1 && tag[n2] >= k+1) ||
                    (tag[n2] >= k+1 && tag[n3] >= k+1) ||
                    (tag[n3] >= k+1 && tag[n0] >= k+1))
                  npri++;
                else
                {
                  fprintf(stderr,"\nQuad with two opposite nodes marked! Not allowed!");
                  fprintf(stderr,"\n   k = %d",k);
                  fprintf(stderr,"\n   tag[n0] = %d",tag[n0]);
                  fprintf(stderr,"\n   tag[n1] = %d",tag[n1]);
                  fprintf(stderr,"\n   tag[n2] = %d",tag[n2]);
                  fprintf(stderr,"\n   tag[n3] = %d\n",tag[n3]);
                  fflush(stderr);
#ifdef PARALLEL
                  MPI_Abort(MPI_COMM_WORLD,0);
#endif
                  exit(0);
                }
              }
              if (j == 1)
              {
                fprintf(stderr,"\nQuad with one node marked! Not allowed!");
                fflush(stderr);
#ifdef PARALLEL
                MPI_Abort(MPI_COMM_WORLD,0);
#endif
                exit(0);
              }
              if (j==0)
                break;
            }
          }
        }

        //P_VLI: realloc for element c_n and maps...we will init maps to -(num_procs+5) since we need neg below num_procs (and num_procs for proc 0) later

        gint = i = ntet-mesh1->ntet;
#ifdef PARALLEL
        MPI_Allreduce(&i,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
        if (gint > 0 && my_rank == 0)
        {
          fprintf(out_f,"\nNumber of new tetrahedra = %d",gint);
          fflush(out_f);
        }
        if (ntet > mesh1->ntet)
        {
          mesh1->tet_n = (int**)realloc((void*)mesh1->tet_n,ntet*sizeof(int*));
          if (num_procs > 1) 
		    mesh1->tet_map = (int*)realloc((void*)mesh1->tet_map,ntet*sizeof(int));
          for (c=mesh1->ntet; c < ntet; c++)
          {
            if (num_procs > 1)
		      mesh1->tet_map[c] = -(num_procs+5);
            mesh1->tet_n[c] = (int*)malloc(4*sizeof(int));
            for (i=0; i < 4; i++)
              mesh1->tet_n[c][i] = -1;
          }
        }
        gint = i = npyr-mesh1->npyr;
#ifdef PARALLEL
        MPI_Allreduce(&i,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
        if (gint > 0 && my_rank == 0)
        {
          fprintf(out_f,"\nNumber of new pyramid = %d",gint);
          fflush(out_f);
        }
        if (npyr > mesh1->npyr)
        {
          mesh1->pyr_n = (int**)realloc((void*)mesh1->pyr_n,npyr*sizeof(int*));
          if (num_procs > 1)
		    mesh1->pyr_map = (int*)realloc((void*)mesh1->pyr_map,npyr*sizeof(int));
          for (c=mesh1->npyr; c < npyr; c++)
          {
            if (num_procs > 1)
		      mesh1->pyr_map[c] = -(num_procs+5);
            mesh1->pyr_n[c] = (int*)malloc(5*sizeof(int));
            for (i=0; i < 5; i++)
              mesh1->pyr_n[c][i] = -1;
          }
        }
        gint = i = npri-mesh1->npri;
#ifdef PARALLEL
        MPI_Allreduce(&i,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
        if (gint > 0 && my_rank == 0)
        {
          fprintf(out_f,"\nNumber of new prisms = %d",gint);
          fflush(out_f);
        }
        if (npri > mesh1->npri)
        {
          mesh1->pri_n = (int**)realloc((void*)mesh1->pri_n,npri*sizeof(int*));
          if (num_procs > 1)
		    mesh1->pri_map = (int*)realloc((void*)mesh1->pri_map,npri*sizeof(int));
          for (c=mesh1->npri; c < npri; c++)
          {
            if (num_procs > 1)
		      mesh1->pri_map[c] = -(num_procs+5);
            mesh1->pri_n[c] = (int*)malloc(6*sizeof(int));
            for (i=0; i < 6; i++)
              mesh1->pri_n[c][i] = -1;
          }
        }
        gint = i = nhex-mesh1->nhex;
#ifdef PARALLEL
        MPI_Allreduce(&i,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
        if (gint > 0 && my_rank == 0)
        {
          fprintf(out_f,"\nNumber of new hexahedra = %d",gint);
          fflush(out_f);
        }
        if (nhex > mesh1->nhex)
        {
          mesh1->hex_n = (int**)realloc((void*)mesh1->hex_n,nhex*sizeof(int*));
          if (num_procs > 1)
		    mesh1->hex_map = (int*)realloc((void*)mesh1->hex_map,nhex*sizeof(int));
          for (c=mesh1->nhex; c < nhex; c++)
          {
            if (num_procs > 1)
		      mesh1->hex_map[c] = -(num_procs+5);
            mesh1->hex_n[c] = (int*)malloc(8*sizeof(int));
            for (i=0; i < 8; i++)
              mesh1->hex_n[c][i] = -1;
          }
        }
        
        //P_VLI:  now we communicate to find highest global element number of each type on each proc,
        //determine how many each will add (letting lowest proc with any nodes in fight set element number), give them a range, 
        //and let them add them all together.  Then, communicate global element numbers to ghost nodes on other procs using stack info

        //
        // create new elements
        //
        if (my_rank == 0) fprintf(out_f,"\nCreating new elements.");
        if (my_rank == 0) fflush(out_f);
        ntet = mesh1->ntet;
        npyr = mesh1->npyr;
        npri = mesh1->npri;
        nhex = mesh1->nhex;
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] != tgn)
            continue;
          for (i=0; i < mesh1->nt[b]; i++)
          {
            b0 = mesh1->t_n[b][i][0];
            b1 = mesh1->t_n[b][i][1];
            b2 = mesh1->t_n[b][i][2];
            for (k=0; k < lay; k++)
            {
              n0 = n3 = b0;
              n1 = n4 = b1;
              n2 = n5 = b2;

              if (tag[b0] > 0 && k < stack[b0]->max)
                n0 = stack[b0]->list[k];
              if (tag[b0] > 0 && k+1 < stack[b0]->max)
                n3 = stack[b0]->list[k+1];

              if (tag[b1] > 0 && k < stack[b1]->max)
                n1 = stack[b1]->list[k];
              if (tag[b1] > 0 && k+1 < stack[b1]->max)
                n4 = stack[b1]->list[k+1];

              if (tag[b2] > 0 && k < stack[b2]->max)
                n2 = stack[b2]->list[k];
              if (tag[b2] > 0 && k+1 < stack[b2]->max)
                n5 = stack[b2]->list[k+1];

              if (n0 == n3 && n1 == n4 && n2 == n5)
                break;

              j = 0;
              if (tag[b0] > 0 && k < stack[b0]->max) j++;
              if (tag[b1] > 0 && k < stack[b1]->max) j++;
              if (tag[b2] > 0 && k < stack[b2]->max) j++;
              if (j==3)
              {
                mesh1->pri_n[npri][0] = n0;
                mesh1->pri_n[npri][1] = n1;
                mesh1->pri_n[npri][2] = n2;
                mesh1->pri_n[npri][3] = n3;
                mesh1->pri_n[npri][4] = n4;
                mesh1->pri_n[npri][5] = n5;
                npri++;
                continue;
              }
              if (j==2)
              {
                if (n3 == n0)
                {
                  mesh1->pyr_n[npyr][0] = n4;
                  mesh1->pyr_n[npyr][1] = n5;
                  mesh1->pyr_n[npyr][2] = n2;
                  mesh1->pyr_n[npyr][3] = n1;
                  mesh1->pyr_n[npyr][4] = n0;
                  npyr++;
                  continue;
                }
                if (n4 == n1)
                {
                  mesh1->pyr_n[npyr][0] = n5;
                  mesh1->pyr_n[npyr][1] = n3;
                  mesh1->pyr_n[npyr][2] = n0;
                  mesh1->pyr_n[npyr][3] = n2;
                  mesh1->pyr_n[npyr][4] = n1;
                  npyr++;
                  continue;
                }
                if (n5 == n2)
                {
                  mesh1->pyr_n[npyr][0] = n3;
                  mesh1->pyr_n[npyr][1] = n4;
                  mesh1->pyr_n[npyr][2] = n1;
                  mesh1->pyr_n[npyr][3] = n0;
                  mesh1->pyr_n[npyr][4] = n2;
                  npyr++;
                  continue;
                }
              }
              if (j==1)
              {
                if (n3 != n0)
                {
                  mesh1->tet_n[ntet][0] = n0;
                  mesh1->tet_n[ntet][1] = n1;
                  mesh1->tet_n[ntet][2] = n2;
                  mesh1->tet_n[ntet][3] = n3;
                  ntet++;
                  continue;
                }
                if (n4 != n1)
                {
                  mesh1->tet_n[ntet][0] = n0;
                  mesh1->tet_n[ntet][1] = n1;
                  mesh1->tet_n[ntet][2] = n2;
                  mesh1->tet_n[ntet][3] = n4;
                  ntet++;
                  continue;
                }
                if (n5 != n2)
                {
                  mesh1->tet_n[ntet][0] = n0;
                  mesh1->tet_n[ntet][1] = n1;
                  mesh1->tet_n[ntet][2] = n2;
                  mesh1->tet_n[ntet][3] = n5;
                  ntet++;
                  continue;
                }
              }
            }
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            b0 = mesh1->q_n[b][i][0];
            b1 = mesh1->q_n[b][i][1];
            b2 = mesh1->q_n[b][i][2];
            b3 = mesh1->q_n[b][i][3];
            for (k=0; k < lay; k++)
            {
              n0 = n4 = b0;
              n1 = n5 = b1;
              n2 = n6 = b2;
              n3 = n7 = b3;

              if (tag[b0] > 0 && k < stack[b0]->max)
                n0 = stack[b0]->list[k];
              if (tag[b0] > 0 && k+1 < stack[b0]->max)
                n4 = stack[b0]->list[k+1];

              if (tag[b1] > 0 && k < stack[b1]->max)
                n1 = stack[b1]->list[k];
              if (tag[b1] > 0 && k+1 < stack[b1]->max)
                n5 = stack[b1]->list[k+1];

              if (tag[b2] > 0 && k < stack[b2]->max)
                n2 = stack[b2]->list[k];
              if (tag[b2] > 0 && k+1 < stack[b2]->max)
                n6 = stack[b2]->list[k+1];

              if (tag[b3] > 0 && k < stack[b3]->max)
                n3 = stack[b3]->list[k];
              if (tag[b3] > 0 && k+1 < stack[b3]->max)
                n7 = stack[b3]->list[k+1];

              if (n0 == n4 && n1 == n5 && n2 == n6 && n3 == n7)
                break;

              j = 0;
              if (tag[b0] > 0 && k < stack[b0]->max) j++;
              if (tag[b1] > 0 && k < stack[b1]->max) j++;
              if (tag[b2] > 0 && k < stack[b2]->max) j++;
              if (tag[b3] > 0 && k < stack[b3]->max) j++;
              if (j==4)
              {
                mesh1->hex_n[nhex][0] = n0;
                mesh1->hex_n[nhex][1] = n1;
                mesh1->hex_n[nhex][2] = n2;
                mesh1->hex_n[nhex][3] = n3;
                mesh1->hex_n[nhex][4] = n4;
                mesh1->hex_n[nhex][5] = n5;
                mesh1->hex_n[nhex][6] = n6;
                mesh1->hex_n[nhex][7] = n7;
                nhex++;
                continue;
              }
              if (j==2)
              {
                if (n2 == n6 && n3 == n7)
                {
                  mesh1->pri_n[npri][0] = n5;
                  mesh1->pri_n[npri][1] = n2;
                  mesh1->pri_n[npri][2] = n1;
                  mesh1->pri_n[npri][3] = n4;
                  mesh1->pri_n[npri][4] = n3;
                  mesh1->pri_n[npri][5] = n0;
                  npri++;
                  continue;
                }
                if (n3 == n7 && n0 == n4)
                {
                  mesh1->pri_n[npri][0] = n6;
                  mesh1->pri_n[npri][1] = n3;
                  mesh1->pri_n[npri][2] = n2;
                  mesh1->pri_n[npri][3] = n5;
                  mesh1->pri_n[npri][4] = n0;
                  mesh1->pri_n[npri][5] = n1;
                  npri++;
                  continue;
                }
                if (n0 == n4 && n1 == n5)
                {
                  mesh1->pri_n[npri][0] = n7;
                  mesh1->pri_n[npri][1] = n0;
                  mesh1->pri_n[npri][2] = n3;
                  mesh1->pri_n[npri][3] = n6;
                  mesh1->pri_n[npri][4] = n1;
                  mesh1->pri_n[npri][5] = n2;
                  npri++;
                  continue;
                }
                if (n1 == n5 && n2 == n6)
                {
                  mesh1->pri_n[npri][0] = n4;
                  mesh1->pri_n[npri][1] = n1;
                  mesh1->pri_n[npri][2] = n0;
                  mesh1->pri_n[npri][3] = n7;
                  mesh1->pri_n[npri][4] = n2;
                  mesh1->pri_n[npri][5] = n3;
                  npri++;
                  continue;
                }
              }
            }
          }
        }

        //P_VLI: need temp holding spot for nt, nq
        int *qn, *tn;
        qn = (int*)calloc(mesh1->nb,sizeof(int));
        tn = (int*)calloc(mesh1->nb,sizeof(int));
        
        //P_VLI:  set as initial numbers of quads, tris
        for (n = 0; n < mesh1->nb; n++)
        {
          qn[n] = mesh1->nq[n];
          tn[n] = mesh1->nt[n];
        }
        
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        // P_VLI: correct boundary definitions was moved after new element creation, since we need to add to tri_map and quad_map
        if (my_rank == 0) fprintf(out_f,"\nCorrecting boundary definitions.");
        if (my_rank == 0) fflush(out_f);
        for (b=0; b < mesh1->nb; b++)
        {
          if (geom->layers[b] == tgn)
          {
            for (i=0; i < mesh1->nt[b]; i++)
            {
              b0 = mesh1->t_n[b][i][0];
              b1 = mesh1->t_n[b][i][1];
              b2 = mesh1->t_n[b][i][2];
              if (tag[b0] > 0 && stack[b0]->max > 0)
                mesh1->t_n[b][i][0] = stack[b0]->list[0];
              if (tag[b1] > 0 && stack[b1]->max > 0)
                mesh1->t_n[b][i][1] = stack[b1]->list[0];
              if (tag[b2] > 0 && stack[b2]->max > 0)
                mesh1->t_n[b][i][2] = stack[b2]->list[0];
            }
            for (i=0; i < mesh1->nq[b]; i++)
            {
              b0 = mesh1->q_n[b][i][0];
              b1 = mesh1->q_n[b][i][1];
              b2 = mesh1->q_n[b][i][2];
              b3 = mesh1->q_n[b][i][3];
              if (tag[b0] > 0 && stack[b0]->max > 0)
                mesh1->q_n[b][i][0] = stack[b0]->list[0];
              if (tag[b1] > 0 && stack[b1]->max > 0)
                mesh1->q_n[b][i][1] = stack[b1]->list[0];
              if (tag[b2] > 0 && stack[b2]->max > 0)
                mesh1->q_n[b][i][2] = stack[b2]->list[0];
              if (tag[b3] > 0 && stack[b3]->max > 0)
                mesh1->q_n[b][i][3] = stack[b3]->list[0];
            }
          }
        }
        // create temporary count of new boundary faces
        int *ntri, *nquad;
        ntri  = new int[mesh1->nb];
        nquad = new int[mesh1->nb];
        for (b=0; b < mesh1->nb; b++)
        {
          ntri[b] = mesh1->nt[b];
          nquad[b] = mesh1->nq[b];
        }
        for (e=0; e < ne; e++)
        {
          b0 = edge[e].nodes[0];
          b1 = edge[e].nodes[1];
          b  = edge[e].boundary;
          for (k=0; k < lay; k++)
          {
            n3 = n0 = b0;
            n2 = n1 = b1;
            if (tag[b0] > 0 && k < stack[b0]->max)
              n0 = stack[b0]->list[k];
            if (tag[b0] > 0 && k+1 < stack[b0]->max)
              n3 = stack[b0]->list[k+1];
            if (tag[b1] > 0 && k < stack[b1]->max)
              n1 = stack[b1]->list[k];
            if (tag[b1] > 0 && k+1 < stack[b1]->max)
              n2 = stack[b1]->list[k+1];
            if (n0 == n3 && n1 == n2)
              break;

            if (n0 != n3 && n1 != n2)
            {
              nquad[b]++;
            }
            if (n1 != n2 && n3 == n0)
            {
              ntri[b]++;
            }
            if (n2 == n1 && n3 != n0)
            {
              ntri[b]++;
            }
          }
        } 
        for (b=0; b < mesh1->nb; b++)
        {
          if (ntri[b] > mesh1->nt[b])
          {
            mesh1->t_n[b]=(int**)realloc((void*)mesh1->t_n[b],ntri[b]*sizeof(int*));
            for (i=mesh1->nt[b]; i < ntri[b]; i++)
              mesh1->t_n[b][i] = (int*)malloc(3*sizeof(int));
            if (num_procs > 1)
              mesh1->tri_map[b]=(int*)realloc((void*)mesh1->tri_map[b],ntri[b]*sizeof(int));
          }
          if (nquad[b] > mesh1->nq[b])
          {
            mesh1->q_n[b]=(int**)realloc((void*)mesh1->q_n[b],nquad[b]*sizeof(int*));
            for (i=mesh1->nq[b]; i < nquad[b]; i++)
              mesh1->q_n[b][i] = (int*)malloc(4*sizeof(int));
            if (num_procs > 1)
              mesh1->quad_map[b]=(int*)realloc((void*)mesh1->quad_map[b],nquad[b]*sizeof(int));
          }
        }
        for (e=0; e < ne; e++)
        {
          b0 = edge[e].nodes[0];
          b1 = edge[e].nodes[1];
          b  = edge[e].boundary;
          for (k=0; k < lay; k++)
          {
            n3 = n0 = b0;
            n2 = n1 = b1;
            if (tag[b0] > 0 && k < stack[b0]->max)
              n0 = stack[b0]->list[k];
            if (tag[b0] > 0 && k+1 < stack[b0]->max)
              n3 = stack[b0]->list[k+1];
            if (tag[b1] > 0 && k < stack[b1]->max)
              n1 = stack[b1]->list[k];
            if (tag[b1] > 0 && k+1 < stack[b1]->max)
              n2 = stack[b1]->list[k+1];
            if (n0 == n3 && n1 == n2)
              break;

            if (n0 != n3 && n1 != n2)
            {
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n0] = geom->closest(mesh1->nodep[n0],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n1] = geom->closest(mesh1->nodep[n1],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n2] = geom->closest(mesh1->nodep[n2],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n3] = geom->closest(mesh1->nodep[n3],n,vec);
              //mesh1->q_n[b]=(int**)realloc((void*)mesh1->q_n[b],(mesh1->nq[b]+1)*sizeof(int*));
              //mesh1->q_n[b][mesh1->nq[b]] = (int*)malloc(4*sizeof(int));
              if (num_procs > 1)
              {
              //mesh1->quad_map[b]=(int*)realloc((void*)mesh1->quad_map[b],(mesh1->nq[b]+1)*sizeof(int));
              mesh1->quad_map[b][mesh1->nq[b]] = -(num_procs+5);//init to neg to catch error
              }
              mesh1->q_n[b][mesh1->nq[b]][0] = n0;
              mesh1->q_n[b][mesh1->nq[b]][1] = n1;
              mesh1->q_n[b][mesh1->nq[b]][2] = n2;
              mesh1->q_n[b][mesh1->nq[b]][3] = n3;
              mesh1->nq[b]++;
            }
            if (n1 != n2 && n3 == n0)
            {
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n0] = geom->closest(mesh1->nodep[n0],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n1] = geom->closest(mesh1->nodep[n1],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n2] = geom->closest(mesh1->nodep[n2],n,vec);
              //mesh1->t_n[b]=(int**)realloc((void*)mesh1->t_n[b],(mesh1->nt[b]+1)*sizeof(int*));
              //mesh1->t_n[b][mesh1->nt[b]] = (int*)malloc(3*sizeof(int));
              if (num_procs > 1)
              {
              //mesh1->tri_map[b]=(int*)realloc((void*)mesh1->tri_map[b],(mesh1->nt[b]+1)*sizeof(int));
              mesh1->tri_map[b][mesh1->nt[b]] = -(num_procs+5);//init to neg to catch error
              }
              mesh1->t_n[b][mesh1->nt[b]][0] = n0;
              mesh1->t_n[b][mesh1->nt[b]][1] = n1;
              mesh1->t_n[b][mesh1->nt[b]][2] = n2;
              mesh1->nt[b]++;
            }
            if (n2 == n1 && n3 != n0)
            {
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n0] = geom->closest(mesh1->nodep[n0],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n1] = geom->closest(mesh1->nodep[n1],n,vec);
              n = b;
              vec = Vector(0.0,0.0,0.0);
              mesh1->nodep[n3] = geom->closest(mesh1->nodep[n3],n,vec);
              //mesh1->t_n[b]=(int**)realloc((void*)mesh1->t_n[b],(mesh1->nt[b]+1)*sizeof(int*));
              //mesh1->t_n[b][mesh1->nt[b]] = (int*)malloc(3*sizeof(int));
              if (num_procs > 1)
              {
              //mesh1->tri_map[b]=(int*)realloc((void*)mesh1->tri_map[b],(mesh1->nt[b]+1)*sizeof(int));
              mesh1->tri_map[b][mesh1->nt[b]] = -(num_procs+5);//init to neg to catch error
              }
              mesh1->t_n[b][mesh1->nt[b]][0] = n0;
              mesh1->t_n[b][mesh1->nt[b]][1] = n1;
              mesh1->t_n[b][mesh1->nt[b]][2] = n3;
              mesh1->nt[b]++;
            }
          }
        } 

        delete[] ntri;
        delete[] nquad;

        if (ne > 0)
        {
          delete[] edge;
          ne=0;
        }
        
        /*for (z = 0; z < mesh1->nb; z++)
          {
          fprintf(out_f,"\nb = %d, nq = %d, nt = %d\n",z,mesh1->nq[z],mesh1->nt[z]);
          fprintf(out_f,"\nmesh1->  b = %d, nq = %d, nt = %d\n",z,qn[z],tn[z]);
          fflush(out_f);
          }*/

        //P_VLI: only need the following for growing maps in parallel...
        //NOTE: instead of whichever proc has the most nodes and if there is a tie it goes to the lowest proc, the elements will be owned by the lowest proc with any node in the race
        if (num_procs > 1)
        {
        //find highest global element number of each type on current proc (and quad and tri per bd)
        int *lmax = 0;
        lmax = (int*)calloc((4 + 2*mesh1->nb),sizeof(int));
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          lmax[n] = -2; //set low for debug
          
        //no need for 0 lmax in nodes since always start mesh with some nodes
        
        for (n = 0; n < 4; n++)
        {
          switch(n)
          {
            case 0:
              for (p = 0; p < ntet; p++)
                lmax[n] = MAX(lmax[n],mesh1->tet_map[p]);
              if (mesh1->ntet < ntet && lmax[n] == -2)
                lmax[n] = -1; //so will increment and begin at 0
              break;
            case 1:
              for (p = 0; p < npri; p++)
                lmax[n] = MAX(lmax[n],mesh1->pri_map[p]);
              if (mesh1->npri < npri && lmax[n] == -2)
                lmax[n] = -1; //so will increment and begin at 0
              break;
            case 2:
              for (p = 0; p < npyr; p++)
                lmax[n] = MAX(lmax[n],mesh1->pyr_map[p]);
              if (mesh1->npyr < npyr && lmax[n] == -2)
                lmax[n] = -1; //so will increment and begin at 0
              break;
            case 3:
              for (p = 0; p < nhex; p++)
                lmax[n] = MAX(lmax[n],mesh1->hex_map[p]);
              if (mesh1->nhex < nhex && lmax[n] == -2)
                lmax[n] = -1; //so will increment and begin at 0
              break;
            default:
              break;
          }
        }
        
        //now, do for quad and tri maps
        for (n = 0; n < 2; n++)
          {
          switch(n)
            {
            case 0: 
              for (p = 0; p < mesh1->nb; p++)
                {
                for (z = 0; z < mesh1->nt[p]; z++)
                  {
                  lmax[4 + p] = MAX(lmax[4 + p],mesh1->tri_map[p][z]);
                  }
                if (tn[p] < mesh1->nt[p] && lmax[4 + p] == -2)
                  lmax[4 + p] = -1; //so will increment and begin at 0
                }
              break;
            case 1:
              for (p = 0; p < mesh1->nb; p++)
                {
                for (z = 0; z < mesh1->nq[p]; z++)
                  {
                  lmax[4 + mesh1->nb + p] = MAX(lmax[4 + mesh1->nb + p],mesh1->quad_map[p][z]);
                  }
                if (qn[p] < mesh1->nq[p] && lmax[4 + mesh1->nb + p] == -2)
                  lmax[4 + mesh1->nb + p] = -1; //so will increment and begin at 0
                }
              break;
            default:
              fprintf(stderr,"\nYou have more than quads and triangles on boundaries!  Exiting....\n");
              fflush(stderr);
#ifdef PARALLEL
              MPI_Abort(MPI_COMM_WORLD,0);
#endif
              exit(0);
            break;
            }
          }
        

        /*//debug
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
        {
          fprintf(out_f,"\n local max %d on proc %d = %d \n",n,my_rank,lmax[n]);
          fflush(out_f);
        }*/
		
        //set up for global max of each type
        glmax = (int*)calloc((4 + 2*mesh1->nb),sizeof(int));
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          glmax[n] = -3; //set low to assure actually working

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD); //to assure all arrays ready before reduce  
 
        //now, comm to find globalmax on all procs
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          MPI_Allreduce(&(lmax[n]),&(glmax[n]),1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
		
        /*//debug
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
        {
          fprintf(out_f,"\n global max %d on proc %d = %d \n",n,my_rank,glmax[n]);
          fflush(out_f);
        }*/
        
        //set up array to hold each range for each proc on each proc
        new_elem = (int**)calloc((4 + 2*mesh1->nb),sizeof(int*));
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          new_elem[n] = (int*)calloc(num_procs,sizeof(int));
 
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          for (p = 0; p < num_procs; p++)
            new_elem[n][p] = -1; //init for debug, since overwritten in gather
		
        //hold list of number of owned elements, based on lowest proc and most nodes
        nelem_owned = (int*)calloc((4 + 2*mesh1->nb),sizeof(int));

        //create counters based on nnodes per proc
        own = (int*)calloc(num_procs,sizeof(int));

        //need to sift through and determine who will name which element (as long as unique, no matter who names)
        for (n = 0; n < 4; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (p = mesh1->ntet; p < ntet; p++)
              {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 4; q++)
                  own[mesh1->pmap[mesh1->tet_n[p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->tet_map[p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->tet_map[p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  }
              }
              break;
            case 1:
              for (p = mesh1->npri; p < npri; p++)
              {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 6; q++)
                  own[mesh1->pmap[mesh1->pri_n[p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->pri_map[p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->pri_map[p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  }
              }
              break;
            case 2:
              for (p = mesh1->npyr; p < npyr; p++)
              {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 5; q++)
                  own[mesh1->pmap[mesh1->pyr_n[p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->pyr_map[p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->pyr_map[p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  }
              }
              break;
            case 3:
              for (p = mesh1->nhex; p < nhex; p++)
              {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 8; q++)
                  own[mesh1->pmap[mesh1->hex_n[p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->hex_map[p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->hex_map[p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[n]++;
                      }
                    flag = 1;
                    }
                  }
              }
              break;
            default:
              break;
          }
        }
        
        for (n = 0; n < 2; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (z = 0; z < mesh1->nb; z++)
              {
              for (p = tn[z]; p < mesh1->nt[z]; p++)
                {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 3; q++)
                  own[mesh1->pmap[mesh1->t_n[z][p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->tri_map[z][p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[4 + z]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->tri_map[z][p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[4 + z]++;
                      }
                    flag = 1;
                    }
                  }
                }
              }
              break;
            case 1:
              for (z = 0; z < mesh1->nb; z++)
              {
              for (p = qn[z]; p < mesh1->nq[z]; p++)
                {
                //reset own
                for (q = 0; q < num_procs; q++)
                  own[q] = 0;
                //set up ownership of elem nodes
                for (q = 0; q < 4; q++)
                  own[mesh1->pmap[mesh1->q_n[z][p][q]][1]]++;
                //determine elem owner
                flag = 0;
                for (q = 0; q < num_procs && !flag; q++)
                  {
                  if (own[q] > 0 && q > 0)
                    {
                    mesh1->quad_map[z][p] = -q;
                    if (q == my_rank)
                      {
                      nelem_owned[4 + mesh1->nb + z]++;
                      }
                    flag = 1;
                    }
                  else if (own[q] > 0 && q == 0)
                    {
                    mesh1->quad_map[z][p] = -num_procs;
                    if (q == my_rank)
                      {
                      nelem_owned[4 + mesh1->nb + z]++;
                      }
                    flag = 1;
                    }
                  }
                }
              }
              break;
            default:
              break;
          }
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD); //to assure all arrays ready before gather  
  
        //now, we do an all gather so procs can determine range for new elem numbers
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          MPI_Allgather(&(nelem_owned[n]),1,MPI_INT,new_elem[n],1,MPI_INT,MPI_COMM_WORLD);
#endif

        /*//debug
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          {
          for (p = 0; p < num_procs; p++)
            {
            fprintf(out_f,"new_elem[%d][%d] = %d\n",n,p,new_elem[n][p]);
            }
          fprintf(out_f,"\n");
          }
        fflush(out_f);*/
    
        //so we start at next element number, will push through all procs and types  
        for (p = 0; p < (4 + 2*mesh1->nb); p++)
          glmax[p]++; //means a -1 has no elem created ever

        //now, we can easily create global elem numbers for owned elem
        //get proc starting number
        for (n = 0; n < my_rank; n++)
          for (p = 0; p < (4 + 2*mesh1->nb); p++)
            glmax[p] += new_elem[p][n]; 
  
        //need to give owned elem global numbers
        for (n = 0; n < 4; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (p = mesh1->ntet; p < ntet; p++)
              {
                if ((mesh1->tet_map[p] == -my_rank && my_rank != 0) || (mesh1->tet_map[p] == -num_procs && my_rank == 0))
                {
                  mesh1->tet_map[p] = glmax[n];
                  glmax[n]++;
                }
              }
              break;
            case 1:
              for (p = mesh1->npri; p < npri; p++)
              {
                if ((mesh1->pri_map[p] == -my_rank && my_rank != 0) || (mesh1->pri_map[p] == -num_procs && my_rank == 0))
                {
                  mesh1->pri_map[p] = glmax[n];
                  glmax[n]++;
                }
              }
              break;
            case 2:
              for (p = mesh1->npyr; p < npyr; p++)
              {
                if ((mesh1->pyr_map[p] == -my_rank && my_rank != 0) || (mesh1->pyr_map[p] == -num_procs && my_rank == 0))
                {
                  mesh1->pyr_map[p] = glmax[n];
                  glmax[n]++;
                }
              }
              break;
            case 3:
              for (p = mesh1->nhex; p < nhex; p++)
              {
                if ((mesh1->hex_map[p] == -my_rank && my_rank != 0) || (mesh1->hex_map[p] == -num_procs && my_rank == 0))
                {
                  mesh1->hex_map[p] = glmax[n];
                  glmax[n]++;
                }
              }
              break;
            default:
              break;
          }
        }
        
        //need to give owned elem global numbers
        for (n = 0; n < 2; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (z = 0; z < mesh1->nb; z++)
              {
              for (p = tn[z]; p < mesh1->nt[z]; p++)
              {
                if ((mesh1->tri_map[z][p] == -my_rank && my_rank != 0) || (mesh1->tri_map[z][p] == -num_procs && my_rank == 0))
                {
                  mesh1->tri_map[z][p] = glmax[4 + z];
                  glmax[4 + z]++;
                }
              }
              }
              break;
            case 1:
              for (z = 0; z < mesh1->nb; z++)
              {
              for (p = qn[z]; p < mesh1->nq[z]; p++)
              {
                if ((mesh1->quad_map[z][p] == -my_rank && my_rank != 0) || (mesh1->quad_map[z][p] == -num_procs && my_rank == 0))
                {
                  mesh1->quad_map[z][p] = glmax[4 + mesh1->nb + z];
                  glmax[4 + mesh1->nb + z]++;
                }
              }
              }
              break;
            default:
              break;
          }
        }
        
#ifdef PARALLEL
        //now, comm to find globalmax on all procs for debug
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          MPI_Allreduce(&(glmax[n]),&(lmax[n]),1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
          
        //now, comm to find globalmax on all procs for debug
        for (n = 0; n < (4 + 2*mesh1->nb); n++)
          glmax[n] = lmax[n];
          
        //debug
        //for (n = 0; n < (4 + 2*mesh1->nb); n++)
          //fprintf(out_f,"\nglmax[%d] = %d, which num_elem, not index!\n",n,glmax[n]);
        //fflush(out_f);

        if (lmax > 0)
          free(lmax);
        if (new_elem > 0 && num_procs > 0)
        {
          for (n = 0; n < (4 + 2*mesh1->nb); n++)
            free(new_elem[n]);
          free(new_elem);
        }
        //will free glmax after check
        if (nelem_owned > 0)
          free(nelem_owned);
        if (own > 0 && num_procs > 0)
          free(own); 

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);

        //finally, get the global number for elem owned by other procs
        sendcnt = new int [num_procs];
        recvcnt = new int [num_procs];
        for (n = 0; n < num_procs; n++)
        {
          sendcnt[n] = 0;
          recvcnt[n] = 0;
        }
        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];

        //cycle thru elem	
        for (n = 0; n < 4; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (p = mesh1->ntet; p < ntet; p++)
                {
                if (mesh1->tet_map[p] < 0 && mesh1->tet_map[p] > -num_procs && mesh1->tet_map[p] != -my_rank)
                  sendcnt[abs(mesh1->tet_map[p])]++;
                if (mesh1->tet_map[p] == -num_procs && my_rank != 0)
                  sendcnt[0]++;
                }
              break;
            case 1:
              for (p = mesh1->npri; p < npri; p++)
                {
                if (mesh1->pri_map[p] < 0 && mesh1->pri_map[p] > -num_procs && mesh1->pri_map[p] != -my_rank)
                  sendcnt[abs(mesh1->pri_map[p])]++;
                else if (mesh1->pri_map[p] == -num_procs && my_rank != 0)
                  sendcnt[0]++;
                }
              break;
            case 2:
              for (p = mesh1->npyr; p < npyr; p++)
                if (mesh1->pyr_map[p] < 0 && mesh1->pyr_map[p] > -num_procs && mesh1->pyr_map[p] != -my_rank)
                  sendcnt[abs(mesh1->pyr_map[p])]++;
                else if (mesh1->pyr_map[p] == -num_procs && my_rank != 0)
                  sendcnt[0]++;
              break;
            case 3:
              for (p = mesh1->nhex; p < nhex; p++)
                if (mesh1->hex_map[p] < 0 && mesh1->hex_map[p] > -num_procs && mesh1->hex_map[p] != -my_rank)
                  sendcnt[abs(mesh1->hex_map[p])]++;
                else if (mesh1->hex_map[p] == -num_procs && my_rank != 0)
                  sendcnt[0]++;
              break;
            default:
              break;
          }
        }
        
        //cycle thru elem	
        for (n = 0; n < 2; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (z = 0; z < mesh1->nb; z++)
                for (p = tn[z]; p < mesh1->nt[z]; p++)
                  if (mesh1->tri_map[z][p] < 0 && mesh1->tri_map[z][p] > -num_procs && mesh1->tri_map[z][p] != -my_rank)
                    sendcnt[abs(mesh1->tri_map[z][p])]++;
                  else if (mesh1->tri_map[z][p] == -num_procs && my_rank != 0)
                    sendcnt[0]++;
              break;
            case 1:
              for (z = 0; z < mesh1->nb; z++)
                for (p = qn[z]; p < mesh1->nq[z]; p++)
                  if (mesh1->quad_map[z][p] < 0 && mesh1->quad_map[z][p] > -num_procs && mesh1->quad_map[z][p] != -my_rank)
                    sendcnt[abs(mesh1->quad_map[z][p])]++;
                  else if (mesh1->quad_map[z][p] == -num_procs && my_rank != 0)
                    sendcnt[0]++;
              break;
            default:
              break;
          }
        }

        //now, send and recv count for all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
 
          MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
 
          MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
 
        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        delete [] srequest;
        delete [] rrequest;
        delete [] statuses;
 
        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];
        
        //allocate buffers to send/recv elem/node indices
        bdim2 = new int[num_procs];
        for (n = 0; n < num_procs; n++)
          bdim2[n] = MAX(sendcnt[n],recvcnt[n])*10*sizeof(int); //indices (2 qualifiers and up to 8 nodes)
        sbuff = new char*[num_procs];
        rbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
        {
          sbuff[n] = new char[bdim2[n]];
          rbuff[n] = new char[bdim2[n]];
        }
#endif
		
        //alloc vars to hold elem type
        int tt = -10;
        int prt = -20;
        int pyt = -30;
        int ht = -40;
        //another benefit here is they can be used on all procs the same
        int *trit = 0;
        int *quadt = 0;
        
        trit = (int*)calloc(mesh1->nb,sizeof(int));
        quadt = (int*)calloc(mesh1->nb,sizeof(int));
        
        q = -50; //start
        //init to higher than four basic
        for (n = 0; n < mesh1->nb; n++)
          {
          trit[n] = q;
          q-=10;
          }
        for (n = 0; n < mesh1->nb; n++)
          {
          quadt[n] = q;
          q-=10;
          }

#ifdef PARALLEL
        //now, pack elem indices as per procs that will be sending them
        for (q = 0; q < num_procs; q++)
        {
        
          if (sendcnt[q] == 0 || q == my_rank)
            continue;

          //set position to 0
          sposition=0;
 
          //pack elem type, index (for return) and nodes
          for (n = 0; n < 4; n++)
          {
          
            //loop over elements
            switch(n)
            {
              case 0:
                for (p = mesh1->ntet; p < ntet; p++)
                {
                  if (mesh1->tet_map[p] < 0 && mesh1->tet_map[p] > -num_procs && abs(mesh1->tet_map[p]) == q)
                  {
                    MPI_Pack(&(tt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 4; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->tet_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->tet_map[p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(tt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 4; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->tet_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                break;
              case 1:
                for (p = mesh1->npri; p < npri; p++)
                {
                  if (mesh1->pri_map[p] < 0 && mesh1->pri_map[p] > -num_procs && abs(mesh1->pri_map[p]) == q)
                  {
                    MPI_Pack(&(prt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 6; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->pri_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->pri_map[p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(prt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 6; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->pri_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                break;
              case 2:
                for (p = mesh1->npyr; p < npyr; p++)
                {
                  if (mesh1->pyr_map[p] < 0 && mesh1->pyr_map[p] > -num_procs && abs(mesh1->pyr_map[p]) == q)
                  {
                    MPI_Pack(&(pyt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 5; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->pyr_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->pyr_map[p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(pyt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 5; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->pyr_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                break;
              case 3:
                for (p = mesh1->nhex; p < nhex; p++)
                {
                  if (mesh1->hex_map[p] < 0 && mesh1->hex_map[p] > -num_procs && abs(mesh1->hex_map[p]) == q)
                  {
                    MPI_Pack(&(ht),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 8; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->hex_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->hex_map[p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(ht),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (z = 0; z < 8; z++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->hex_n[p][z]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                break;
              default:
                break;
            }
          }
          
          //pack elem type, index (for return) and nodes
          for (n = 0; n < 2; n++)
          {
            //loop over elements
            switch(n)
            {
              case 0:
                for (z = 0; z < mesh1->nb; z++)
                {
                for (p = tn[z]; p < mesh1->nt[z]; p++)
                {
                  if (mesh1->tri_map[z][p] < 0 && mesh1->tri_map[z][p] > -num_procs && abs(mesh1->tri_map[z][p]) == q)
                  {
                    MPI_Pack(&(trit[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (r = 0; r < 3; r++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->t_n[z][p][r]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->tri_map[z][p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(trit[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (r = 0; r < 3; r++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->t_n[z][p][r]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                }
                break;
              case 1:
                for (z = 0; z < mesh1->nb; z++)
                {
                for (p = qn[z]; p < mesh1->nq[z]; p++)
                {
                  if (mesh1->quad_map[z][p] < 0 && mesh1->quad_map[z][p] > -num_procs && abs(mesh1->quad_map[z][p]) == q)
                  {
                    MPI_Pack(&(quadt[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (r = 0; r < 4; r++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->q_n[z][p][r]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                  if (mesh1->quad_map[z][p] == -num_procs && q == 0)
                  {
                    MPI_Pack(&(quadt[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    for (r = 0; r < 4; r++)
                    {
                      MPI_Pack(&(mesh1->pmap[mesh1->q_n[z][p][r]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                    }
                  }
                }
                }
                break;
              default:
                break;
            }
          }  
        }
#endif

        free(trit);
        free(quadt);
        
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
   
        //now, send and recv packets from all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (sendcnt[p] > 0 && p != my_rank) 
          {
            MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
 
          if (recvcnt[p] > 0 && p != my_rank) 
          {
            MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }

        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);
 
        delete [] srequest;
        delete [] rrequest;
        delete [] statuses;
 
        srequest = new MPI_Request[num_procs];
        rrequest = new MPI_Request[num_procs];
        statuses = new MPI_Status[num_procs];
 
        //now, delete sbuff and resize
        for (n = 0; n < num_procs; n++)
          if (sbuff[n] != 0) delete [] sbuff[n];
        delete [] sbuff;
 
        //resize
        bdim21 = new int[num_procs];
        for (n = 0; n < num_procs; n++)
          bdim21[n] = MAX(sendcnt[n],recvcnt[n])*3*sizeof(int); //nodes and indices
        sbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
          sbuff[n] = new char[bdim21[n]];
 
        temp_node = temp_node2 = 0;

        List check; //will hold nodes we wish to check

        //init check
        check.Redimension(0);

        int el = 0; //will hold matching element 
 
        //now, unpack type, index, nodes, replace with type, index, global index
        for (p = 0; p < num_procs; p++)
        {
          if (recvcnt[p] == 0 || p == my_rank)
            continue;
            
          sposition = rposition = 0; //reset pos
          for (n = 0; n < recvcnt[p]; n++)
          {
            //unpack type initially
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

            //switch
            switch(temp_node)
            {
              case -10:
                //unpack index to send back
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                for (q = 0; q < 4; q++)
                {
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                  check.Add_To_List(el);
                }
                //check.print(out_f);
                flag = 0;
                for (q = mesh1->ntet; q < ntet && !flag; q++)
                {
                  el = 0; //reset
                  for (z = 0; z < 4; z++)
                  {
                    if (check.Is_In_List(mesh1->pmap[mesh1->tet_n[q][z]][0]))
                    {
                      el++;
                    }
                  }
                  if (el == 4)
                  {
                    el = mesh1->tet_map[q];
                    flag = 1;
                  } else if (el != 4 && q == ntet-1)
                  {
                    fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                    fflush(stderr);
                  }
                }
                if (flag)
                {
                MPI_Pack(&(tt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                }
                check.Redimension(0); //reset list
                break;
              case -20:
                //unpack index to send back
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                for (q = 0; q < 6; q++)
                {
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                  check.Add_To_List(el);
                }
                //check.print(out_f);
                flag = 0;
                for (q = mesh1->npri; q < npri && !flag; q++)
                {
                  el = 0; //reset
                  for (z = 0; z < 6; z++)
                  {
                    if (check.Is_In_List(mesh1->pmap[mesh1->pri_n[q][z]][0]))
                    {
                      el++;
                    }
                  }
                  if (el == 6)
                  {
                    el = mesh1->pri_map[q];
                    flag = 1;
                  } else if (el != 6 && q == npri-1)
                  {
                    fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                    fflush(stderr);
                  }
                }
                if (flag)
                {
                MPI_Pack(&(prt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                }
                check.Redimension(0); //reset list
                break;
              case -30:
                //unpack index to send back
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                for (q = 0; q < 5; q++)
                {
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                  check.Add_To_List(el);
                }
                //check.print(out_f);
                flag = 0;
                for (q = mesh1->npyr; q < npyr && !flag; q++)
                {
                  el = 0; //reset
                  for (z = 0; z < 5; z++)
                  {
                    if (check.Is_In_List(mesh1->pmap[mesh1->pyr_n[q][z]][0]))
                    {
                      el++;
                    }
                  }
                  if (el == 5)
                  {
                    el = mesh1->pyr_map[q];
                    flag = 1;
                  } else if (el != 5 && q == npyr-1)
                  {
                    fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                    fflush(stderr);
                  }
                }
                if (flag)
                {
                MPI_Pack(&(pyt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                }
                check.Redimension(0); //reset list
                break;
              case -40:
                //unpack index to send back
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                for (q = 0; q < 8; q++)
                {
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                  check.Add_To_List(el);
                }
                //check.print(out_f);
                flag = 0;
                for (q = mesh1->nhex; q < nhex && !flag; q++)
                {
                  el = 0; //reset
                  for (z = 0; z < 8; z++)
                  {
                    if (check.Is_In_List(mesh1->pmap[mesh1->hex_n[q][z]][0]))
                    {
                      el++;
                    }
                  }
                  if (el == 8)
                  {
                    el = mesh1->hex_map[q];
                    flag = 1;
                  } else if (el != 8 && q == nhex-1)
                  {
                    fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                    fflush(stderr);
                  }
                }
                if (flag)
                {
                MPI_Pack(&(ht),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                }
                check.Redimension(0); //reset list
                break;
              //will use this for boundary element types
              default:
                if ((temp_node <= -50) && (temp_node > (-50 - 10*mesh1->nb)))
                  {
                  //unpack index to send back
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  for (q = 0; q < 3; q++)
                    {
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                    check.Add_To_List(el);
                    }
                  //check.print(out_f);
                  flag = 0;
                  for (q = tn[(-50 - temp_node)/10]; q < mesh1->nt[(-50 - temp_node)/10] && !flag; q++)
                    {
                    el = 0; //reset
                    for (z = 0; z < 3; z++)
                      {
                      if (check.Is_In_List(mesh1->pmap[mesh1->t_n[(-50 - temp_node)/10][q][z]][0]))
                        {
                        el++;
                        }
                      }
                    if (el == 3)
                      {
                      el = mesh1->tri_map[(-50 - temp_node)/10][q];
                      flag = 1;
                      } 
                    else if (el != 3 && q == mesh1->nt[(-50 - temp_node)/10]-1)
                      {  
                      fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                      fflush(stderr);
                      }
                    }
                  if (flag)
                  {
                  MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  }
                  check.Redimension(0); //reset list
                  }
                if ((temp_node <= (-50 - 10*mesh1->nb)) && (temp_node > (-50 - 20*mesh1->nb)))
                  {
                  //unpack index to send back
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  for (q = 0; q < 4; q++)
                    {
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                    check.Add_To_List(el);
                    }
                  //check.print(out_f);
                  flag = 0;
                  for (q = qn[(-50 - 10*mesh1->nb - temp_node)/10]; q < mesh1->nq[(-50 - 10*mesh1->nb - temp_node)/10] && !flag; q++)
                    {
                    el = 0; //reset
                    for (z = 0; z < 4; z++)
                      {
                      if (check.Is_In_List(mesh1->pmap[mesh1->q_n[(-50 - 10*mesh1->nb - temp_node)/10][q][z]][0]))
                        {
                        el++;
                        }
                      }
                    if (el == 4)
                      {
                      el = mesh1->quad_map[(-50 - 10*mesh1->nb - temp_node)/10][q];
                      flag = 1;
                      } 
                    else if (el != 4 && q == mesh1->nq[(-50 - 10*mesh1->nb - temp_node)/10]-1)
                      {  
                      fprintf(stderr,"\nWARNING: You have a bad element sent that cannot be found!\n");
                      fflush(stderr);
                      }
                    }
                  if (flag)
                  {
                  MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  }
                  check.Redimension(0); //reset list
                  }
                break;
            }
          }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        //now, delete rbuff and bdim2, resize
        delete [] bdim2;
        for (n = 0; n < num_procs; n++)
          if (rbuff[n] != 0) delete [] rbuff[n];
            delete [] rbuff;

        //resize
        rbuff = new char*[num_procs];
        for (n = 0; n < num_procs; n++)
          rbuff[n] = new char[bdim21[n]];
  	
        //now, send and recv packets from all procs
        nreq_s = nreq_r = 0;
        for (p = 0; p < num_procs; p++)
        {
          if (recvcnt[p] > 0 && p != my_rank) 
          {
            MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
  
          if (sendcnt[p] > 0 && p != my_rank) 
          {
            MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
  
        //now, wait for finish before unpacking
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        //finally, unpack nodes in proper position
        for (p = 0; p < num_procs; p++)
        {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;
        
          rposition = 0; //reset pos
          for (n = 0; n < sendcnt[p]; n++)
          {
            //unpack type initially
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

            //switch
            switch(temp_node)
            {
              case -10:
                //unpack index use
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //unpack global elem number
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->tet_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                break;
              case -20:
                //unpack index use
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //unpack global elem number
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->pri_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                break;
              case -30:
                //unpack index use
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //unpack global elem number
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->pyr_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                break;
              case -40:
                //unpack index use
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //unpack global elem number
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->hex_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                break;
              //will use this for boundary element types
              default:
                if ((temp_node <= -50) && (temp_node > (-50 - 10*mesh1->nb)))
                  {
                  //unpack index to use
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack global elem number
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->tri_map[(-50 - temp_node)/10][temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                  }
                if ((temp_node <= (-50 - 10*mesh1->nb)) && (temp_node > (-50 - 20*mesh1->nb)))
                  {
                  //unpack index to use
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack global elem number
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(mesh1->quad_map[(-50 - 10*mesh1->nb - temp_node)/10][temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
                  }
                break;
            }
          }
        }	

        //finally, free MPI mem
        delete[] sendcnt;
        delete[] recvcnt;
        delete[] srequest;
        delete[] rrequest;
        delete[] statuses;
        for (n = 0; n < num_procs; n++)
        {
          if (sbuff[n] != 0) delete [] sbuff[n];
          if (rbuff[n] != 0) delete [] rbuff[n];
        }
        delete[] sbuff;
        delete[] rbuff;
        delete[] bdim21;
#endif
        }

        //free orig qn, tn
        free(qn);
        free(tn);

        //finally, reset element counts
        mesh1->ntet = ntet;
        mesh1->npyr = npyr;
        mesh1->npri = npri;
        mesh1->nhex = nhex;
        
#ifdef PARALLEL
        if (num_procs > 1)
        {
        bool eflag = false;
        //check for negatives and numbers too large
        for (n = 0; n < 4; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (p = 0; p < mesh1->ntet; p++)
              {
                if (mesh1->tet_map[p] < 0 || mesh1->tet_map[p] >= glmax[n])
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has map %d which is out of range.\n",p,mesh1->tet_map[p]);
                  fflush(stderr);
                  eflag = true;
                }
                for (q = 0; q < 4; q++)
                {
                //already checked all nodes, so looking for negative indices or nodes that don't match up
                  if (mesh1->tet_n[p][q] < 0)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has node %d(L = %d) which is -1.\n",p,q,mesh1->tet_n[p][q]);
                    fflush(stderr);
                    eflag = true;
                  }
                  if (mesh1->pmap[mesh1->tet_n[p][q]][0] < 0 || mesh1->pmap[mesh1->tet_n[p][q]][0] >= globalmax)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,mesh1->tet_n[p][q],mesh1->pmap[mesh1->tet_n[p][q]][0]);
                    fflush(stderr);
                    eflag = true;
                  }
                }
              }
              break;
            case 1:
              for (p = 0; p < mesh1->npri; p++)
              {
                if (mesh1->pri_map[p] < 0 || mesh1->pri_map[p] >= glmax[n])
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has map %d which is out of range.\n",p,mesh1->pri_map[p]);
                  fflush(stderr);
                  eflag = true;
                }
                for (q = 0; q < 6; q++)
                {
                  //already checked all nodes, so looking for negative indices or nodes that don't match up
                  if (mesh1->pri_n[p][q] < 0)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has node %d(L = %d) which is -1.\n",p,q,mesh1->pri_n[p][q]);
                    fflush(stderr);
                    eflag = true;
                  }
                  if (mesh1->pmap[mesh1->pri_n[p][q]][0] < 0 || mesh1->pmap[mesh1->pri_n[p][q]][0] >= globalmax)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,mesh1->pri_n[p][q],mesh1->pmap[mesh1->pri_n[p][q]][0]);
                    fflush(stderr);
                    eflag = true;
                  }
                }
              }
              break;
            case 2:
              for (p = 0; p < mesh1->npyr; p++)
              {
                if (mesh1->pyr_map[p] < 0 || mesh1->pyr_map[p] >= glmax[n])
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has map %d which is out of range.\n",p,mesh1->pyr_map[p]);
                  fflush(stderr);
                  eflag = true;
                }
                for (q = 0; q < 5; q++)
                {
                //already checked all nodes, so looking for negative indices or nodes that don't match up
                  if (mesh1->pyr_n[p][q] < 0)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has node %d(L = %d) which is -1.\n",p,q,mesh1->pyr_n[p][q]);
                    fflush(stderr);
                    eflag = true;
                  }
                  if (mesh1->pmap[mesh1->pyr_n[p][q]][0] < 0 || mesh1->pmap[mesh1->pyr_n[p][q]][0] >= globalmax)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,mesh1->pyr_n[p][q],mesh1->pmap[mesh1->pyr_n[p][q]][0]);
                    fflush(stderr);
                    eflag = true;
                  }
                }
              }
              break;
            case 3:
              for (p = 0; p < mesh1->nhex; p++)
              {
                if (mesh1->hex_map[p] < 0 || mesh1->hex_map[p] >= glmax[n])
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has map %d which is out of range.\n",p,mesh1->hex_map[p]);
                  fflush(stderr);
                  eflag = true;
                }
                for (q = 0; q < 8; q++)
                {
                //already checked all nodes, so looking for negative indices or nodes that don't match up
                  if (mesh1->hex_n[p][q] < 0)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has node %d(L = %d) which is -1.\n",p,q,mesh1->hex_n[p][q]);
                    fflush(stderr);
                    eflag = true;
                  }
                  if (mesh1->pmap[mesh1->hex_n[p][q]][0] < 0 || mesh1->pmap[mesh1->hex_n[p][q]][0] >= globalmax)
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,mesh1->hex_n[p][q],mesh1->pmap[mesh1->hex_n[p][q]][0]);
                    fflush(stderr);
                    eflag = true;
                  }
                }
              }
              break;
            default:
              break;
          }
        }
        
        //bd elem
        for (n = 0; n < 2; n++)
        {
          //loop over elements
          switch(n)
          {
            case 0:
              for (z = 0; z < mesh1->nb; z++)
              {
                for (p = 0; p < mesh1->nt[z]; p++)
                {
                  if (mesh1->tri_map[z][p] < 0 || mesh1->tri_map[z][p] >= glmax[4 + z])
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has map %d which is out of range.\n",p,z,mesh1->tri_map[z][p]);
                    fflush(stderr);
                    eflag = true;
                  }
                  for (q = 0; q < 3; q++)
                  {
                  //already checked all nodes, so looking for negative indices or nodes that don't match up
                    if (mesh1->t_n[z][p][q] < 0)
                    {
                      fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has node %d(L = %d) which is -1.\n",p,z,q,mesh1->t_n[z][p][q]);
                      fflush(stderr);
                      eflag = true;
                    }
                    if (mesh1->pmap[mesh1->t_n[z][p][q]][0] < 0 || mesh1->pmap[mesh1->t_n[z][p][q]][0] >= globalmax)
                    {
                      fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has node %d(L = %d) which is out of range globally @ %d.\n",p,z,q,mesh1->t_n[z][p][q],mesh1->pmap[mesh1->t_n[z][p][q]][0]);
                      fflush(stderr);
                      eflag = true;
                    }
                  }
                }
              }
              break;
            case 1:
              for (z = 0; z < mesh1->nb; z++)
              {
                for (p = 0; p < mesh1->nq[z]; p++)
                {
                  if (mesh1->quad_map[z][p] < 0 || mesh1->quad_map[z][p] >= glmax[4 + mesh1->nb + z])
                  {
                    fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has map %d which is out of range.\n",p,z,mesh1->quad_map[z][p]);
                    fflush(stderr);
                    eflag = true;
                  }
                  for (q = 0; q < 4; q++)
                  {
                  //already checked all nodes, so looking for negative indices or nodes that don't match up
                    if (mesh1->q_n[z][p][q] < 0)
                    {
                      fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has node %d(L = %d) which is -1.\n",p,z,q,mesh1->q_n[z][p][q]);
                      fflush(stderr);
                      eflag = true;
                    }
                  if (mesh1->pmap[mesh1->q_n[z][p][q]][0] < 0 || mesh1->pmap[mesh1->q_n[z][p][q]][0] >= globalmax)
                    {
                      fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has node %d(L = %d) which is out of range globally @ %d.\n",p,z,q,mesh1->q_n[z][p][q],mesh1->pmap[mesh1->q_n[z][p][q]][0]);
                      fflush(stderr);
                      eflag = true;
                    }
                  }
                }
              }
              break;
            default:
              break;
          }
        }
        if (eflag)
        {
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        }
        
        }
#endif
        //since used for debug, freed here
        if (glmax > 0)
          free(glmax);
        
        
        //P_VLI:  this can be done after all communication since nothing is added, only moved.
        // force boundary nodes of adjacent boundaries to lie on geometry chain (P_VLI: does same on all procs, owned and no owned)
        for (b=0; b < mesh1->nb-1; b++)
        {
          if (!bft[b])
            continue;

          for (n=0; n < mesh1->nn; n++)
            tag[n] = 0;
          for (i=0; i < mesh1->nt[b]; i++)
          {
            n0 = mesh1->t_n[b][i][0];
            n1 = mesh1->t_n[b][i][1];
            n2 = mesh1->t_n[b][i][2];
            tag[n0] = (b+1);
            tag[n1] = (b+1);
            tag[n2] = (b+1);
          }
          for (i=0; i < mesh1->nq[b]; i++)
          {
            n0 = mesh1->q_n[b][i][0];
            n1 = mesh1->q_n[b][i][1];
            n2 = mesh1->q_n[b][i][2];
            n3 = mesh1->q_n[b][i][3];
            tag[n0] = (b+1);
            tag[n1] = (b+1);
            tag[n2] = (b+1);
            tag[n3] = (b+1);
          }
          for (j=b+1; j < mesh1->nb; j++)
          {
            if (!bft[j])
              continue;
            for (i=0; i < mesh1->nt[j]; i++)
            {
              n0 = mesh1->t_n[j][i][0];
              n1 = mesh1->t_n[j][i][1];
              n2 = mesh1->t_n[j][i][2];
              if (tag[n0] == (b+1)) tag[n0] = (j+1);
              if (tag[n1] == (b+1)) tag[n1] = (j+1);
              if (tag[n2] == (b+1)) tag[n2] = (j+1);
            }
            for (i=0; i < mesh1->nq[j]; i++)
            {
              n0 = mesh1->q_n[j][i][0];
              n1 = mesh1->q_n[j][i][1];
              n2 = mesh1->q_n[j][i][2];
              n3 = mesh1->q_n[j][i][3];
              if (tag[n0] == (b+1)) tag[n0] = (j+1);
              if (tag[n1] == (b+1)) tag[n1] = (j+1);
              if (tag[n2] == (b+1)) tag[n2] = (j+1);
              if (tag[n3] == (b+1)) tag[n3] = (j+1);
            }
            for (n=0; n < mesh1->nn; n++)
            {
              if (tag[n] != (j+1))
                continue;
              i = b;
              m = j;
              mesh1->nodep[n] = geom->closest_chain(mesh1->nodep[n],i,m);
            }
          }
        }

        n=0;
        mesh1->exchange_tag(tag,n);

        // replace nodes
        // since this is inserting under the existing mesh, the computational mesh
        // is no longer valid, so overwrite computational mesh with current physical mesh
        for (n=0; n < mesh1->nn; n++)
        {
          mesh1->nodec[n] = mesh1->nodep[n];
        }

        for (n=0; n < mesh1->nn; n++)
          if (stack[n] != 0)
            delete stack[n];
        delete[] stack;
      } else
      {
        delete[] nrm;
        delete[] desired;
      }

      free(tag);
      
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif

    }  // end of boundary grouping loop
    
    if (restart == 1)
    {
      if (fflag == 1)
      {
        mesh1->smooth_io(1,geom,fname);
        delete mesh1;
        exit(0);
      } else
        mesh1->smooth_io(1,geom,sname);
    }

    if (v_layers > 0)
      if (my_rank == 0) fprintf(out_f,"\nInsertion of layer %d completed.\n",lay);
    else
      if (my_rank == 0) fprintf(out_f,"\nInsertion of layers completed.\n");
    if (my_rank == 0) fflush(out_f);
    
    if (mxnew == 0 && mxlay < lay) // reset loop counter to skip empty layers
    {
      lay = mxlay+1;
      if (my_rank == 0) fprintf(out_f,"\nRESETTING LAYER COUNTER TO %d TO SKIP EMPTY LAYERS!",lay);
      if (my_rank == 0) fflush(out_f);
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

  } // end of global layer loop
  
  if (fflag == 1)
  {
    mesh1->smooth_io(1,geom,fname);
    delete mesh1;
    exit(0);
  } else if (restart == 0)
    mesh1->smooth_io(1,geom,sname);
  
  delete mesh1;

  delete[] bft;

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  return(1);
}
