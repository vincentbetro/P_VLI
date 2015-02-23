#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "geometry.h"
#include "Point.h"
#include "Vector.h"
#include "Util.h"
#include "smooth.h"
#include "HUGG.h"
#include "CGNS.h"
#include "Linked_List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/
extern int displace;

int Smesh_obj::displaced(geometry *geom, char dname[])
{
  FILE* fp;
  int m, n, nr, flag, offset, count;
  int b, i, j, n0, n1, n2, n3;
  const int bdim = 200;
  char buff[bdim];
  double x, y, z, dx, dy, dz, s, tol;

  flag = 0;

  // store current location in computational coordinate
  for (n=0; n < nn; n++)
    nodec[n] = nodep[n];

  if ((fp = fopen(dname,"r")) != 0)
  {
    tol = 1.0e+20;
    for (b=0; b < nb; b++)
    {
      for (j=0; j < nt[b]; j++)
      {
        n0 = t_n[b][j][0];
        n1 = t_n[b][j][1];
        n2 = t_n[b][j][2];
        tol = MIN(tol,distance(nodep[n0],nodep[n1]));
        tol = MIN(tol,distance(nodep[n0],nodep[n2]));
        tol = MIN(tol,distance(nodep[n1],nodep[n2]));
      }
      for (j=0; j < nq[b]; j++)
      {
        n0 = q_n[b][j][0];
        n1 = q_n[b][j][1];
        n2 = q_n[b][j][2];
        n3 = q_n[b][j][3];
        tol = MIN(tol,distance(nodep[n0],nodep[n1]));
        tol = MIN(tol,distance(nodep[n1],nodep[n2]));
        tol = MIN(tol,distance(nodep[n2],nodep[n3]));
        tol = MIN(tol,distance(nodep[n3],nodep[n0]));
      }
    }
    tol *= 0.00001;
    if (my_rank == 0) fprintf(out_f,"\nComputed tolerance for node comparison = %lg",tol);
    if (my_rank == 0) fflush(out_f);

    dx = dy = dz = 0.0;
    if (my_rank == 0) fprintf(out_f,"\nReading displacements from file <%s>.",dname);
    fgets(buff,bdim,fp);
    sscanf(buff,"%d",&offset);
    //sscanf(buff,"%d %lf",&offset,&s);
    if (my_rank == 0) fprintf(out_f,"\nNode offset in file = %d",offset);
    //if (my_rank == 0) fprintf(out_f,"\nScale factor in file = %lf",s);
    if (my_rank == 0) fflush(out_f);
    i=0;
    while (fgets(buff,bdim,fp) != NULL)
    {
      // for m >= 0 an individual node displacement is defined.
      // for m < 0 and dx,dy,dz specified the displacement of body -m is specified
      // for m = -1 and one double the rotation angle about X is specified, additional lines indicate bodies
      // for m = -2 and one double the rotation angle about Y is specified, additional lines indicate bodies
      // for m = -3 and one double the rotation angle about Z is specified, additional lines indicate bodies

      nr=sscanf(buff,"%d %lf %lf %lf",&m,&dx,&dy,&dz);
      if (nr == 4 && m >= 0)
      {
        if (sqrt(dx*dx+dy*dy+dz*dz) < tol)
          continue;

        m -= offset;
        if (num_procs == 1)
          n = m;
        else
        {
          // determine local node number
          n = -1;
          for (j=0; j < nn && n < 0; j++)
            if (pmap[j][0] == m)
              n = j;
        }

        if (n < 0) continue;

        nodec[n] = nodec[n] + Point(dx,dy,dz);

        i++;
      } else if (nr == 4 && m < 0)
      {
        b = (-m)-1;
        if (my_rank == 0) fprintf(out_f,"\nTranslating body %d by (%lf, %lf, %lf)",b+1,dx,dy,dz);
        if (my_rank == 0) fflush(out_f);
        int *tag = new int[nn];
        for (n=0; n < nn; n++)
          tag[n] = 0;
        for (j=0; j < nt[b]; j++)
        {
          n0 = t_n[b][j][0];
          n1 = t_n[b][j][1];
          n2 = t_n[b][j][2];
          tag[n0] = 1;
          tag[n1] = 1;
          tag[n2] = 1;
        }
        for (j=0; j < nq[b]; j++)
        {
          n0 = q_n[b][j][0];
          n1 = q_n[b][j][1];
          n2 = q_n[b][j][2];
          n3 = q_n[b][j][3];
          tag[n0] = 1;
          tag[n1] = 1;
          tag[n2] = 1;
          tag[n3] = 1;
        }
        for (n=0; n < nn; n++)
        {
          if (tag[n])
          {
            nodec[n] = nodec[n] + Point(dx,dy,dz);
            i++;
          }
        }
        delete[] tag;
      } else if (nr == 2 && m < 0)
      {
        Point cg;
        m = abs(m);
        double dang = dx;
        if (my_rank == 0 && m == 1) fprintf(out_f,"\nRotating %lg degrees about X axis",dang);
        if (my_rank == 0 && m == 2) fprintf(out_f,"\nRotating %lg degrees about Y axis",dang);
        if (my_rank == 0 && m == 3) fprintf(out_f,"\nRotating %lg degrees about Z axis",dang);
        double pi = 4.0*atan(1.0);
        dang *= pi/180.0;
        fgets(buff,bdim,fp);
        nr=sscanf(buff,"%lf %lf %lf",&dx,&dy,&dz);
        cg = Point(dx,dy,dz);
        if (my_rank == 0) fprintf(out_f,"\nRotation point = (%lf, %lf, %lf)",cg[0],cg[1],cg[2]);
        if (my_rank == 0) fflush(out_f);
        int *tag = new int[nn];
        for (n=0; n < nn; n++)
          tag[n] = 0;
        
        while (fgets(buff,bdim,fp) != NULL)
        {
          sscanf(buff,"%d",&b);
          if (my_rank == 0 && b > 0) fprintf(out_f,"\nRotating body %d",b);
          if (b > 0)
          {
            b--;
            for (j=0; j < nt[b]; j++)
            {
              n0 = t_n[b][j][0];
              n1 = t_n[b][j][1];
              n2 = t_n[b][j][2];
              tag[n0] = 1;
              tag[n1] = 1;
              tag[n2] = 1;
            }
            for (j=0; j < nq[b]; j++)
            {
              n0 = q_n[b][j][0];
              n1 = q_n[b][j][1];
              n2 = q_n[b][j][2];
              n3 = q_n[b][j][3];
              tag[n0] = 1;
              tag[n1] = 1;
              tag[n2] = 1;
              tag[n3] = 1;
            }
          } else
            break;
        }
        switch (m)
        {
          case 1: // Rotate about X axis
            for (n=0; n < nn; n++)
            {
              if (tag[n] == 0) continue;
              Vector dr = Vector(cg,nodec[n]);
              double ds = sqrt(dr[1]*dr[1]+dr[2]*dr[2]);
              double theta = atan2(dr[2],dr[1]);
              theta += dang;
              nodec[n] = Point(nodec[n][0],cg[1]+ds*cos(theta),cg[2]+ds*sin(theta));
              i++;
            }
            break;
          case 2: // Rotate about Y axis
            for (n=0; n < nn; n++)
            {
              if (tag[n] == 0) continue;
              Vector dr = Vector(cg,nodec[n]);
              double ds = sqrt(dr[0]*dr[0]+dr[2]*dr[2]);
              double theta = atan2(dr[0],dr[2]);
              theta += dang;
              nodec[n] = Point(cg[0]+ds*sin(theta),nodec[n][1],cg[2]+ds*cos(theta));
              i++;
            }
            break;
          case 3: // Rotate about Z axis
            for (n=0; n < nn; n++)
            {
              if (tag[n] == 0) continue;
              Vector dr = Vector(cg,nodec[n]);
              double ds = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
              double theta = atan2(dr[1],dr[0]);
              theta += dang;
              nodec[n] = Point(cg[0]+ds*cos(theta),cg[1]+ds*sin(theta),nodec[n][2]);
              i++;
            }
            break;
          default:
            break;
        }

        delete[] tag;
      }
    }
    MPI_Allreduce(&i,&n,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (my_rank == 0) fprintf(out_f,"\nTotal number of nodes displaced = %d",n);
    if (my_rank == 0) fflush(out_f);
    fclose(fp);

    flag = 1;
  } else
  {
    Point p0, p1, p2, q0, q1;
    Vector vec1, vec2, vec3;
    double ds, ang;
    double *mag = new double[nn];
    double *wgt = new double[nn];
    Vector *nrm = new Vector[nn];
    for (n=0; n < nn; n++)
    {
      nrm[n] = Vector(0.0,0.0,0.0);
      mag[n] = 0.0;
      wgt[n] = 0.0;
    }

    tol = 1.0e+20;
    for (b=0; b < nb; b++)
    {
      for (j=0; j < nt[b]; j++)
      {
        n0 = t_n[b][j][0];
        n1 = t_n[b][j][1];
        n2 = t_n[b][j][2];
        ds = distance(nodep[n0],nodep[n1]);
        tol = MIN(tol,ds);
        mag[n0] += ds;
        wgt[n0] += 1.0;
        mag[n1] += ds;
        wgt[n1] += 1.0;
        ds = distance(nodep[n1],nodep[n2]);
        tol = MIN(tol,ds);
        mag[n1] += ds;
        wgt[n1] += 1.0;
        mag[n2] += ds;
        wgt[n2] += 1.0;
        ds = distance(nodep[n2],nodep[n0]);
        tol = MIN(tol,ds);
        mag[n2] += ds;
        wgt[n2] += 1.0;
        mag[n0] += ds;
        wgt[n0] += 1.0;
        vec1 = Vector(nodep[n0],nodep[n1]);
        vec2 = Vector(nodep[n0],nodep[n2]);
        vec3 = (vec1 % vec2)*0.5;
        vec3.normalize();
        vec1.normalize();
        vec2.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n0] += vec3*ang;
        vec1 = Vector(nodep[n1],nodep[n2]);
        vec2 = Vector(nodep[n1],nodep[n0]);
        vec1.normalize();
        vec2.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n1] += vec3*ang;
        vec1 = Vector(nodep[n2],nodep[n0]);
        vec2 = Vector(nodep[n2],nodep[n1]);
        vec1.normalize();
        vec2.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n2] += vec3*ang;
      }
      for (j=0; j < nq[b]; j++)
      {
        n0 = q_n[b][j][0];
        n1 = q_n[b][j][1];
        n2 = q_n[b][j][2];
        n3 = q_n[b][j][3];
        ds = distance(nodep[n0],nodep[n1]);
        tol = MIN(tol,ds);
        mag[n0] += ds;
        wgt[n0] += 1.0;
        mag[n1] += ds;
        wgt[n1] += 1.0;
        ds = distance(nodep[n1],nodep[n2]);
        tol = MIN(tol,ds);
        mag[n1] += ds;
        wgt[n1] += 1.0;
        mag[n2] += ds;
        wgt[n2] += 1.0;
        ds = distance(nodep[n2],nodep[n3]);
        tol = MIN(tol,ds);
        mag[n2] += ds;
        wgt[n2] += 1.0;
        mag[n3] += ds;
        wgt[n3] += 1.0;
        ds = distance(nodep[n3],nodep[n0]);
        tol = MIN(tol,ds);
        mag[n3] += ds;
        wgt[n3] += 1.0;
        mag[n0] += ds;
        wgt[n0] += 1.0;
        vec1 = Vector(nodep[n0],nodep[n1]);
        vec2 = Vector(nodep[n0],nodep[n3]);
        vec3 = vec1 % vec2;
        vec1.normalize();
        vec2.normalize();
        vec3.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n0] += vec3*ang;
        vec1 = Vector(nodep[n1],nodep[n2]);
        vec2 = Vector(nodep[n1],nodep[n0]);
        vec3 = vec1 % vec2;
        vec1.normalize();
        vec2.normalize();
        vec3.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n1] += vec3*ang;
        vec1 = Vector(nodep[n2],nodep[n3]);
        vec2 = Vector(nodep[n2],nodep[n1]);
        vec3 = vec1 % vec2;
        vec1.normalize();
        vec2.normalize();
        vec3.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n2] += vec3*ang;
        vec1 = Vector(nodep[n3],nodep[n0]);
        vec2 = Vector(nodep[n3],nodep[n2]);
        vec3 = vec1 % vec2;
        vec1.normalize();
        vec2.normalize();
        vec3.normalize();
        ang = acos(MIN(1.0,MAX(0.0,vec1*vec2)));
        nrm[n3] += vec3*ang;
      }
    }
    for (n=0; n < nn; n++)
      if (wgt[n] > 1.0e-10)
      {
        nrm[n].normalize();
        mag[n] = mag[n]/wgt[n]*3.0;
      }

    double ltol = 0.001*tol;
    MPI_Allreduce(&ltol,&tol,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    if (my_rank == 0) fprintf(out_f,"\nComputed tolerance for node comparison = %lg",tol);
    if (my_rank == 0) fflush(out_f);

    for (n=0; n < nn; n++)
      nodec[n] = nodep[n];

    for (b=0; b < nb; b++)
    {
      if (my_rank == 0) fprintf(out_f,"\nProjecting boundary %d nodes to surface.",b+1);
      if (my_rank == 0) fflush(out_f);
      for (i=0; i < nt[b]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[b][i][j];
          m = b;
          p0 = nodep[n] - Point(nrm[n][0],nrm[n][1],nrm[n][2])*tol;
          p1 = nodep[n] + Point(nrm[n][0],nrm[n][1],nrm[n][2])*mag[n];
          q0 = nodep[n] + Point(nrm[n][0],nrm[n][1],nrm[n][2])*tol;
          q1 = nodep[n] - Point(nrm[n][0],nrm[n][1],nrm[n][2])*mag[n];
          if (geom->Intersect(p0,p1,p2,m) == 1)
            nodec[n] = p2;
          else if (geom->Intersect(q0,q1,p2,m) == 1)
            nodec[n] = p2;
          else
          {
            vec1 = Vector(0.0,0.0,0.0);
            nodec[n] = geom->closest(nodep[n],m,vec1);
          }
        }
      }
      for (i=0; i < nq[b]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[b][i][j];
          m = b;
          p0 = nodep[n] - Point(nrm[n][0],nrm[n][1],nrm[n][2])*tol;
          p1 = nodep[n] + Point(nrm[n][0],nrm[n][1],nrm[n][2])*mag[n];
          q0 = nodep[n] + Point(nrm[n][0],nrm[n][1],nrm[n][2])*tol;
          q1 = nodep[n] - Point(nrm[n][0],nrm[n][1],nrm[n][2])*mag[n];
          if (geom->Intersect(p0,p1,p2,m) == 1)
            nodec[n] = p2;
          else if (geom->Intersect(q0,q1,p2,m) == 1)
            nodec[n] = p2;
          else
          {
            vec1 = Vector(0.0,0.0,0.0);
            nodec[n] = geom->closest(nodep[n],m,vec1);
          }
        }
      }
    }
    // force boundary nodes of adjacent boundaries to lie on geometry chain
    for (int g=0; g < nb; g++)
    {
      for (n=0; n < nn; n++)
        wgt[n] = 0.0;
      for (b=0; b < nb; b++)
      {
        if (b==g)
          continue;
        for (i=0; i < nt[b]; i++)
        {
          n0 = t_n[b][i][0];
          n1 = t_n[b][i][1];
          n2 = t_n[b][i][2];
          wgt[n0] = 1.0;
          wgt[n1] = 1.0;
          wgt[n2] = 1.0;
        }
        for (i=0; i < nq[b]; i++)
        {
          n0 = q_n[b][i][0];
          n1 = q_n[b][i][1];
          n2 = q_n[b][i][2];
          n3 = q_n[b][i][3];
          wgt[n0] = 1.0;
          wgt[n1] = 1.0;
          wgt[n2] = 1.0;
          wgt[n3] = 1.0;
        }
      }
      for (i=0; i < nt[g]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[g][i][j];
          if (fabs(wgt[n]-1.0) < 1.e-10) wgt[n] = 2.0;
        }
      }
      for (i=0; i < nq[g]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[g][i][j];
          if (fabs(wgt[n]-1.0) < 1.e-10) wgt[n] = 2.0;
        }
      }
      for (n=0; n < nn; n++)
      {
        if (fabs(wgt[n]-2.0) > 1.e-10)
          continue;
        i = g;
        j = b;
        p1 = nodep[n];
        p2 = geom->closest_chain(p1,i,j);
        nodec[n] = p2;
      }
    }
    
    flag=0;
    for (n=0; n < nn; n++)
    {
      dx = nodec[n][0]-nodep[n][0];
      dy = nodec[n][1]-nodep[n][1];
      dz = nodec[n][2]-nodep[n][2];
      if (fabs(dx) > tol || fabs(dy) > tol || fabs(dz) > tol)
        flag++;
    }

    int global;
    MPI_Allreduce(&flag,&global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (my_rank == 0) fprintf(out_f,"\nTotal number of nodes displaced = %d",global);
    if (my_rank == 0) fflush(out_f);

    flag = MIN(1,global);
    delete[] mag;
    delete[] wgt;
    delete[] nrm;
  }

  return (flag);
}

