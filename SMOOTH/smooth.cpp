#ifdef PARALLEL
#include "mpi.h"
#endif
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

void Smesh_obj::exchange_tag(int *tag, int min_flag)
{
#ifdef PARALLEL
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  int temp_node, temp_node2, temp_tag, n, p;
  
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
  for (n = 0; n < nn; n++)
    if (pmap[n][1] != my_rank)
      sendcnt[pmap[n][1]]++; 

  MPI_Barrier(MPI_COMM_WORLD);

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

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*(2+min_flag)*sizeof(int); //indices
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (pmap[n][1] == p)
      {
        MPI_Pack(&(pmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        if (min_flag) MPI_Pack(&(tag[n]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  

  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      if (min_flag)
      {
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_tag,1,MPI_INT,MPI_COMM_WORLD);
        tag[temp_node] = MIN(tag[temp_node],temp_tag);
      }
      //pack node 
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(tag[temp_node]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(tag[temp_node]),1,MPI_INT,MPI_COMM_WORLD);
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
  delete[] bdim;
#endif

  return;
}

void Smesh_obj::exchange_double(double *fun)
{
#ifdef PARALLEL
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  int temp_node, temp_node2, n, p;
  
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
  for (n = 0; n < nn; n++)
    if (pmap[n][1] != my_rank)
      sendcnt[pmap[n][1]]++; 

  MPI_Barrier(MPI_COMM_WORLD);

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

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*(sizeof(int)+sizeof(double));
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (pmap[n][1] == p)
      {
        MPI_Pack(&(pmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
 
  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  
  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      //pack node and nrm
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(fun[temp_node]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(fun[temp_node]),1,MPI_DOUBLE,MPI_COMM_WORLD);
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
  delete[] bdim;
#endif

  return;
}

void Smesh_obj::exchange_Vector(Vector *vec)
{
#ifdef PARALLEL
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  int temp_node, temp_node2, n, p;
  
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
  for (n = 0; n < nn; n++)
    if (pmap[n][1] != my_rank)
      sendcnt[pmap[n][1]]++; 

  MPI_Barrier(MPI_COMM_WORLD);

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

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*(sizeof(int)+3*sizeof(double));
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (pmap[n][1] == p)
      {
        MPI_Pack(&(pmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
 
  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  
  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      //pack node and nrm
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(vec[temp_node][0]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(vec[temp_node][1]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(vec[temp_node][2]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(vec[temp_node][0]),1,MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(vec[temp_node][1]),1,MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(vec[temp_node][2]),1,MPI_DOUBLE,MPI_COMM_WORLD);
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
  delete[] bdim;
#endif

  return;
}

int Smesh_obj::smooth_io(int mode, geometry *geom, char sname[])
{
  int n, flag;
  
  flag = 0;

  switch (mode)
  {
    case -1:
      //read in comp and physical meshes
      if (read_comp_physical(sname) != 0)
      {
        fprintf(stderr,"\nCouldn't open file %s",sname);
        fflush(stderr);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }

      if (nb != geom->ngb)
      {
        fprintf(stderr,"\nNumber of boundaries in mesh does not match geometry file!");
        fflush(stderr);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
      break;

    case 1:

      // write new physical mesh file

      mesh_stats();
  
      if (write_mesh(nodep,sname) != 0)
      {
        fprintf(stderr,"\nCouldn't open file %s",sname);
        fflush(stderr);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
      break;
  }

  return flag;
}

int Smesh_obj::read_comp_physical(char sname[])
{
  int b, n, i;
  //open and read computational mesh CGNS file
  int flag = 0;
  int parallel = 1;
  //flag not needed anymore, but these are reset to nonsensical values to assure read went well
  nn = nb = ntet = npyr = npri = nhex = -1;
  //actually reads in files and allocates memory for the object
  if (num_procs > 1)
  {
    flag = P_CGNS_read(parallel,sname,nn,&nodep,nb,&bname,&nt,&t_n,&nq,&q_n,ntet,&tet_n,
                       npyr,&pyr_n,npri,&pri_n,nhex,&hex_n,&pmap,&tri_map,&quad_map,
                       &tet_map,&pyr_map,&pri_map,&hex_map);
  } else
  {
    //no need for maps since no need to reassemble
    flag = CGNS_read(sname,nn,&nodep,nb,&bname,&nt,&t_n,&nq,&q_n,ntet,&tet_n,npyr,&pyr_n,
                     npri,&pri_n,nhex,&hex_n);
    //allocate for pmap
    pmap = (int**)calloc(nn,sizeof(int*));
    for (i = 0; i < nn; i++)
      pmap[i] = (int*)calloc(3,sizeof(int));
    //fill out serially
    for (i = 0; i < nn; i++)
    {
      pmap[i][0] = i;
      pmap[i][1] = 0;
      pmap[i][2] = i;
    }
  }

  int gnn, gntet, gnpyr, gnpri, gnhex, gnt, gnq;

  gnn=nn;
  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;

  int local;

  if (num_procs > 1)
  {
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    gnn=local=0;
    for (n=0; n < nn; n++)
      local = MAX(local,pmap[n][0]+1);
    gnn=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnn,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
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
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of nodes      = %d",gnn);
    fprintf(out_f,"\nNumber of tetrahedra = %d",gntet);
    fprintf(out_f,"\nNumber of pyramid    = %d",gnpyr);
    fprintf(out_f,"\nNumber of prism      = %d",gnpri);
    fprintf(out_f,"\nNumber of hexahedra  = %d",gnhex);
    fprintf(out_f,"\nNumber of boundaries = %d",nb);
    fflush(out_f);
  }

  for (b=0; b < nb; b++)
  {
    gnt=nt[b];
    gnq=nq[b];
#ifdef PARALLEL
    if (num_procs > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      gnt=local=0;
      for (n=0; n < nt[b]; n++)
        local = MAX(local,tri_map[b][n]+1);
      MPI_Allreduce(&local,&gnt,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      gnq=local=0;
      for (n=0; n < nq[b]; n++)
        local = MAX(local,quad_map[b][n]+1);
      MPI_Allreduce(&local,&gnq,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    }
#endif
    if (my_rank == 0)
    {
      fprintf(out_f,"\nBoundary %d, name = %s",b+1,bname[b]);
      fprintf(out_f,"\nTotal number of triangles      = %d",gnt);
      fprintf(out_f,"\nTotal number of quadrilaterals = %d",gnq);
      fflush(out_f);
    }
  }
  
  //for opt-based smooth, physical and computational match, so we just duplicate nodes
  nodec = (Point*)malloc(nn*sizeof(Point));
  
  for (i = 0; i < nn; i++)
    nodec[i] = nodep[i];

  //return flag for success
  return (flag);
}


int Smesh_obj::write_mesh(Point *node, char sname[])
{
  int flag = 0; 
  int parallel = 1;
  
  if (num_procs > 1)
  {
    flag=P_CGNS_write(parallel,sname,nn,node,nb,bname,nt,t_n,nq,q_n,ntet,tet_n,npyr,pyr_n,
                      npri,pri_n,nhex,hex_n,pmap,tri_map,quad_map,tet_map,pyr_map,pri_map,hex_map);
  } else
  {
    //no need for maps since no need to reassemble
    flag=CGNS_write(sname,nn,node,nb,bname,nt,t_n,nq,q_n,ntet,tet_n,npyr,pyr_n,
                    npri,pri_n,nhex,hex_n); 
  }

  return (flag);
}

void Smesh_obj::mesh_stats()
{
  int b, c, i, j, s, n;
  int *tmp;
  int gnn, gntet, gnpyr, gnpri, gnhex, gnt, gnq;

#ifdef PARALLEL
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  gnn=nn;
  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;

  int local;

  if (num_procs > 1)
  {
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    gnn=local=0;
    for (n=0; n < nn; n++)
      local = MAX(local,pmap[n][0]+1);
    gnn=local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnn,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
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
  if (my_rank == 0)
  {
    fprintf(out_f,"\nTotal number of nodes = %d",gnn);
    fprintf(out_f,"\nTotal number of tetrahedra = %d",gntet);
    fprintf(out_f,"\nTotal number of pyramids = %d",gnpyr);
    fprintf(out_f,"\nTotal number of prisms = %d",gnpri);
    fprintf(out_f,"\nTotal number of hexahedra = %d",gnhex);
    fprintf(out_f,"\nTotal number of boundaries = %d",nb);
    fflush(out_f);
  }
  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;
#ifdef PARALLEL
  // reset global numbers for averaging operations below
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntet,&gntet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npyr,&gnpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npri,&gnpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nhex,&gnhex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
#endif

  int n0, n1, n2, n3, n4, n5, n6, n7, gint, dest, tag, source;
  double vmn, vmx, avg, vol, jmn, jmx, jvg, gdouble;
  Point p0, p1, p2, p3, p4, p5, p6, p7;
  Point pmn, pmx, qmn, qmx;
  double amn, amx, area, ds;
  Vector v1, v2, v3, v4, nrm;
  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 

  tmp = new int[nn];

  for (b=0; b < nb; b++)
  {
    if (my_rank == 0) fprintf(out_f,"\n\nStatistics for boundary %d : %s",b+1,bname[b]);

    for (n=0; n < nn; n++)
      tmp[n] = 0;

    for (c=0; c < nt[b]; c++)
      for (i=0; i < 3; i++)
        tmp[t_n[b][c][i]] = 1;
    for (c=0; c < nq[b]; c++)
      for (i=0; i < 4; i++)
        tmp[q_n[b][c][i]] = 1;

    amn = 1.0e20;
    amx = -1.0e20;
    avg = 0.0;
    jmn = 1.0e20;
    jmx = -1.0e20;
    jvg = 0.0;
    i=j=0;

    for (c=0; c < ntet; c++)
    {
      for (s=0; s < 4; s++)
      {
        switch (s)
        {
          case 0:
            n0 = tet_n[c][0];
            n1 = tet_n[c][1];
            n2 = tet_n[c][2];
            n3 = tet_n[c][3];
            break;
          case 1:
            n0 = tet_n[c][0];
            n1 = tet_n[c][3];
            n2 = tet_n[c][1];
            n3 = tet_n[c][2];
            break;
          case 2:
            n0 = tet_n[c][1];
            n1 = tet_n[c][3];
            n2 = tet_n[c][2];
            n3 = tet_n[c][0];
            break;
          case 3:
            n0 = tet_n[c][0];
            n1 = tet_n[c][2];
            n2 = tet_n[c][3];
            n3 = tet_n[c][1];
            break;
          default:
            break;
        }
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        
        if (tmp[n0]+tmp[n1]+tmp[n2] == 3)
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p0,p2);
          v3 = Vector(p1,p2);
          ds = v1.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p1)/2.0;
          }
          jvg += ds;
          j++;
          ds = v2.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p2)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p2)/2.0;
          }
          jvg += ds;
          j++;
          ds = v3.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p1+p2)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p1+p2)/2.0;
          }
          jvg += ds;
          j++;
          nrm = v1 % v2;
          nrm.normalize();
          v3 = Vector(p0,p3);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2)/3.0;
          }
          avg += ds;
          i++;
        }
      }
    }

    for (c=0; c < npyr; c++)
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
      if (tmp[n0]+tmp[n1]+tmp[n2]+tmp[n3] == 4)
      {
        v1 = Vector(p0,p1);
        v2 = Vector(p1,p2);
        v3 = Vector(p2,p3);
        v4 = Vector(p3,p0);
        ds = v1.magnitude();
        if (ds < jmn)
        {
          jmn = ds;
          qmn = (p0+p1)/2.0;
        }
        if (ds > jmx)
        {
          jmx = ds;
          qmx = (p0+p1)/2.0;
        }
        jvg += ds;
        j++;
        ds = v2.magnitude();
        if (ds < jmn)
        {
          jmn = ds;
          qmn = (p1+p2)/2.0;
        }
        if (ds > jmx)
        {
          jmx = ds;
          qmx = (p1+p2)/2.0;
        }
        jvg += ds;
        j++;
        ds = v3.magnitude();
        if (ds < jmn)
        {
          jmn = ds;
          qmn = (p2+p3)/2.0;
        }
        if (ds > jmx)
        {
          jmx = ds;
          qmx = (p2+p3)/2.0;
        }
        jvg += ds;
        j++;
        ds = v4.magnitude();
        if (ds < jmn)
        {
          jmn = ds;
          qmn = (p0+p3)/2.0;
        }
        if (ds > jmx)
        {
          jmx = ds;
          qmx = (p0+p3)/2.0;
        }
        jvg += ds;
        j++;
        v1 = Vector(p0,p2);
        v2 = Vector(p1,p3);
        nrm = v1 % v2;
        nrm.normalize();
        v3 = Vector((p0+p1+p2+p3)*0.25,p4);
        ds = fabs(v3 * nrm);
        if (ds < amn)
        {
          amn = ds;
          pmn = (p0+p1+p2+p3)/4.0;
        }
        if (ds > amx)
        {
          amx = ds;
          pmx = (p0+p1+p2+p3)/4.0;
        }
        avg += ds;
        i++;
      }
      for (s=1; s < 5; s++)
      {
        switch (s)
        {
          case 1:
            n0 = pyr_n[c][0];
            n1 = pyr_n[c][1];
            n2 = pyr_n[c][2];
            n3 = pyr_n[c][3];
            n4 = pyr_n[c][4];
            break;
          case 2:
            n0 = pyr_n[c][1];
            n1 = pyr_n[c][2];
            n2 = pyr_n[c][3];
            n3 = pyr_n[c][0];
            n4 = pyr_n[c][4];
            break;
          case 3:
            n0 = pyr_n[c][2];
            n1 = pyr_n[c][3];
            n2 = pyr_n[c][0];
            n3 = pyr_n[c][1];
            n4 = pyr_n[c][4];
            break;
          case 4:
            n0 = pyr_n[c][3];
            n1 = pyr_n[c][0];
            n2 = pyr_n[c][1];
            n3 = pyr_n[c][2];
            n4 = pyr_n[c][4];
            break;
          default:
            break;
        }
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        if (tmp[n0]+tmp[n1]+tmp[n4] == 3)
        {
          v1 = Vector(p4,p1);
          v2 = Vector(p4,p0);
          v3 = Vector(p0,p1);
          ds = v1.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p4+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p4+p1)/2.0;
          }
          jvg += ds;
          j++;
          ds = v2.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p4+p0)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p4+p0)/2.0;
          }
          jvg += ds;
          j++;
          ds = v3.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p1)/2.0;
          }
          jvg += ds;
          j++;

          nrm = v1 % v2;
          nrm.normalize();
          v3 = Vector((p0+p1+p4)/3.0,p2);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p4)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p4)/3.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p4)/3.0,p3);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p4)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p4)/3.0;
          }
          avg += ds;
          i++;
        }
      }
    }

    for (c=0; c < npri; c++)
    {
      for (s=0; s < 3; s++)
      {
        switch (s)
        {
          case 0:
            n0 = pri_n[c][0];
            n1 = pri_n[c][1];
            n2 = pri_n[c][2];
            n3 = pri_n[c][3];
            n4 = pri_n[c][4];
            n5 = pri_n[c][5];
            break;
          case 1:
            n0 = pri_n[c][1];
            n1 = pri_n[c][2];
            n2 = pri_n[c][0];
            n3 = pri_n[c][4];
            n4 = pri_n[c][5];
            n5 = pri_n[c][3];
            break;
          case 2:
            n0 = pri_n[c][2];
            n1 = pri_n[c][0];
            n2 = pri_n[c][1];
            n3 = pri_n[c][5];
            n4 = pri_n[c][3];
            n5 = pri_n[c][4];
            break;
          default:
            break;
        }
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        p5 = nodep[n5];
        if (tmp[n0]+tmp[n1]+tmp[n3]+tmp[4] == 4)
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p1,p4);
          v3 = Vector(p4,p3);
          v4 = Vector(p3,p0);
          ds = v1.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p1)/2.0;
          }
          jvg += ds;
          j++;
          ds = v2.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p1+p4)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p1+p4)/2.0;
          }
          jvg += ds;
          j++;
          ds = v3.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p3+p4)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p3+p4)/2.0;
          }
          jvg += ds;
          j++;
          ds = v4.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p3)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p3)/2.0;
          }
          jvg += ds;
          j++;
          v1 = Vector(p0,p4);
          v2 = Vector(p3,p1);
          nrm = v1 % v2;
          nrm.normalize();
          v3 = Vector((p0+p1+p3+p4)/4.0,p2);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p3+p4)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p3+p4)/4.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p3+p4)/4.0,p5);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p3+p4)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p3+p4)/4.0;
          }
          avg += ds;
          i++;
        }
      }
      for (s=3; s < 5; s++)
      {
        switch (s)
        {
          case 3:
            n0 = pri_n[c][0];
            n1 = pri_n[c][1];
            n2 = pri_n[c][2];
            n3 = pri_n[c][3];
            n4 = pri_n[c][4];
            n5 = pri_n[c][5];
            break;
          case 4:
            n0 = pri_n[c][5];
            n1 = pri_n[c][4];
            n2 = pri_n[c][3];
            n3 = pri_n[c][2];
            n4 = pri_n[c][1];
            n5 = pri_n[c][0];
            break;
          default:
            break;
        }
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        p5 = nodep[n5];
        if (tmp[n0]+tmp[n1]+tmp[n2] == 3)
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p0,p2);
          v3 = Vector(p1,p2);
          ds = v1.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p1)/2.0;
          }
          jvg += ds;
          j++;
          ds = v2.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p2)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p2)/2.0;
          }
          jvg += ds;
          j++;
          ds = v3.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p1+p2)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p1+p2)/2.0;
          }
          jvg += ds;
          j++;
          nrm = v1 % v2;
          nrm.normalize();
          v3 = Vector((p0+p1+p2)/3.0,p3);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2)/3.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p2)/3.0,p4);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2)/3.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p2)/3.0,p5);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2)/3.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2)/3.0;
          }
          avg += ds;
          i++;
        }
      }
    }
    for (c=0; c < nhex; c++)
    {
      for (s=0; s < 6; s++)
      {
        switch (s)
        {
          case 0:
            n0 = hex_n[c][0];
            n1 = hex_n[c][1];
            n2 = hex_n[c][2];
            n3 = hex_n[c][3];
            n4 = hex_n[c][4];
            n5 = hex_n[c][5];
            n6 = hex_n[c][6];
            n7 = hex_n[c][7];
            break;
          case 1:
            n0 = hex_n[c][0];
            n1 = hex_n[c][4];
            n2 = hex_n[c][5];
            n3 = hex_n[c][1];
            n4 = hex_n[c][3];
            n5 = hex_n[c][7];
            n6 = hex_n[c][6];
            n7 = hex_n[c][2];
            break;
          case 2:
            n0 = hex_n[c][1];
            n1 = hex_n[c][5];
            n2 = hex_n[c][6];
            n3 = hex_n[c][2];
            n4 = hex_n[c][0];
            n5 = hex_n[c][4];
            n6 = hex_n[c][7];
            n7 = hex_n[c][3];
            break;
          case 3:
            n0 = hex_n[c][2];
            n1 = hex_n[c][6];
            n2 = hex_n[c][7];
            n3 = hex_n[c][3];
            n4 = hex_n[c][1];
            n5 = hex_n[c][5];
            n6 = hex_n[c][4];
            n7 = hex_n[c][0];
            break;
          case 4:
            n0 = hex_n[c][3];
            n1 = hex_n[c][7];
            n2 = hex_n[c][4];
            n3 = hex_n[c][0];
            n4 = hex_n[c][2];
            n5 = hex_n[c][6];
            n6 = hex_n[c][5];
            n7 = hex_n[c][1];
            break;
          case 5:
            n0 = hex_n[c][7];
            n1 = hex_n[c][6];
            n2 = hex_n[c][5];
            n3 = hex_n[c][4];
            n4 = hex_n[c][3];
            n5 = hex_n[c][2];
            n6 = hex_n[c][1];
            n7 = hex_n[c][0];
            break;
          default:
            break;
        }
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        p5 = nodep[n5];
        p6 = nodep[n6];
        p7 = nodep[n7];
        if (tmp[n0]+tmp[n1]+tmp[n2]+tmp[n3] == 4)
        {
          v1 = Vector(p0,p1);
          v2 = Vector(p1,p2);
          v3 = Vector(p2,p3);
          v4 = Vector(p3,p0);
          ds = v1.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p0+p1)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p0+p1)/2.0;
          }
          jvg += ds;
          j++;
          ds = v2.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p1+p2)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p1+p2)/2.0;
          }
          jvg += ds;
          j++;
          ds = v3.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p2+p3)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p2+p3)/2.0;
          }
          jvg += ds;
          j++;
          ds = v4.magnitude();
          if (ds < jmn)
          {
            jmn = ds;
            qmn = (p3+p0)/2.0;
          }
          if (ds > jmx)
          {
            jmx = ds;
            qmx = (p3+p0)/2.0;
          }
          jvg += ds;
          j++;
          v1 = Vector(p0,p2);
          v2 = Vector(p1,p3);
          nrm = v1 % v2;
          nrm.normalize();
          v3 = Vector((p0+p1+p2+p3)/4.0,p4);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2+p3)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2+p3)/4.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p2+p3)/4.0,p5);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2+p3)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2+p3)/4.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p2+p3)/4.0,p6);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2+p3)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2+p3)/4.0;
          }
          avg += ds;
          i++;
          v3 = Vector((p0+p1+p2+p3)/4.0,p7);
          ds = fabs(v3 * nrm);
          if (ds < amn)
          {
            amn = ds;
            pmn = (p0+p1+p2+p3)/4.0;
          }
          if (ds > amx)
          {
            amx = ds;
            pmx = (p0+p1+p2+p3)/4.0;
          }
          avg += ds;
          i++;
        }
      }
    }

    gdouble=avg;
    gint=i;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&i,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/MAX(1,gint);
    gdouble=jvg;
    gint=j;
#ifdef PARALLEL
    MPI_Allreduce(&jvg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&j,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    jvg = gdouble/MAX(1,gint);

#ifdef PARALLEL
    MPI_Allreduce(&jmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = jmn;
    cmn_g_in.proc = my_rank;
    jmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&qmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&qmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
    MPI_Allreduce(&jmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = jmx;
    cmx_g_in.proc = my_rank;
    jmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&qmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&qmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&qmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&qmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tangential spacing = %g @ (%g, %g, %g)",jmn,qmn[0],qmn[1],qmn[2]);
      fprintf(out_f,"\nAverage tangential spacing = %g",jvg);
      fprintf(out_f,"\nMaximum tangential spacing = %g @ (%g, %g, %g)",jmx,qmx[0],qmx[1],qmx[2]);
    }

#ifdef PARALLEL
    MPI_Allreduce(&amn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = amn;
    cmn_g_in.proc = my_rank;
    amn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
    MPI_Allreduce(&amx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = amx;
    cmx_g_in.proc = my_rank;
    amx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum normal spacing = %g @ (%g, %g, %g)",amn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage normal spacing = %g",avg);
      fprintf(out_f,"\nMaximum normal spacing = %g @ (%g, %g, %g)",amx,pmx[0],pmx[1],pmx[2]);
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    amn = 1.0e20;
    amx = -1.0e20;
    avg = 0.0;
    if (nt[b] > 0)
    {
      for (c=0; c < nt[b]; c++)
      {
        n0 = t_n[b][c][0];
        n1 = t_n[b][c][1];
        n2 = t_n[b][c][2];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        area = triangle_area(p0,p1,p2);
        if (area < amn)
        {
          amn = area;
          pmn = (p0+p1+p2)/3.0;
        }
        if (area > amx)
        {
          amx = area;
          pmx = (p0+p1+p2)/3.0;
        }
        avg += area;
      }
    }
    gdouble=avg;
    gint=nt[b];
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nt[b],&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (gint > 0)
    {
      avg = gdouble/MAX(1,gint);

#ifdef PARALLEL
      MPI_Allreduce(&amn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
      cmn_g_in.val = amn;
      cmn_g_in.proc = my_rank;
      amn = gdouble;
      MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (cmn_g_out.proc != 0)
      {
        if (my_rank == cmn_g_out.proc)
        {
          tag = my_rank;
          dest = 0;
          MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        }
        if (my_rank == 0)
        {
          tag = source = cmn_g_out.proc;
          MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        }
      }
      MPI_Allreduce(&amx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      cmx_g_in.val = amx;
      cmx_g_in.proc = my_rank;
      amx = gdouble;
      MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (cmx_g_out.proc != 0)
      {
        if (my_rank == cmx_g_out.proc)
        {
          tag = my_rank;
          dest = 0;
          MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        }
        if (my_rank == 0)
        {
          tag = source = cmx_g_out.proc;
          MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        }
      }
#endif

      if (my_rank == 0)
      {
        fprintf(out_f,"\n");
        fprintf(out_f,"\nMinimum triangle area = %g @ (%g, %g, %g)",amn,pmn[0],pmn[1],pmn[2]);
        fprintf(out_f,"\nAverage triangle area = %g",avg);
        fprintf(out_f,"\nMaximum triangle area = %g @ (%g, %g, %g)",amx,pmx[0],pmx[1],pmx[2]);
      }
    }

    amn = 1.0e20;
    amx = -1.0e20;
    avg = 0.0;
    if (nq[b] > 0)
    {
      for (c=0; c < nq[b]; c++)
      {
        n0 = q_n[b][c][0];
        n1 = q_n[b][c][1];
        n2 = q_n[b][c][2];
        n3 = q_n[b][c][3];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        area = quadrilateral_area(p0,p1,p2,p3);
        if (area < amn)
        {
          amn = area;
          pmn = (p0+p1+p2+p3)/4.0;
        }
        if (area > amx)
        {
          amx = area;
          pmx = (p0+p1+p2+p3)/4.0;
        }
        avg += area;
      }
    }
    gdouble=avg;
    gint=nq[b];
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nq[b],&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (gint > 0)
    {
      avg = gdouble/MAX(1,gint);

#ifdef PARALLEL
      MPI_Allreduce(&amn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
      cmn_g_in.val = amn;
      cmn_g_in.proc = my_rank;
      amn = gdouble;
      MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (cmn_g_out.proc != 0)
      {
        if (my_rank == cmn_g_out.proc)
        {
          tag = my_rank;
          dest = 0;
          MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        }
        if (my_rank == 0)
        {
          tag = source = cmn_g_out.proc;
          MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        }
      }
      MPI_Allreduce(&amx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      cmx_g_in.val = amx;
      cmx_g_in.proc = my_rank;
      amx = gdouble;
      MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (cmx_g_out.proc != 0)
      {
        if (my_rank == cmx_g_out.proc)
        {
          tag = my_rank;
          dest = 0;
          MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
          MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        }
        if (my_rank == 0)
        {
          tag = source = cmx_g_out.proc;
          MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        }
      }
#endif

      if (my_rank == 0)
      {
        fprintf(out_f,"\n");
        fprintf(out_f,"\nMinimum quadrilateral area = %g @ (%g, %g, %g)",amn,pmn[0],pmn[1],pmn[2]);
        fprintf(out_f,"\nAverage quadrilateral area = %g",avg);
        fprintf(out_f,"\nMaximum quadrilateral area = %g @ (%g, %g, %g)",amx,pmx[0],pmx[1],pmx[2]);
      }
    }
  }
  delete[] tmp;

  double *hist;

  c=MAX(MAX(MAX(ntet*4,npyr*5),npri*6),nhex*8);
  c=MAX(c,1);
  hist = new double[c];

  if (gntet > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < ntet; c++)
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
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gntet;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tetrahedral volume = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage tetrahedral volume = %g",avg);
      fprintf(out_f,"\nMaximum tetrahedral volume = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(ntet,hist);
  }
    
  if (gnpyr > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < npyr; c++)
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
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4)/5.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4)/5.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnpyr;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum pyramid volume = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage pyramid volume = %g",avg);
      fprintf(out_f,"\nMaximum pyramid volume = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(npyr,hist);
  }

  if (gnpri > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < npri; c++)
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
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4+p5)/6.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4+p5)/6.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnpri;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum prism volume = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage prism volume = %g",avg);
      fprintf(out_f,"\nMaximum prism volume = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(npri,hist);
  }

  if (gnhex > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
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
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      vol = hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnhex;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum hexahderal volume = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage hexahderal volume = %g",avg);
      fprintf(out_f,"\nMaximum hexahderal volume = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(nhex,hist);
  }

  if (gntet > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < ntet; c++)
    {
      n0 = tet_n[c][0];
      n1 = tet_n[c][1];
      n2 = tet_n[c][2];
      n3 = tet_n[c][3];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      vol = tetrahedral_aspect_ratio(p0,p1,p2,p3);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gntet;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tetrahedral aspect ratio = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage tetrahedral aspect ratio = %g",avg);
      fprintf(out_f,"\nMaximum tetrahedral aspect ratio = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(ntet,hist);
  }
    
  if (gnpyr > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < npyr; c++)
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
      vol = pyramid_aspect_ratio(p0,p1,p2,p3,p4);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4)/5.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4)/5.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnpyr;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum pyramid aspect ratio = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage pyramid aspect ratio = %g",avg);
      fprintf(out_f,"\nMaximum pyramid aspect ratio = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(npyr,hist);
  }

  if (gnpri > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < npri; c++)
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
      vol = prism_aspect_ratio(p0,p1,p2,p3,p4,p5);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4+p5)/6.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4+p5)/6.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnpri;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum prism aspect ratio = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage prism aspect ratio = %g",avg);
      fprintf(out_f,"\nMaximum prism aspect ratio = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(npri,hist);
  }

  if (gnhex > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
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
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      vol = hexahedral_aspect_ratio(p0,p1,p2,p3,p4,p5,p6,p7);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gnhex;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum hexahderal aspect ratio = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage hexahderal aspect ratio = %g",avg);
      fprintf(out_f,"\nMaximum hexahderal aspect ratio = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(nhex,hist);
  }

  if (gntet > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
    for (c=0; c < ntet; c++)
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
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p0;
      }
      avg += vol;
      v1 = Vector(p1,p0);
      v2 = Vector(p1,p2);
      v3 = Vector(p1,p3);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p1;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p1;
      }
      avg += vol;
      v1 = Vector(p2,p0);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p3);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p2;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p2;
      }
      avg += vol;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p1);
      v3 = Vector(p3,p2);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p3;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p3;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(4*gntet);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tetrahedral solid angle = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage tetrahedral solid angle = %g",avg);
      fprintf(out_f,"\nMaximum tetrahedral solid angle = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gnpyr > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
    for (c=0; c < npyr; c++)
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
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p0;
      }
      avg += vol;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p1;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p1;
      }
      avg += vol;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p4);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p2;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p2;
      }
      avg += vol;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p4);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p3;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p3;
      }
      avg += vol;
      v1 = Vector(p4,p0);
      v2 = Vector(p4,p1);
      v3 = Vector(p4,p2);
      vol = solid_angle(v1,v2,v3);
      v1 = Vector(p4,p0);
      v2 = Vector(p4,p2);
      v3 = Vector(p4,p3);
      vol += solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p4;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p4;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(5*gnpyr);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum pyramid solid angle = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage pyramid solid angle = %g",avg);
      fprintf(out_f,"\nMaximum pyramid solid angle = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gnpri > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
    for (c=0; c < npri; c++)
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
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p0;
      }
      avg += vol;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p1;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p1;
      }
      avg += vol;
      v1 = Vector(p2,p0);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p5);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p2;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p2;
      }
      avg += vol;
      v1 = Vector(p3,p5);
      v2 = Vector(p3,p4);
      v3 = Vector(p3,p0);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p3;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p3;
      }
      avg += vol;
      v1 = Vector(p4,p3);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p1);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p4;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p4;
      }
      avg += vol;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p3);
      v3 = Vector(p5,p2);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p5;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p5;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(6*gnpri);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum prism solid angle = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage prism solid angle = %g",avg);
      fprintf(out_f,"\nMaximum prism solid angle = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gnhex > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
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
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p0;
      }
      avg += vol;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p5);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p1;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p1;
      }
      avg += vol;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p6);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p2;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p2;
      }
      avg += vol;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p7);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p3;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p3;
      }
      avg += vol;
      v1 = Vector(p4,p7);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p0);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p4;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p4;
      }
      avg += vol;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p6);
      v3 = Vector(p5,p1);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p5;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p5;
      }
      avg += vol;
      v1 = Vector(p6,p5);
      v2 = Vector(p6,p7);
      v3 = Vector(p6,p2);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p6;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p6;
      }
      avg += vol;
      v1 = Vector(p7,p6);
      v2 = Vector(p7,p4);
      v3 = Vector(p7,p3);
      vol = solid_angle(v1,v2,v3);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = p7;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = p7;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(8*gnhex);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum hexahderal solid angle = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage hexahderal solid angle = %g",avg);
      fprintf(out_f,"\nMaximum hexahderal solid angle = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gntet > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    for (c=0; c < ntet; c++)
    {
      n0 = tet_n[c][0];
      n1 = tet_n[c][1];
      n2 = tet_n[c][2];
      n3 = tet_n[c][3];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      vol = tetrahedral_condition_number(p0,p1,p2,p3,1);
      hist[c] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/gntet;

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tetrahedral condition number = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage tetrahedral condition number = %g",avg);
      fprintf(out_f,"\nMaximum tetrahedral condition number = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(ntet,hist);
  }

  if (gnpyr > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
    for (c=0; c < npyr; c++)
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
      vol = tetrahedral_condition_number(p0,p1,p3,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p3+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p3+p4)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p1,p2,p0,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p4)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p2,p3,p1,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p1+p2+p3+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p1+p2+p3+p4)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p3,p0,p2,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p2+p3+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p2+p3+p4)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(4*gnpyr);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum pyramid corner condition number = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage pyramid corner condition number = %g",avg);
      fprintf(out_f,"\nMaximum pyramid corner condition number = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gnpri > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
    for (c=0; c < npri; c++)
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
      vol = tetrahedral_condition_number(p0,p1,p2,p3,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p3)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p3)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p1,p2,p0,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p4)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p2,p0,p1,p5,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p5)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p5)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p3,p5,p4,p0,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p3+p4+p5)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p3+p4+p5)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p4,p3,p5,p1,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p1+p3+p4+p5)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p1+p3+p4+p5)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p5,p4,p3,p2,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p2+p3+p4+p5)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p2+p3+p4+p5)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(6*gnpri);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum prism corner condition number = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage prism corner condition number = %g",avg);
      fprintf(out_f,"\nMaximum prism corner condition number = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  if (gnhex > 0)
  {
    vmn = 1.0e20;
    vmx = -1.0e20;
    avg = 0.0;
    i=0;
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
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      vol = tetrahedral_condition_number(p0,p1,p3,p4,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p3+p4)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p3+p4)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p1,p2,p0,p5,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p1+p2+p5)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p1+p2+p5)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p2,p3,p1,p6,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p1+p2+p3+p6)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p1+p2+p3+p6)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p3,p0,p2,p7,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p2+p3+p7)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p2+p3+p7)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p4,p7,p5,p0,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p0+p4+p5+p7)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p0+p4+p5+p7)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p5,p4,p6,p1,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p1+p4+p5+p6)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p1+p4+p5+p6)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p6,p5,p7,p2,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p2+p5+p6+p7)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p2+p5+p6+p7)/4.0;
      }
      avg += vol;
      vol = tetrahedral_condition_number(p7,p6,p4,p3,0);
      hist[i++] = vol;
      if (vol < vmn)
      {
        vmn = vol;
        pmn = (p3+p4+p6+p7)/4.0;
      }
      if (vol > vmx)
      {
        vmx = vol;
        pmx = (p3+p4+p6+p7)/4.0;
      }
      avg += vol;
    }
    gdouble=avg;
#ifdef PARALLEL
    MPI_Allreduce(&avg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    avg = gdouble/(8*gnhex);

#ifdef PARALLEL
    MPI_Allreduce(&vmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmn;
    cmn_g_in.proc = my_rank;
    vmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&vmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = vmx;
    cmx_g_in.proc = my_rank;
    vmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum hexahedral corner condition number = %g @ (%g, %g, %g)",
             vmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage hexahedral corner condition number = %g",avg);
      fprintf(out_f,"\nMaximum hexahedral corner condition number = %g @ (%g, %g, %g)",
             vmx,pmx[0],pmx[1],pmx[2]);
    }

    if (num_procs == 1) histogram(i,hist);
  }

  delete[] hist;

  int neg, neg_skw, pos_skw, pos;
  Point cg;
  double jac;

  if (gntet > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;

    jmn = 1.0e20;
    jmx = -1.0e20;
    jvg = 0.0;
    for (c=0; c < ntet; c++)
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
      if (jac < 0.0)
        neg++;
      else
        pos++;
      if (jac < jmn)
      {
        jmn = jac;
        pmn = (p0+p1+p2+p3)/4.0;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = (p0+p1+p2+p3)/4.0;
      }
      jvg += jac;
    }
#ifdef PARALLEL
    MPI_Allreduce(&neg,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg = gint;
    MPI_Allreduce(&neg_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg_skw = gint;
    MPI_Allreduce(&pos_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos_skw = gint;
    MPI_Allreduce(&pos,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos = gint;
#endif
    
    gdouble=jvg;
#ifdef PARALLEL
    MPI_Allreduce(&jvg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    jvg = gdouble/gntet;

#ifdef PARALLEL
    MPI_Allreduce(&jmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = jmn;
    cmn_g_in.proc = my_rank;
    jmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&jmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = jmx;
    cmx_g_in.proc = my_rank;
    jmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum tetrahedral Jacobian = %g @ (%g, %g, %g)",jmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage tetrahedral Jacobian = %g",jvg);
      fprintf(out_f,"\nMaximum tetrahedral Jacobian = %g @ (%g, %g, %g)",jmx,pmx[0],pmx[1],pmx[2]);
      fprintf(out_f,"\n");
      fprintf(out_f,"\nTetrahedral Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
      fprintf(out_f,"\n");
      fflush(out_f);
    }
  }

  if (gnpyr > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;

    jmn = 1.0e20;
    jmx = -1.0e20;
    jvg = 0.0;
    for (c=0; c < npyr; c++)
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
      avg = vmn = vmx = (jac = scalar_triple_product(v1,v2,v3));
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p0;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p0;
      }
      jvg += jac;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p1;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p1;
      }
      jvg += jac;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p2;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p2;
      }
      jvg += jac;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p3;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p3;
      }
      jvg += jac;
      avg *= 0.25;

      if (vmx < 0.0)
        neg++;
      else if (avg < 0.0 && vmx > 0.0)
        neg_skw++;
      else if (avg > 0.0 && vmn < 0.0)
        pos_skw++;
      else
        pos++;
    }
#ifdef PARALLEL
    MPI_Allreduce(&neg,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg = gint;
    MPI_Allreduce(&neg_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg_skw = gint;
    MPI_Allreduce(&pos_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos_skw = gint;
    MPI_Allreduce(&pos,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos = gint;
#endif

    gdouble=jvg;
#ifdef PARALLEL
    MPI_Allreduce(&jvg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    jvg = gdouble/(4*gnpyr);

#ifdef PARALLEL
    MPI_Allreduce(&jmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = jmn;
    cmn_g_in.proc = my_rank;
    jmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&jmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = jmx;
    cmx_g_in.proc = my_rank;
    jmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum pyramid Jacobian = %g @ (%g, %g, %g)",jmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage pyramid Jacobian = %g",jvg);
      fprintf(out_f,"\nMaximum pyramid Jacobian = %g @ (%g, %g, %g)",jmx,pmx[0],pmx[1],pmx[2]);
      fprintf(out_f,"\n");
      fprintf(out_f,"\nPyramid Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
      fprintf(out_f,"\n");
      fflush(out_f);
    }
  }

  if (gnpri > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;

    jmn = 1.0e20;
    jmx = -1.0e20;
    jvg = 0.0;
    for (c=0; c < npri; c++)
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
      avg = vmn = vmx = (jac = scalar_triple_product(v1,v2,v3));
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p0;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p0;
      }
      jvg += jac;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p1;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p1;
      }
      jvg += jac;
      v1 = Vector(p2,p0);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p2;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p2;
      }
      jvg += jac;
      v1 = Vector(p3,p5);
      v2 = Vector(p3,p4);
      v3 = Vector(p3,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p3;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p3;
      }
      jvg += jac;
      v1 = Vector(p4,p3);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p4;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p4;
      }
      jvg += jac;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p3);
      v3 = Vector(p5,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p5;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p5;
      }
      jvg += jac;

      avg /= 6.0;

      if (vmx < 0.0)
        neg++;
      else if (avg < 0.0 && vmx > 0.0)
        neg_skw++;
      else if (avg > 0.0 && vmn < 0.0)
        pos_skw++;
      else
        pos++;
    }
#ifdef PARALLEL
    MPI_Allreduce(&neg,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg = gint;
    MPI_Allreduce(&neg_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg_skw = gint;
    MPI_Allreduce(&pos_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos_skw = gint;
    MPI_Allreduce(&pos,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos = gint;
#endif

    gdouble=jvg;
#ifdef PARALLEL
    MPI_Allreduce(&jvg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    jvg = gdouble/(6*gnpri);

#ifdef PARALLEL
    MPI_Allreduce(&jmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = jmn;
    cmn_g_in.proc = my_rank;
    jmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&jmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = jmx;
    cmx_g_in.proc = my_rank;
    jmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum prism Jacobian = %g @ (%g, %g, %g)",jmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage prism Jacobian = %g",jvg);
      fprintf(out_f,"\nMaximum prism Jacobian = %g @ (%g, %g, %g)",jmx,pmx[0],pmx[1],pmx[2]);
      fprintf(out_f,"\n");
      fprintf(out_f,"\nPrism Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
      fprintf(out_f,"\n");
      fflush(out_f);
    }
  }

  if (gnhex > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;

    jmn = 1.0e20;
    jmx = -1.0e20;
    jvg = 0.0;
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
      avg = vmn = vmx = (jac = scalar_triple_product(v1,v2,v3));
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p0;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p0;
      }
      jvg += jac;
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p1;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p1;
      }
      jvg += jac;
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p6);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p2;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p2;
      }
      jvg += jac;
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p7);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p3;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p3;
      }
      jvg += jac;
      v1 = Vector(p4,p7);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p4;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p4;
      }
      jvg += jac;
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p6);
      v3 = Vector(p5,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p5;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p5;
      }
      jvg += jac;
      v1 = Vector(p6,p5);
      v2 = Vector(p6,p7);
      v3 = Vector(p6,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p6;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p6;
      }
      jvg += jac;
      v1 = Vector(p7,p6);
      v2 = Vector(p7,p4);
      v3 = Vector(p7,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      avg += (jac = scalar_triple_product(v1,v2,v3));
      vmn = MIN(vmn,jac);
      vmx = MAX(vmx,jac);
      if (jac < jmn)
      {
        jmn = jac;
        pmn = p7;
      }
      if (jac > jmx)
      {
        jmx = jac;
        pmx = p7;
      }
      jvg += jac;

      avg *= 0.125;

      if (vmx < 0.0)
        neg++;
      else if (avg < 0.0 && vmx > 0.0)
        neg_skw++;
      else if (avg > 0.0 && vmn < 0.0)
        pos_skw++;
      else
        pos++;
    }
#ifdef PARALLEL
    MPI_Allreduce(&neg,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg = gint;
    MPI_Allreduce(&neg_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    neg_skw = gint;
    MPI_Allreduce(&pos_skw,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos_skw = gint;
    MPI_Allreduce(&pos,&gint,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    pos = gint;
#endif

    gdouble=jvg;
#ifdef PARALLEL
    MPI_Allreduce(&jvg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    jvg = gdouble/(8*gnhex);

#ifdef PARALLEL
    MPI_Allreduce(&jmn,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = jmn;
    cmn_g_in.proc = my_rank;
    jmn = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmn[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmn[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&pmn[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmn[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
     
    MPI_Allreduce(&jmx,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    cmx_g_in.val = jmx;
    cmx_g_in.proc = my_rank;
    jmx = gdouble;
    MPI_Allreduce(&cmx_g_in,&cmx_g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (cmx_g_out.proc != 0)
    {
      if (my_rank == cmx_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&pmx[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&pmx[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmx_g_out.proc;
        MPI_Recv(&pmx[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&pmx[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif

    if (my_rank == 0)
    {
      fprintf(out_f,"\n");
      fprintf(out_f,"\nMinimum hexahedral Jacobian = %g @ (%g, %g, %g)",jmn,pmn[0],pmn[1],pmn[2]);
      fprintf(out_f,"\nAverage hexahedral Jacobian = %g",jvg);
      fprintf(out_f,"\nMaximum hexahedral Jacobian = %g @ (%g, %g, %g)",jmx,pmx[0],pmx[1],pmx[2]);
      fprintf(out_f,"\n");
      fprintf(out_f,"\nHexahedral Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
      fprintf(out_f,"\n");
      fflush(out_f);
    }
  }

}

