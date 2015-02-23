#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "geometry.h"
#include "HUGG.h"
#include "Vector.h"
#include "viscous_insert.h"
#include "journal.h"
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

//global variables
int my_rank; /*rank of process*/
int num_procs; /*number of processes*/
FILE *in_f, *jou_f, *out_f; /*input output journal files global*/
int displace = 0;
int E_mode = 1;
double kappa = 0.0;

int main(int argcs, char* pArgs[])
{
  int source; //rank of sender
  int dest; //rank of receiver
  int tag = 0; //tag for messages
  MPI_Status status; //return status for receive
  int i, j, k, nf, lsmoo, osmoo;
  int digits, cnvrg, nsearch, restart;
  const int bdim = 132; //buffer dim
  char extension[bdim]; //file extension to allow padded digits
  char** fnames; //file name storage
  char buff[bdim]; //buffer
  char sname[bdim]; //mesh name storage
  char dname[bdim]; //displacement file name storage
  char param_file[bdim]; // parameter file name
  time_t tm;
  char *t_char;
 
  //Start up MPI
  MPI_Init(&argcs, &pArgs);
   
  //Find out process rank of current instance
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
  //Find out number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  in_f = stdin;
  out_f = stdout;

  if (my_rank == 0)
  {
    //create journal file
    if ((jou_f=fopen("P_Move.jou","w")) == NULL)
    {
      printf("\nCouldn't open file journal file");
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  
    //check for standard input (if not used as a library)
    if (--argcs < 1)
    {
      printf("\nNo input file specified!");
      printf("\nUsing standard input!");
    } else
    {
      if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file <%s>\n",pArgs[argcs]);
        fprintf(stderr,"\n\nUsage: P_Move.XXX [batch_input_file]\n");
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }

  if (my_rank == 0 && in_f != stdin)
  {
    sprintf(buff,"P_Move.out");

    if ((out_f=fopen(buff,"w")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file output file %s",buff);
      if (my_rank == 0)
        fclose(jou_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //perform time check (if not used as a library)
  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0)
  {
    fprintf(out_f,"\nP_Move run started at %s",t_char);

    //print program info 
    fprintf(out_f,"\n======================================================================");
    fprintf(out_f,"\n COPYRIGHT 2003-2012 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA       ");
    fprintf(out_f,"\n                                                                      ");
    fprintf(out_f,"\n                    RIGHTS IN DATA                                    ");
    fprintf(out_f,"\n                                                                      ");
    fprintf(out_f,"\n THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT   ");
    fprintf(out_f,"\n    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO         ");
    fprintf(out_f,"\n         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL        ");
    fprintf(out_f,"\n              RIGHTS IN DATA CONTRACT CLAUSES.                        ");
    fprintf(out_f,"\n       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY     ");
    fprintf(out_f,"\n              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA              ");
    fprintf(out_f,"\n                                                                      ");
    fprintf(out_f,"\n Parallel Mesh Movement (P_Move)                                      ");
    fprintf(out_f,"\n NOTE: This data includes the UT SimCenter at Chattanooga P_OPT       ");
    fprintf(out_f,"\n code, which was developed under private non-government funding.      ");
    fprintf(out_f,"\n This software is submitted with limited rights to use, reproduce,    ");
    fprintf(out_f,"\n and disclose this data for Government Purposes only.                 ");
    fprintf(out_f,"\n Requests for access to the software for non-governmental purposes    ");
    fprintf(out_f,"\n should be referrred to                                               "); 
    fprintf(out_f,"\n                                                                      ");
    fprintf(out_f,"\n    Dr. Steve Karman                                                  "); 
    fprintf(out_f,"\n    Steve-Karman@utc.edu                                              "); 
    fprintf(out_f,"\n    423-425-5492  or  423-425-5470                                    "); 
    fprintf(out_f,"\n                                                                      ");
    fprintf(out_f,"\n    SimCenter: National Center for Computational Engineering          "); 
    fprintf(out_f,"\n    701 East M. L. King Boulevard                                     "); 
    fprintf(out_f,"\n    Chattanooga, TN 37403                                             "); 
    fprintf(out_f,"\n======================================================================\n");
    fflush(out_f);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    journal(in_f, out_f, jou_f, "#Number of geometry files ->",nf);
    //check for adequate number of geom files
    if (nf <= 0)
    {
      fprintf(out_f,"\nNumber of geometry files must be > 0\n");
      fflush(out_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  
    // inform slave processes
    tag = 0;
    for (dest = 1; dest < num_procs; dest++)
      MPI_Send(&nf, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

    fnames = (char**)malloc(nf*sizeof(char*));
  
    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));
      sprintf(buff,"#File %d ->",i+1);
      journal(in_f, out_f, jou_f, buff, fnames[i]);
      // inform slave processes of file names, using strlen for count
      tag = i+1;
      for (dest = 1; dest < num_procs; dest++)
        MPI_Send(fnames[i], strlen(fnames[i])+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
  } else
  {
    // Obtain number of geometry files from master
    source = 0; //reset to master
    tag = 0;
    MPI_Recv(&nf, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    //read geometry (if not used as a library...normally already done)
    fnames = (char**)malloc(nf*sizeof(char*));
  
    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));

      tag = i+1; //based on master
      MPI_Recv(fnames[i], bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
    }
  }

  int *facet = 0;
  int *mat = 0;
  double *gvert = 0;
  int ngn = 0;
  int nfacets = 0;
  read_geometry(nf, fnames, ngn, nfacets, gvert, facet, mat, out_f);
  
  //create geometry
  sprintf(param_file,"%s","geometry.params");
  geometry* geom = new geometry();
  geom->process_geom(ngn, nfacets, gvert, facet, mat, 0, 0, 0, param_file, out_f);
  
  //clean up fnames
  for (i=0; i < nf; i++)
    free(fnames[i]);
  free(fnames);
  
  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    //THIS IS A GENERIC FILE NAME...proc nuumbers in actual input added later!
    //read in grid file name
    journal(in_f, out_f, jou_f, "#Enter physical grid file name, serial version ->",sname);
  
    // inform slave processes of file name, using strlen for count
    tag = 0;
    for (dest = 1; dest < num_procs; dest++)
      MPI_Send(sname, strlen(sname)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);	
    
    //read in displacement file name
    dname[0] = '\0';
    journal(in_f, out_f, jou_f, "#Enter displacement file name ->",dname);
  
    // inform slave processes of file name, using strlen for count
    tag = 1;
    for (dest = 1; dest < num_procs; dest++)
      MPI_Send(dname, strlen(dname)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);	
    
    journal(in_f, out_f, jou_f, "#Enter # of outer smoothing sweeps > ",osmoo);
  
    journal(in_f, out_f, jou_f, "#Enter # of Linear-Elastic smoothing sweeps > ",lsmoo);
  
    journal(in_f, out_f, jou_f, "#Enter stiffness model (Bitwise flags bit 3-1/ds, bit 2-1/vol, bit 1-AR ) > [0 - 7] >",E_mode);

    journal(in_f, out_f, jou_f, "#Enter kappa [ >= 0.0] > ",kappa);

    journal(in_f, out_f, jou_f, "#Enter convergence OoM [ > 0] > ",cnvrg);
  
    journal(in_f, out_f, jou_f, "#Enter # of GMRES search directions [ > 0] > ",nsearch);
  
    fclose(jou_f);
    if (in_f != stdin) fclose(in_f);
  } else
  {
    source = 0; //reset to master
    tag = 0;
    MPI_Recv(sname, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  
    tag = 1;
    MPI_Recv(dname, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  }

  E_mode = MAX(0,MIN(7,E_mode));
  MPI_Bcast(&E_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&osmoo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lsmoo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cnvrg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nsearch, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kappa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if (my_rank == 0)
  {
    fprintf(out_f,"\n\nUser supplied input parameters:");
    fprintf(out_f,"\nNumber of outer smoothing sweeps = %d",osmoo);
    fprintf(out_f,"\nNumber of Linear-Elastic smoothing sweeps = %d",lsmoo);
    fprintf(out_f,"\nStiffness model (Bitwise flags bit 3-1/ds, bit 2-1/vol, bit 1-AR ) = %d",E_mode);
    fprintf(out_f,"\nKappa [ >= 0.0] = %g",kappa);
    fprintf(out_f,"\nConvergence OoM = %d",cnvrg);
    fprintf(out_f,"\nNumber of GMRES search directions = %d",nsearch);
    fflush(out_f);
  }

  mesh_deform(geom,sname,dname,osmoo,lsmoo,nsearch,cnvrg);
 
  //perform end time check (if not used as a library)
  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0) fprintf(out_f,"\nP_Move run completed at %s",t_char);
  
  MPI_Barrier(MPI_COMM_WORLD);

  //delete geom
  delete geom;

  if (my_rank == 0 && out_f != stdout)
    fclose(out_f);

  MPI_Finalize(); 
  
  return(0);
}
