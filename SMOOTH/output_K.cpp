#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"
#include "smooth.h"

#include <hdf5.h>

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

void Smesh_obj::output_K(int mdim, double K[][3][3], int ia[], int iau[], int ja[])
{
  char objname[32], filename[32], buff[32];
  int *itmp, i, j, n;
  double *dtmp;
  hid_t file_id, Obj_group_id;
  hid_t dspace_id, dset_id;
  hsize_t dim;
  herr_t status;

  int isize = sizeof(int);
  int dsize = sizeof(double);
  itmp=(int*)malloc(isize);
  dtmp=(double*)malloc(dsize);

  // Initialize HDF
  status = H5open();

  // create object name entry for later use
  sprintf(objname,"HUGG_K");
  
  // create file name
  filename[0]='\0';
  strcat(filename,objname);
  strcat(filename,".hdf");

  // open hdf5 file, overwriting if necessary.
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0)
  {
    fprintf(out_f,"\nError: unable to create file <%s>\n",filename);
    exit(0);
  }

  Obj_group_id = H5Gcreate(file_id, objname, 0);

  dim = 2;
  itmp=(int*)realloc((void*)itmp,(int)dim*isize);
  itmp[0] = nn;
  itmp[1] = mdim;
  fprintf(out_f,"\nWriting stiffness matrix:");
  fprintf(out_f,"\nNumber of nodes  = %d",itmp[0]);
  fprintf(out_f,"\n1st dimension of K = %d",itmp[1]);
  dspace_id = H5Screate_simple(1,&dim,NULL);
  dset_id = H5Dcreate(Obj_group_id,"Parameters",H5T_NATIVE_INT,dspace_id,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,itmp);
  status=H5Sclose(dspace_id);
  status=H5Dclose(dset_id);
  free(itmp);

  dim = nn+1;
  dspace_id=H5Screate_simple(1,&dim,NULL);
  dset_id = H5Dcreate(Obj_group_id,"ia",H5T_NATIVE_INT,dspace_id,H5P_DEFAULT);
  if ((status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ia)) < 0)
  {
    fprintf(out_f,"\nError writing ia array.");
    exit(0);
  }
  status=H5Dclose(dset_id);

  dim = nn;
  dspace_id=H5Screate_simple(1,&dim,NULL);
  dset_id = H5Dcreate(Obj_group_id,"iau",H5T_NATIVE_INT,dspace_id,H5P_DEFAULT);
  if ((status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,iau)) < 0)
  {
    fprintf(out_f,"\nError writing iau array.");
    exit(0);
  }
  status=H5Dclose(dset_id);

  dim = mdim;
  dspace_id=H5Screate_simple(1,&dim,NULL);
  dset_id = H5Dcreate(Obj_group_id,"ja",H5T_NATIVE_INT,dspace_id,H5P_DEFAULT);
  if ((status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ja)) < 0)
  {
    fprintf(out_f,"\nError writing ja array.");
    exit(0);
  }
  status=H5Dclose(dset_id);

  dtmp=(double*)realloc((void*)dtmp,(int)mdim*dsize);
  for (j=0; j < 3; j++)
  {
    for (i=0; i < 3; i++)
    {
      sprintf(buff,"Entry%d_%d",i,j);
      dset_id = H5Dcreate(Obj_group_id,buff,H5T_NATIVE_DOUBLE,dspace_id,H5P_DEFAULT);
      for (n=0; n < mdim; n++)
        dtmp[n]=K[n][i][j];
      if ((status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,dtmp))<0)
      {
        fprintf(out_f,"\nError writing matrix entry %d %d.", i, j);
        exit(0);
      }
      status=H5Dclose(dset_id);
    }
  }
  free(dtmp);

  status = H5Sclose(dspace_id);

  status = H5Gclose(Obj_group_id);

  status = H5Fclose(file_id);

  status = H5close();

}
