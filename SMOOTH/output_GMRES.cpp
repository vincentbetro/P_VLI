#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"
#include "smooth.h"

void Smesh_obj::output_GMRES(int mdim, double K[][3][3], int ia[], int iau[], int ja[],
                             double u[], double v[], double w[])
{
  FILE *fp;
  const int bdim = 80;
  char filename[bdim];
  int i, j, n;

  // create file name
  filename[0]='\0';
  strcat(filename,"gmres.out");

  if ((fp = fopen(filename,"w")) == NULL)
  {
    printf("\nError opening file <%s>.",filename);
    exit(0);
  }

  fprintf(fp,"# Dimension of matrix\n");
  fprintf(fp,"%d\n",nn);
  fprintf(fp,"# Number of non-zero matrices\n");
  fprintf(fp,"%d\n",mdim);

  fprintf(fp,"# ia matrix\n");
  for (n=0; n < nn+1; n++)
    fprintf(fp,"%d\n",ia[n]);

  fprintf(fp,"# iau matrix\n");
  for (n=0; n < nn; n++)
    fprintf(fp,"%d\n",iau[n]);

  fprintf(fp,"# ja matrix\n");
  for (n=0; n < mdim; n++)
    fprintf(fp,"%d\n",ja[n]);

  fprintf(fp,"# A matrix\n");
  for (n=0; n < mdim; n++)
  {
    for (i=0; i < 3; i++)
    {
      for (j=0; j < 3; j++)
        fprintf(fp,"%22.15e ",K[n][i][j]);
      fprintf(fp,"\n");
    }
  }

  fprintf(fp,"# b matrix\n");
  for (n=0; n < nn; n++)
    fprintf(fp,"%22.15e %22.15e %22.15e\n",u[n],v[n],w[n]);

  fprintf(fp,"# x matrix\n");
  for (n=0; n < nn; n++)
    fprintf(fp,"0.0 0.0 0.0\n");

  fclose(fp);

}
