#include <stdio.h>

#ifndef Bar_2_h
#define Bar_2_h
class Bar_2
{
 public:

  void init()
  {
    for (int i=0; i < 2; i++)
      nodes[i]= -1;
  }
  void print(FILE *fp)
  {
    fprintf(fp,"\nEdge nodes:");
    int i;
    for (i=0; i < 2; i++)
      fprintf(fp,"\nNode %d = %d",i,nodes[i]);
  }

  int nodes[2];
};
#endif

#ifndef Bar_3_h
#define Bar_3_h
class Bar_3
{
 public:

  void init()
  {
    for (int i=0; i < 3; i++)
      nodes[i]= -1;
  }
  void print(FILE *fp)
  {
    fprintf(fp,"\nEdge nodes:");
    int i;
    for (i=0; i < 3; i++)
      fprintf(fp,"\nNode %d = %d",i,nodes[i]);
  }

  int nodes[3];
};
#endif
