#include <stdio.h>
#include <omp.h>
int main()
{ int i;
 #pragma omp parallel for
 for (i=0;i<10;i++)
 printf("%d ",i);
 return 0;
}