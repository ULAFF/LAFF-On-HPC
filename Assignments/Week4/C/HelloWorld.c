#include <stdio.h>
#include <stdlib.h>

#include "omp.h"

int main(int argc, char *argv[])
{
  
#pragma omp parallel 
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    
    printf( "Hello World! from %d of %d\n", tid, nthreads );
  }
}
