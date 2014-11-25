/*********************************************\
| Distributed algorithm to compute X^t+1 with |
|             grid of processors              |
\*********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

typedef struct { 
  int row; // number of lines
  int col; // number of columns
  double* data; // values
} matrix_t;

/*********************\
| Read the input data |
\*********************/



void set_arguments(int widht, int height, int p, int t, matrix_t *mat, /* file */ ) {

}

/***************************\
| The distributed algorithm |
\***************************/

/* average function */

double average(double center, double north, double south, double east, double west, double p) {
  return (1-p)*center + (north + south + east + west) / 4;
}
