/*********************************************\
| Distributed algorithm to compute X^t+1 with |
|             grid of processors              |
\*********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

typedef struct { 
  int row; // number of lines
  int col; // number of columns
  double* data; // values
} matrix_t;

int my_id, nb_proc;      // the id of the processor, and the number of processors
int nb_col,nb_row;       // number of rows (and columns) of the grid
MPI_Comm MPI_HORIZONTAL; // communicator for horizontal broadcast
MPI_Comm MPI_VERTICAL;   // communicator for vertical broadcast
int i_col,i_row;         // position of the processor in the grid

/*********************\
| Read the input data |
\*********************/

void set_arguments(int widht, int height, int p, int t, matrix_t *mat, /* file */ )
{
 /* Simple parser */
}

/***************************\
| The distributed algorithm |
\***************************/

/* init communicators */

void init_communicators()
{
  nb_col = (int)sqrt(nb_proc);
  nb_row = nb_col; // assert square grid
  MPI_Comm_split(MPI_COMM_WORLD, my_id / nb_row, my_id % nb_col, &MPI_HORIZONTAL);
  MPI_Comm_split(MPI_COMM_WORLD, my_id % nb_col, my_id / nb_col, &MPI_VERTICAL);
}

/* average function */

double average(double center, double north, double south, double east, double west, double p)
{
  return (1-p)*center + (north + south + east + west) / 4;
}

/*************\
| The program |
\*************/

int main(int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  init_communicators();

  MPI_Finalize();

  return 0;
}


