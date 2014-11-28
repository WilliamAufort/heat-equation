/*********************************************\
| Distributed algorithm to compute X^t+1 with |
|             grid of processors              |
\*********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

/*************************************************\
| One processor works on several buffers :
| - one matrix which contains the initial
| values of X which are not on the border;
| - 4 "vectors" which corresponds to the 2 
| lines and 2 columns at the border (used to
| send and to do the computations;
| - 4 other vectors which are used for 
| receving data from the neighbourgs.
|
| We use these vectors instead of global matrix,
| because we don't need to copy data (especiallty 
| before SEND) during the execution.
\*************************************************/

int my_id, nb_proc;      // the id of the processor, and the number of processors
int nb_col,nb_row;       // number of rows (and columns) of the grid
MPI_Comm MPI_HORIZONTAL; // communicator for horizontal broadcast
MPI_Comm MPI_VERTICAL;   // communicator for vertical broadcast
int i_col,i_row;         // position of the processor in the grid

/* Processors data */

double* matrix, first_row, last_row, first_col, last_col;

/* initialize data */

void init_datas()
{
	nb_col = (int)sqrt(nb_proc);
	nb_row = nb_col; // assert square grid
	
	matrix = malloc(sizeof(double)*nb_col*nb_row);
	first_row = malloc(sizeof(double)*nb_col);
	last_row = malloc(sizeof(double)*nb_col);
	first_col = malloc(sizeof(double)*nb_row);
	last_col = malloc(sizeof(double)*nb_row);
}

/* free data */

void free_datas() 
{
	free(first_row);
	free(last_row);
	free(first_col);
	free(last_col);
}

/*********************\
| Read the input data |
\*********************/

void set_arguments(int widht, int height, int p, int t, matrix_t *mat, /* file */ )
{
 /* Simple parser TODO */
}

/***************************\
| The distributed algorithm |
\***************************/

/* init communicators */

void init_communicators()
{
  MPI_Comm_split(MPI_COMM_WORLD, my_id / nb_col, my_id % nb_col, &MPI_HORIZONTAL);
  MPI_Comm_split(MPI_COMM_WORLD, my_id % nb_col, my_id / nb_col, &MPI_VERTICAL);
}

/* average function */

double average(double center, double north, double south, double east, double west, double p)
{
  return (1-p)*center + (north + south + east + west) / 4;
}

/* the core of the parallel algorithm */

void compute_image() 
{
	// usefull variables
	int my_col = my_id % nb_col;
	int my_row = my_id / nb_col;

	// init buffers for receiving
	double* tmp_first_row = malloc(sizeof(double)*/*size*/);
	double* tmp_last_row = malloc(sizeof(double)*/*size*/);
	double* tmp_first_col = malloc(sizeof(double)*/*size*/);
	double* tmp_last_col = malloc(sizeof(double)*/*size*/);

	// Scatter data

	// Do the sends
    MPI_Send(first_row, nb_col, MPI_DOUBLE, (my_row - 1) % nb_row, 1, MPI_VERTICAL);
    MPI_Send(last_row, nb_col, MPI_DOUBLE, (my_row + 1) % nb_row, 1, MPI_VERTICAL);

    MPI_Send(first_col, nb_row, MPI_DOUBLE, (my_col - 1) % nb_col, 1, MPI_HORIZONTAL);
    MPI_Send(last_col, nb_row, MPI_DOUBLE, (my_col + 1) % nb_col, 1, MPI_VERTICAL);
    
    // Do the computations
	int i, j;
	for(i = 0; i < /*size*/, i++) {
		for(i = 0; j < /*size*/, j++) {
			// update TODO
        }
	}

    // Do the receives 
    MPI_Recv(tmp_first_row, nb_col, MPI_DOUBLE, (my_row + 1) % nb_row, 1, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(tmp_last_row, nb_col, MPI_DOUBLE, (my_row - 1) % nb_row, 1, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(tmp_first_col, nb_row, MPI_DOUBLE, (my_col + 1) % nb_col, 1, MPI_HORIZONTAL, MPI_STATUS_IGNORE);
	MPI_Recv(tmp_last_col, nb_row, MPI_DOUBLE, (my_row - 1) % nb_col, 1, MPI_HORIZONTAL, MPI_STATUS_IGNORE);

	// Do the lasts computations TODO

	// Free the buffers
 
	free(tmp_first_row);
	free(tmp_first_col);
	free(tmp_last_row);
	free(tmp_last_col);

	// Reconstruct the matrix
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
  init_datas();
  // init_matrix(); TODO

  compute_image();

  MPI_Finalize();

  return 0;
}


