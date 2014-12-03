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
double* work_matrix, neighbor_first_row, neighbor_last_row, neighbor_first_col, neighbor_last_col;

/* initialize data */

void init_datas()
{
	nb_col = (int)sqrt(nb_proc);
	nb_row = nb_col; // assert square grid
	
	matrix = malloc(sizeof(double)*(nb_col-2)*(nb_row-2));
	work_matrix = malloc(sizeof(double)*(nb_col-2)*(nb_row-2));
	first_row = malloc(sizeof(double)*nb_col);
	last_row = malloc(sizeof(double)*nb_col);
	first_col = malloc(sizeof(double)*(nb_row));
	last_col = malloc(sizeof(double)*(nb_row));
	neighbor_first_row = malloc(sizeof(double)*nb_col);
	neighbor_last_row = malloc(sizeof(double)*nb_col);
	neighbor_first_col = malloc(sizeof(double)*(nb_row));
	neighbor_last_col = malloc(sizeof(double)*(nb_row));
}

/* free data */

void free_datas() 
{
	free(matrix);
	free(first_row);
	free(last_row);
	free(first_col);
	free(last_col);
	free(work_matrix);
	free(neighbor_first_row);
	free(neighbor_last_row);
	free(neighbor_first_col);
	free(neighbor_last_col);
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
  return (1.0-p)*center + p*(north + south + east + west) / 4.0;
}

/* the core of the parallel algorithm */

void compute_image(double p) 
{
	// usefull variables
	int my_col = my_id % nb_col;
	int my_row = my_id / nb_col;
	int my_col_mid = my_col-2;
	int my_row_mid = my_row-2;

	// Scatter data

	// Do the sends
	int west,east,north,south;
	north=(my_row + 1) % nb_row;
	south=(my_row - 1 + nb_row) % nb_row;//ensure index >0
	east=(my_col - 1 + nb_col) % nb_col;
	west=(my_col + 1) % nb_col;
	MPI_Send(first_row, nb_col, MPI_DOUBLE, norh, 1, MPI_VERTICAL);
	MPI_Send(last_row, nb_col, MPI_DOUBLE, south, 1, MPI_VERTICAL);
	MPI_Send(first_col, nb_row, MPI_DOUBLE, west, 1, MPI_HORIZONTAL);
	MPI_Send(last_col, nb_row, MPI_DOUBLE, east, 1, MPI_HORIZONTAL);
    
    // Do the computations
	int i, j;
	//middle
	for(i = 1; i < nb_col_mid-1, i++) {
		for(j = 1; j < nb_row_mid-1, j++) {
			// update TODO
			work_matrix[i+nb_col*j]=average(matrix[i+nb_col*j],matrix[i+nb_col*(j-1)],matrix[i+nb_col*(j+1)],matrix[i+1+nb_col*j]+matrix[i-1+nb_col*j]);	
        	}
	}
	work_matrix[0]=average(
		matrix[0],
		first_row[1],
		matrix[nb_col_mid],
		matrix[1],
		first_col[1]
	);
	work_matrix[nb_col_mid-1]=average(
		matrix[nb_col_mid-1],
		first_row[nb_col_mid],
		matrix[nb_col_mid+nb_col_mid-1],
		last_col[1],
		matrix[nb_col_mid-2]
	);
	work_matrix[nb_col_mid*(nb_row_mid-1)]=average(
		matrix[nb_col_mid*(nb_row_mid-1)],
		matrix[nb_col_mid*(nb_row_mid-1)-nb_col_mid],
		last_row[1],
		matrix[nb_col*(nb_row_mid-1)+1],
		first_col[nb_row_mid-2]
	);
	[nb_col_mid*nb_row_mid-1]=average(
		matrix[nb_col_mid*nb_row_mid-1],
		matrix[nb_col_mid*nb_row_mid-1-nb_col_mid],
		last_row[nb_col_mid],
		last_col[nb_row_mid],
		matrix[nb_col_mid*nb_row_mid-2]
	);

	//lines
	for(i=1;i < nb_col_mid-1;i++)
	{
		work_matrix[i]=average(
			matrix[i],
			first_row[i+1],
			matrix[i+nb_col_mid],
			matrix[i+1],
			matrix[i-1]
		);
		work_matrix[nb_col_mid*(nb_row_mid-1)+i]=average(
			matrix[nb_col_mid*(nb_row_mid-1)+i],
			matrix[nb_col_mid*(nb_row_mid-2)+i],
			last_row[i+1],
			matrix[nb_col_mid*(nb_row_mid-1)+i+1],
			matrix[nb_col_mid*(nb_row_mid-1)+i-1]
		);
	}
	//cols
	for(j=1;i < nb_row_mid-1;j++)
	{
		work_matrix[nb_col_mid*j]=average(
			matrix[nb_col_mid*j],
			matrix[nb_col_mid*(j-1)],
			matrix[nb_col_mid*(j+1)],
			matrix[nb_col_mid*j+1],
			first_col[j+1]
		);
		work_matrix[nb_col_mid*(j+1)-1]=average(
			matrix[nb_col_mid*(j+1)-1],
			matrix[nb_col_mid*j-1],
			matrix[nb_col_mid*(j+2)-1],
			last_col[j+1],
			matrix[nb_col_mid*(j+1)-2]
		);
	}


	// Do the receives 
	MPI_Recv(neighbor_first_row, nb_col, MPI_DOUBLE, north, 1, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_last_row, nb_col, MPI_DOUBLE, south, 1, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_first_col, nb_row, MPI_DOUBLE, south, 1, MPI_HORIZONTAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_last_col, nb_row, MPI_DOUBLE, north, 1, MPI_HORIZONTAL, MPI_STATUS_IGNORE);

	// Do the lasts computations TODO



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


