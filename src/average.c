/*********************************************\
| Distributed algorithm to compute X^t+1 with |
|             grid of processors              |
\*********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <math.h>
#include "gfx.h"

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
int nb_col,nb_row;       // number of rows (and columns) of the subgrid handeld by the proc
int size;				 // size of the grid (NxN where the grid has size NxN)
MPI_Comm MPI_HORIZONTAL; // communicator for horizontal broadcast
MPI_Comm MPI_VERTICAL;   // communicator for vertical broadcast
int my_col,my_row;         // position of the processor in the grid


/* Processors data */

double* matrix;
double* first_row; 
double* last_row;
double* first_col;
double* last_col;
double* work_matrix;
double* work_first_row; 
double* work_last_row;
double* work_first_col;
double* work_last_col;
double* neighbor_first_row;
double* neighbor_last_row;
double* neighbor_first_col;
double* neighbor_last_col;

/* initialize data */

void init_datas()
{
	nb_col = (int)sqrt(size / nb_proc); // assume nb_proc is a square and nb_proc divides NxN
	nb_row = nb_col; // assert square grid
	int number_proc = (int)sqrt(nb_proc);
	my_col = my_id % number_proc;
	my_row = my_id / number_proc;
	
	matrix = calloc((nb_col-2)*(nb_row-2), sizeof(double));
	work_matrix = calloc((nb_col-2)*(nb_row-2),sizeof(double));
	first_row = calloc(nb_col, sizeof(double));
	last_row = calloc(nb_col, sizeof(double));
	first_col = calloc(nb_row, sizeof(double));
	last_col = calloc(nb_row, sizeof(double));
	neighbor_first_row = calloc(nb_col,sizeof(double));
	neighbor_last_row = calloc(nb_col,sizeof(double));
	neighbor_first_col = calloc(nb_row,sizeof(double));
	neighbor_last_col = calloc(nb_row,sizeof(double));
	work_first_row = calloc(nb_col,sizeof(double));
	work_last_row = calloc(nb_col,sizeof(double));
	work_first_col = calloc(nb_row,sizeof(double));
	work_last_col = calloc(nb_row,sizeof(double));

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
	free(work_first_row);
	free(work_last_row);
	free(work_first_col);
	free(work_last_col);
	free(neighbor_first_row);
	free(neighbor_last_row);
	free(neighbor_first_col);
	free(neighbor_last_col);
}

/***************************\
| The distributed algorithm |
\***************************/

/* init communicators */

void init_communicators()
{
	int p = (int)sqrt(nb_proc);
	MPI_Comm_split(MPI_COMM_WORLD, my_id / p, my_id % p, &MPI_HORIZONTAL);
	MPI_Comm_split(MPI_COMM_WORLD, my_id % p, my_id / p, &MPI_VERTICAL);
}

void swap(double** l1,double** l2)
{
	double* l3=*l1;
	*l1=*l2;
	*l2=l3;
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
	int nb_col_mid = nb_col-2;
	int nb_row_mid = nb_row-2;

	// Scatter data
	// Do the sends
	int N=(int)sqrt(nb_proc);
	int west,east,north,south;
	north=(my_row + 1) % N;
	south=(my_row - 1 + N) % N;
	east=(my_col - 1 + N) % N;
	west=(my_col + 1) % N;
	MPI_Send(first_row, nb_col, MPI_DOUBLE, north, 1, MPI_VERTICAL);
	MPI_Send(last_row, nb_col, MPI_DOUBLE, south, 2, MPI_VERTICAL);
	MPI_Send(first_col, nb_row, MPI_DOUBLE, west, 3, MPI_HORIZONTAL);
	MPI_Send(last_col, nb_row, MPI_DOUBLE, east, 4, MPI_HORIZONTAL);
    
    /* Do the computations */

	int i, j;
	// middle values
	for(i = 1; i < nb_col_mid-1; i++) {
		for(j = 1; j < nb_row_mid-1; j++) {
			work_matrix[i+nb_col_mid*j]=average(
				matrix[i+nb_col_mid*j],
				matrix[i+nb_col_mid*(j-1)],
				matrix[i+nb_col_mid*(j+1)],
				matrix[i+1+nb_col_mid*j],
				matrix[i-1+nb_col_mid*j],
				p);	
        	}
	}
	// Corners of work_matrix
	work_matrix[0]=average( // ok
		matrix[0],
		first_row[1],
		matrix[nb_col_mid],
		matrix[1],
		first_col[1],
		p);
	work_matrix[nb_col_mid-1]=average( // ok
		matrix[nb_col_mid-1],
		first_row[nb_col_mid],
		matrix[nb_col_mid+nb_col_mid-1],
		last_col[1],
		matrix[nb_col_mid-2],
		p);
	work_matrix[nb_col_mid*(nb_row_mid-1)]=average( // ok
		matrix[nb_col_mid*(nb_row_mid-1)],
		matrix[nb_col_mid*(nb_row_mid-1)-nb_col_mid],
		last_row[1],
		matrix[nb_col_mid*(nb_row_mid-1)+1],
		first_col[nb_row_mid],
		p);
	work_matrix[nb_col_mid*nb_row_mid-1]=average( // ok
		matrix[nb_col_mid*nb_row_mid-1],
		matrix[nb_col_mid*nb_row_mid-1-nb_col_mid],
		last_row[nb_col_mid],
		last_col[nb_row_mid],
		matrix[nb_col_mid*nb_row_mid-2],
		p);
	// first and last lines of work_matrix
	for(i=1;i < nb_col_mid-1;i++)
	{
		work_matrix[i]=average(
			matrix[i],
			first_row[i+1],
			matrix[i+nb_col_mid],
			matrix[i+1],
			matrix[i-1],
			p);
		work_matrix[nb_col_mid*(nb_row_mid-1)+i]=average(
			matrix[nb_col_mid*(nb_row_mid-1)+i],
			matrix[nb_col_mid*(nb_row_mid-2)+i],
			last_row[i+1],
			matrix[nb_col_mid*(nb_row_mid-1)+i+1],
			matrix[nb_col_mid*(nb_row_mid-1)+i-1],
			p);
	}

	//first and last columns of work_matrix
	for(j=1;i < nb_row_mid-1;j++)
	{
		work_matrix[nb_col_mid*j]=average(
			matrix[nb_col_mid*j],
			matrix[nb_col_mid*(j-1)],
			matrix[nb_col_mid*(j+1)],
			matrix[nb_col_mid*j+1],
			first_col[j+1],
			p);
		work_matrix[nb_col_mid*(j+1)-1]=average(
			matrix[nb_col_mid*(j+1)-1],
			matrix[nb_col_mid*j-1],
			matrix[nb_col_mid*(j+2)-1],
			last_col[j+1],
			matrix[nb_col_mid*(j+1)-2],
			p);
	}
	// Receive datas for neighbors process
	MPI_Recv(neighbor_first_row, nb_col, MPI_DOUBLE, north, 1, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_last_row, nb_col, MPI_DOUBLE, south, 2, MPI_VERTICAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_first_col, nb_row, MPI_DOUBLE, west, 3, MPI_HORIZONTAL, MPI_STATUS_IGNORE);
	MPI_Recv(neighbor_last_col, nb_row, MPI_DOUBLE, east, 4, MPI_HORIZONTAL, MPI_STATUS_IGNORE);

	/* Do the lasts computations */

	// first and last rows
	for(i=1;i < nb_col-1;i++)
	{
		work_first_row[i]=average(
			first_row[i],
			neighbor_last_row[i],
			matrix[i-1],
			first_row[i+1],
			first_row[i-1],
			p);
		work_last_row[i]=average(
			last_row[i],
			matrix[nb_col_mid*(nb_row_mid-1)+i-1],
			neighbor_first_row[i],
			last_row[i+1],
			last_row[i-1],
			p);
	}

	// first and last columns
	for(j=1;j < nb_row-1;j++)
	{
		work_first_col[j]=average(
			first_col[j],
			first_col[j+1],
			first_col[j-1],
			matrix[nb_col_mid*(j-1)],
			neighbor_last_col[j],
			p);
		work_last_col[j]=average(
			last_col[j],
			last_col[j+1],
			last_col[j-1],
			neighbor_first_col[j],
			matrix[nb_col_mid*j-1],
			p);
	}
	// The 4 last corners
	work_first_row[0]=average(first_row[0],
		neighbor_last_row[0],
		first_col[1],
		first_row[1],
		neighbor_last_col[0],
		p);
	work_first_col[0]=work_first_row[0];

	work_first_row[nb_col-1]=average(first_row[nb_col-1],
		neighbor_last_row[nb_col-1],
		last_col[1],
		neighbor_first_col[0],
		first_row[nb_col-2],
		p);
	work_last_col[0]=work_first_row[nb_col-1];

	work_first_col[nb_col-1]=average(first_col[nb_col-1],
		first_col[nb_row-2],
		neighbor_first_row[0],
		last_row[1],
		neighbor_last_col[nb_row-1],
		p);
	work_last_row[0]=work_first_col[nb_col-1];

	work_last_col[nb_row-1]=average(last_col[nb_row-1],
		last_col[nb_row-2],
		neighbor_first_row[nb_col-1],
		neighbor_first_col[nb_row-1],
		last_row[nb_col-2],
		p);
	work_last_row[nb_col-1]=work_last_col[nb_row-1];
	
	/* Reconstruct the matrix (swaping) */

	swap(&work_matrix,&matrix);
	swap(&work_first_row,&first_row);
	swap(&work_last_row,&last_row);
	swap(&work_first_col,&first_col);
	swap(&work_last_col,&last_col);
}

//set data
void setheat(int i,int j,double t)
{
	i-=my_col*nb_col;
	j-=my_row*nb_row;

	if((0<=i)&&(i<nb_col)
	&& (0<=j)&&(j<nb_row))
	{
		//center
		if((0<i)&&(i<nb_col-1)
		&& (0<j)&&(j<nb_row-1))
		{
			i--;
			j--;
			matrix[i+(nb_col-2)*j]=t;
		}
		//border (no else for corners)
		else 
		{
			if(j==0)
				first_row[i]=t;
			if(j==nb_row-1)
				last_row[i]=t;
			if(i==0)
				first_col[j]=t;
			if(i==nb_col-1)
				last_col[j]=t;
		}
	}
	
}

void printheat(int i,int j)
{
	i-=my_col*nb_col;
	j-=my_row*nb_row;

	if((0<=i)&&(i<nb_col)
	&& (0<=j)&&(j<nb_row))
	{
		printf("temperature of case %i,%i is ",i+my_col*nb_col,j+my_row*nb_row);
		//center
		if((0<i)&&(i<nb_col-1)
		&& (0<j)&&(j<nb_row-1))
		{
			i--;
			j--;
			printf("%e",matrix[i+(nb_col-2)*j]);
		}
		//border
		else 
		{
			if(j==0)
				printf("%e",first_row[i]);
			else if(j==nb_row-1)
				printf("%e",last_row[i]);
			else if(i==0)
				printf("%e",first_col[j]);
			else if(i==nb_col-1)
				printf("%e",last_col[j]);
		}
		printf(".\n");
	}
	
}

/* Divisibility test */

int divides(int a, int b) {
	return ((b / a) * a == b);
}

//drawing
double Gamma = 0.80;
double IntensityMax = 255;

/** Taken from Earl F. Glynn's web page:
* <a href="http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm">Spectra Lab Report</a>
* */
void waveLengthToRGB(double Wavelength){
	double factor;
	double Red,Green,Blue;

	if((Wavelength >= 380) && (Wavelength<440)){
		Red = -(Wavelength - 440) / (440 - 380);
		Green = 0.0;
		Blue = 1.0;
	}else if((Wavelength >= 440) && (Wavelength<490)){
		Red = 0.0;
		Green = (Wavelength - 440) / (490 - 440);
		Blue = 1.0;
	}else if((Wavelength >= 490) && (Wavelength<510)){
		Red = 0.0;
		Green = 1.0;
		Blue = -(Wavelength - 510) / (510 - 490);
	}else if((Wavelength >= 510) && (Wavelength<580)){
		Red = (Wavelength - 510) / (580 - 510);
		Green = 1.0;
		Blue = 0.0;
	}else if((Wavelength >= 580) && (Wavelength<645)){
		Red = 1.0;
		Green = -(Wavelength - 645) / (645 - 580);
		Blue = 0.0;
	}else if((Wavelength >= 645) && (Wavelength<781)){
		Red = 1.0;
		Green = 0.0;
		Blue = 0.0;
	}else{
		Red = 0.0;
		Green = 0.0;
		Blue = 0.0;
	};

	// Let the intensity fall off near the vision limits

	if((Wavelength >= 380) && (Wavelength<420)){
		factor = 0.3 + 0.7*(Wavelength - 380) / (420 - 380);
	}else if((Wavelength >= 420) && (Wavelength<701)){
		factor = 1.0;
	}else if((Wavelength >= 701) && (Wavelength<781)){
		factor = 0.3 + 0.7*(780 - Wavelength) / (780 - 700);
	}else{
		factor = 0.0;
	};
	gfx_color(
		Red==0.0 ? 0 : (int) (IntensityMax * pow(Red * factor, Gamma)),
		Green==0.0 ? 0 : (int) (IntensityMax * pow(Green * factor, Gamma)),
		Blue==0.0 ? 0 : (int) (IntensityMax * pow(Blue * factor, Gamma)));

}
//380-781

void setcolor(double x,int i,int j)
{
	waveLengthToRGB(380+x*(781-380));
	gfx_point(i,j);
}

void drawbloc(double* data,int bloci,int blocj)
{
	bloci*=nb_col;
	blocj*=nb_row;
	int i,j;
	for(i=1;i<nb_col-1;i++)
		for(j=1;j<nb_row-1;j++)
			setcolor(data[(j-1)*(nb_col-2)+(i-1)],bloci+i,blocj+j);	
}

void drawfirstrow(double* data,int bloci,int blocj)
{
	bloci*=nb_col;
	blocj*=nb_row;
	int i;
	for(i=0;i<nb_col;i++)
		setcolor(data[i],bloci+i,blocj);	
}

void drawlastrow(double* data,int bloci,int blocj)
{
	bloci*=nb_col;
	blocj++;
	blocj*=nb_row;
	blocj--;
	int i;
	for(i=0;i<nb_col;i++)
		setcolor(data[i],bloci+i,blocj);	
}

void drawfirstcol(double* data,int bloci,int blocj)
{
	bloci*=nb_col;
	blocj*=nb_row;
	int i;
	for(i=0;i<nb_col;i++)
		setcolor(data[i],bloci,blocj+i);	
}

void drawlastcol(double* data,int bloci,int blocj)
{
	bloci++;
	bloci*=nb_col;
	blocj*=nb_row;
	bloci--;
	int i;
	for(i=0;i<nb_col;i++)
		setcolor(data[i],bloci,blocj+i);	
}

void senddata()
{
	MPI_Send(first_row, nb_col, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
	MPI_Send(last_row, nb_col, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
	MPI_Send(first_col, nb_row, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
	MPI_Send(last_col, nb_row, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
	MPI_Send(matrix, (nb_row-2)*(nb_col-2), MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
}

void receiveandshow()
{
	MPI_Status status;
	int c,r,i;
	int N=(int)sqrt(nb_proc);
	for(i=0;i<nb_proc;i++)
	{
		MPI_Recv(work_first_row, nb_col, MPI_DOUBLE, i, 5, MPI_COMM_WORLD,&status);
		MPI_Recv(work_last_row, nb_col, MPI_DOUBLE, i, 6, MPI_COMM_WORLD,&status);
		MPI_Recv(work_first_col, nb_row, MPI_DOUBLE, i, 7, MPI_COMM_WORLD,&status);
		MPI_Recv(work_last_col, nb_row, MPI_DOUBLE, i, 8, MPI_COMM_WORLD,&status);
		MPI_Recv(work_matrix, (nb_row-2)*(nb_col-2), MPI_DOUBLE, i, 9, MPI_COMM_WORLD,&status);
		c = i % N;
		r = i / N;
		drawbloc(work_matrix,c,r);
		drawfirstcol(work_last_col,c,r);
		drawlastcol(work_first_col,c,r);
		drawfirstrow(work_first_row,c,r);
		drawlastrow(work_last_row,c,r);
	}		
}


/*************\
| The program |
\*************/

int main(int argc, char* argv[])
{
	nb_proc=0;
	my_id=-1;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);


	/* Read the input file */

	int width, height, t;
	double p;
	//Read in file
	assert(argc>=2);
	FILE* f=fopen(argv[1],"r");
	assert(f!=NULL);
	fscanf(f,"%d %d %lf %d \n", &width, &height, &p, &t);
	size=width*height;
	init_datas();
	init_communicators();

  	assert(width == height); // The grid have to be a square grid
	assert(size==nb_row*nb_col*nb_proc); // We consider that the number of processors is a square"
	assert(divides(nb_proc,size)); 
	if (my_id == 0) {
		printf("%d %d %lf %d \n", width, height, p, t);
	}
	int cas, i, j;
	double value; 
	int stop=1;
	while(stop && (EOF != fscanf(f,"%d %d %d %lf", &cas, &i, &j, &value))) {
		switch (cas) { 
		case 0: 
			setheat(i,j,value);
			break;
		case 1:
			fprintf(stderr, "Error : we don't consider constants here \n");
			free_datas();
			exit(1);
			break;
		case 2: // keep request and stop
			stop=0;
			break;
		default:
			fprintf(stderr, "Error : incorrect option : %d \n", cas);
			free_datas();
			break;
		}
	}
	int step;
	for(step=0;step<t;step++)
	{
		compute_image(p);
		//printf("(%d,%d) step %d\n",my_col,my_row,step);
	}
	if(cas==2)
	{
		printheat(i,j);
		while((EOF != fscanf(f,"%d %d %d %lf", &cas, &i, &j, &value))) {
			switch (cas) { 
			case 0: 
			
				break;
			case 1:
			
				break;
			case 2: 
				printheat(i,j);
				break;
			default:
				fprintf(stderr, "Error : incorrect option : %d \n", cas);
				free_datas();
				break;
			}
		}
	}
	fclose(f);

	//showing
	senddata();
	if(my_id==0)
	{
		gfx_open(width,height,"Output");
		receiveandshow();
	}

	free_datas();
	MPI_Finalize();	
	if(my_id==0)
	{
		char c;
		while(1) 
		{
			// Wait for the user to press a character.
			c = gfx_wait();

			// Quit if it is the letter q.
			if(c=='q') break;
		}
	}

	return 0;
}


