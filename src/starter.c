#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(int argc, char* argv[])
{

	/* Read the input file */

	int width, height, t;
	double p;
	//Read in file
	FILE* f=fopen("input.txt","w");
	assert(f!=NULL);
	scanf("%d %d %lf %d \n", &width, &height, &p, &t);
	int nbp=(int)sqrt(width);
	while(width%nbp!=0)nbp--;
	fprintf(f,"%d %d %lf %d \n", width, height, p, t);
	int cas, i,j;
	double value; 
	while((EOF != scanf("%d %d %d %lf", &cas, &i, &j, &value))) 
	{
		fprintf(f,"%d %d %d %lf", cas, i, j, value);
	}
	fclose(f);

	//showing
	char cmd[1024];
	sprintf(cmd,"mpirun -np %d constant input.txt",nbp);
	system(cmd);
	return 0;
}
