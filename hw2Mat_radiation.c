/*
CPSC 521 assignment 2: matrix multiplication.

My official version is hw2Mat.c which, similar to Nbody problem, uses a ring structure.

However, I feel the ring structure is not very "natural" in this case and I implement another version based on radiation structure with reference to a previous code from Blaise Barney.

Besides, here the matrix is not assumed to be symmetric so that any matrix A and B could be as input and multiply with each other.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define NRA 15                // number of rows in matrix A
#define NCA 10                 // number of columns in matrix A & number of rows in matrix B
#define NCB 5                  // number of columns in matrix B


int head();
int worker();

//main function
int main(int argc,char *argv[]) {
    int test;
	int rank; //
	int size; //
	MPI_Status status;
	
	//initialize MPI
	MPI_Init(NULL,NULL);
	//get the current process's rank, to rank
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the total number of processes, to size
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	//---------
	//exit if the number of processes is not satisfied
	if(size<2)
	{			
		printf("Need at least 2 processes!\n");
		MPI_Abort(MPI_COMM_WORLD,99);
	}
	//---------
	
	if(rank==0)
	{
		//set the head process, which is the one with id==0
		head();
	}
	else
	{
		//set the other ring processes, which are those with id!=0
		worker();
	}
	
	//finalize MPI
	MPI_Finalize();
	return EXIT_SUCCESS;
}


//for the head node
int head()
{
	int rank; //
	int size; //
	MPI_Status status;
	
	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int numworkers = size-1;

	double a[NRA][NCA], b[NCA][NCB], c[NRA][NCB];
	//initialize matrix to unit matrix
	for (i=0; i<NRA; i++)
	   for (j=0; j<NCA; j++)
		 a[i][j]= 1.0;

	for (i=0; i<NCA; i++)
	   for (j=0; j<NCB; j++)
		 b[i][j]= 1.0;

	 // send matrix data to the worker tasks
	 int averow = NRA/numworkers;
	 int extra = NRA%numworkers;
	 int offset = 0;
	 int mtype = 1;
	 int rows, i, j, k, rc, dest;
	 for (dest=1; dest<=numworkers; dest++)
	 {
	   rows = (dest <= extra) ? averow+1 : averow;
	   printf("Send %d rows to task %d\n",rows,dest);

	   //offset in matrix A
	   MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

	   //send the number of rows each process is required to compute
	   MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

	   //send each process rows*NCA bits of data starting at offset
	   MPI_Send(&a[offset][0], rows*NCA, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);


	   //send each process the matrix B
	   MPI_Send(&b, NCA*NCB, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
	   offset = offset + rows;
	 }

	 mtype = 2;
	 for (i=1; i<=numworkers; i++)
	 {
	   //recieve the offset value the sending process ended with
	   MPI_Recv(&offset, 1, MPI_INT, i, mtype, MPI_COMM_WORLD, &status);

	   //receive the number of rows the sending process computed
	   MPI_Recv(&rows, 1, MPI_INT, i, mtype, MPI_COMM_WORLD, &status);

	   //receive the final values of matrix C starting at the corresponding offset values                          */
	   MPI_Recv(&c[offset][0], rows*NCB, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
	 }

	 //print results
	 for (i=0; i<NRA; i++)
	 {
	   printf("\n");
	   for (j=0; j<NCB; j++)
		 printf("%f ", c[i][j]);
	 }
	 printf ("\n");

	return 0;
}

//the other worker processes
int worker()
{
	int rank; //
	int size; //
	MPI_Status status;
	
	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int mtype = 1;
	int rows, i, j, k, rc, dest;
	int offset = 0;

	double a[NRA][NCA], b[NCA][NCB], c[NRA][NCB];

	//receive the initial offset position of the matrix
	 MPI_Recv(&offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);

	 //receive the number of rows to compute/
	 MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);

	 //receive the matrix A starting at offset
	 MPI_Recv(&a, rows*NCA, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,                   &status);

	 //receive the matrix B
	 MPI_Recv(&b, NCA*NCB, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,                    &status);

	 for (k=0; k<NCB; k++)
	 for (i=0; i<rows; i++)
	 {
		 c[i][k] = 0.0;
		 for (j=0; j<NCA; j++)
		   c[i][k] = c[i][k] + a[i][j] * b[j][k];
	 }
	 mtype = 1;

	 //send the offset value
	 MPI_Send(&offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);

	 //send the number of rows
	 MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);

	 //send the final portion of C
	 MPI_Send(&c, rows*NCB, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);

	return 0;
}

