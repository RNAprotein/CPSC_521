/*
CPSC 521 assignment 2: matrix multiplication.

The original version is based on N-body script: do functional decomposition (pipeline: Initialization, simulation and termination)

The later version, according to Prof. Wagner's suggestion on Oct. 10th class: do node decomposition (root, other) which is better when adding new ring nodes

To build the program:
	mpicc -o hw2Mat hw2Mat.c
	or (create a log file)
	mpicc -o hw2Mat hw2Mat.c -mpe=mpilog

To run the program:
	mpiexec -n 5 ./hw2Mat 3 matrix.txt
	or (create a log file)
	mpiexec -n 5 ./hw2Mat 30 matrix.txt (31 rounds at the maximum due to the datatype (long int) limitation)
	clog2TOslog2 hw2Mat.clog2
	jumpshot hw2Mat.slog2


  (c) Shu Yang 2012
  Email: syang11@cs.ubc.ca

*/

#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

const int D=4;	//dimension; i.e. take 4*4 matrix as an example

int head(int round, long *M, char *fileName);
int ring(int round, long *M);

//main function
//argv: (rounds, filename)
int main(int argc,char *argv[]) {
    int test;
	int rank; //
	int size; //
	//MPI_Status status;
	
	//initialize MPI
	MPI_Init(&argc,&argv);
	//get the current process's rank, to rank
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the total number of processes, to size
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	//---------
	//exit if the number of arguments is not satisfied
	if (argc != 3) {
		printf("Please provide two arguments.\n %s rounds filename\n",argv[0]);
		MPI_Finalize();
		return 1;
	}
	//exit if the number of processes is not satisfied
	if(size!=D+1)
	{			
		printf("Need %d processes!\n",D+1);
		MPI_Finalize();
		return 1;
	}
	//---------
	
	int rounds = atoi(argv[1]);
	char *fileName=argv[2];

	long M[D*D];	//resultant matrix for each round
	long raw[D*D];	//raw matrix from the file
	int i;
	for (i = 0; i < rounds; i++) {
		if (rank == 0) {
			//set the head process, which is the one with id==0
			head(i, M, fileName); //head(fileName);
		} else {
			//set the other ring processes, which are those with id!=0
			ring(i, raw);
		}
	}
	
	//finalize MPI
	MPI_Finalize();
	return 1;
}


//for the head node
int head(int round, long *M, char *fileName)
{
	int rank; //
	int size; //
	MPI_Status status;
	
	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	//raw matrix from the file
	long raw[D*D];
    //int C[D*D];
    long V[D];
    long *p=M;

    //initialize the matrix raw and C
    if (round == 0) {
		FILE *fp;
		fp = fopen(fileName, "r");
		int i;
		for (i = 0; i < D * D; i++) {
			fscanf(fp, "%lu", p);
			raw[i]=*(p);
			p++;
		}
	}

	//send matrix to the first ring node
	int tag_code=round;//int tag_code=1;
	MPI_Send(M, D * D, MPI_LONG, 1, tag_code, MPI_COMM_WORLD);

	//receive matrix from the last ring node
	MPI_Recv(M, D * D, MPI_LONG, size-1, tag_code, MPI_COMM_WORLD, &status);

	//print the results
	printf("\n-------For the %d round:--------\n", round);
	int i;
	for(i=0;i<D;i++)
	{
		int j;
		for(j=0;j<D;j++)
		{
			printf ("%lu\t",M[i*D+j]);
		}
		printf ("\n");
	}
	return 0;
}

//the other ring processes
int ring(int round, long *raw)
{
	int rank; //
	int size; //
	MPI_Status status;
	long M[D*D];	//resultant matrix for each round
	long T[D*D];	//temporary resultant matrix for each round
	int tag_code=round;

	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	//find the neighbors
	int previous=rank-1;
	int next=rank+1;
	if(next==D+1){	//next==size
		next=0;
	}

	//raw matrix from the file
	if (round==0){
		if(rank==1){
			MPI_Recv(raw, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD, &status);
			int i;
			for (i = 0; i < D * D; i++) {
				M[i]=*(raw+i);
			}
		}else{
			MPI_Recv(raw, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD, &status);
			MPI_Recv(M, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD, &status);
		}

		int row=rank-1;
		int i;
		for(i=0;i<D;i++)
		{
			M[row*D+i]=0;
			int j;
			for(j=0;j<D;j++)
			{
				M[row*D+i]+=(*(raw+row*D+j)) * (*(raw+row*j+i));
			}
		}

		if(next !=0){
			MPI_Send(raw, D * D, MPI_LONG, next, tag_code, MPI_COMM_WORLD);
		}
		MPI_Send(M, D * D, MPI_LONG, next, tag_code, MPI_COMM_WORLD);
	}else{
		if (rank == 1) {
			MPI_Recv(M, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD,	&status);
			int i;
			for (i = 0; i < D * D; i++) {
				T[i] = M[i];
			}
		}else{
			MPI_Recv(T, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD, &status);
			MPI_Recv(M, D * D, MPI_LONG, previous, tag_code, MPI_COMM_WORLD, &status);
		}

		int row = rank - 1;
		int i;
		for (i = 0; i < D; i++) {
			M[row * D + i] = 0;
			int j;
			for (j = 0; j < D; j++) {
				M[row * D + i] += (T[row * D + j]) * (*(raw+row*j+i));
			}
		}

		if(next !=0){
			MPI_Send(T, D * D, MPI_LONG, next, tag_code, MPI_COMM_WORLD);
		}
		MPI_Send(M, D * D, MPI_LONG, next, tag_code, MPI_COMM_WORLD);
	}

	return 0;
}

