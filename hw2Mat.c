/*
CPSC 521 assignment 2: matrix multiplication.

The original version is based on N-body script: do functional decomposition (pipeline: Initialization, simulation and termination)

The later version, according to Prof. Wagner's suggestion on Oct. 10th class: do node decomposition (root, other) which is better when adding new ring nodes

  (c) Shu Yang 2012
  Email: syang11@cs.ubc.ca

*/

#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

const int D=4;	//dimension; i.e. take 4*4 matrix as an example

int head(int round, int *M, char *fileName);
int ring(int round, int *M);

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
	//exit if the number of processes is not satisfied
	if(size<2)
	{			
		printf("Need at least 2 processes!\n");
		MPI_Abort(MPI_COMM_WORLD,99);
	}
	//---------
	
	int rounds = atoi(argv[1]);
	char *fileName=argv[2];

	int M[D*D];	//resultant matrix for each round
	int raw[D*D];	//raw matrix from the file
	for (int i = 0; i < rounds; i++) {
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
int head(int round, int *M, char *fileName)
{
	int rank; //
	int size; //
	MPI_Status status;
	
	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	//raw matrix from the file
	int raw[D*D];
    //int C[D*D];
    int V[D];
    int *p=M;

    //initialize the matrix raw and C
    if (round == 0) {
		FILE *fp;
		fp = fopen(fileName, "r");
		for (int i = 0; i < D * D; i++) {
			fscanf(fp, "%d", p);
			raw[i]=*(p);
			p++;
		}
	}

	//send matrix to the first ring node
	int tag_code=round;//int tag_code=1;
	MPI_Send(M, D * D, MPI_INT, 1, tag_code, MPI_COMM_WORLD);

	//receive matrix from the last ring node
	MPI_Recv(M, D * D, MPI_INT, size-1, tag_code, MPI_COMM_WORLD, &status);

	//print the results
	printf("\n-------For the %d round:--------\n", round);
	for(int i=0;i<D;i++)
	{
		for(int j=0;j<D;j++)
		{
			printf ("%d\t",M[i*D+j]);
		}
		printf ("\n");
	}
	return 0;
}

//the other ring processes
int ring(int round, int *raw)
{
	int rank; //
	int size; //
	MPI_Status status;
	int M[D*D];	//resultant matrix for each round
	int T[D*D];	//temporary resultant matrix for each round
	int tag_code=round;

	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	//find the neighbors
	int previous=rank-1;
	int next=rank+1;
	if(next==size){
		next=0;
	}

	//raw matrix from the file
	if (round==0){
		if(rank==1){
			MPI_Recv(raw, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD, &status);
			for (int i = 0; i < D * D; i++) {
				M[i]=*(raw+i);
			}
		}else{
			MPI_Recv(raw, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD, &status);
			MPI_Recv(M, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD, &status);
		}

		int row=rank-1;
		for(int i=0;i<D;i++)
		{
			M[row*D+i]=0;
			for(int j=0;j<D;j++)
			{
				M[row*D+i]+=(*(raw+row*D+j)) * (*(raw+row*j+i));
			}
		}

		if(next !=0){
			MPI_Send(raw, D * D, MPI_INT, next, tag_code, MPI_COMM_WORLD, &status);
		}
		MPI_Send(M, D * D, MPI_INT, next, tag_code, MPI_COMM_WORLD, &status);
	}else{
		if (rank == 1) {
			MPI_Recv(M, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD,	&status);
			for (int i = 0; i < D * D; i++) {
				T[i] = M[i];
			}
		}else{
			MPI_Recv(T, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD, &status);
			MPI_Recv(M, D * D, MPI_INT, previous, tag_code, MPI_COMM_WORLD, &status);
		}

		int row = rank - 1;
		for (int i = 0; i < D; i++) {
			M[row * D + i] = 0;
			for (int j = 0; j < D; j++) {
				M[row * D + i] += (T[row * D + j]) * (*(raw+row*j+i));
			}
		}

		if(next !=0){
			MPI_Send(T, D * D, MPI_INT, next, tag_code, MPI_COMM_WORLD, &status);
		}
		MPI_Send(M, D * D, MPI_INT, next, tag_code, MPI_COMM_WORLD, &status);
	}

	return 0;
}

