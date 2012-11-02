#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mpi.h"


//-------------
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
//-------------



int head();
int ring();

//main function
int main(void) {
	
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
		ring();
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

	const int M=4;
    const int K=3;
    const int N=2;

	//result matrix C
	int C[M*N];
	int V[M];
	//initialize every cell in C to be 0
	for(int i=0;i<M*N;i++)
	{
		C[i]=0;
	}
	//for each column of C (N columns)
	for(int i=0;i<N;i++)
	{		
		//get the data from the other ring processes
		for(int j=1;j<size;j++)
		{
			MPI_Recv(V, M, MPI_INT, j, i, MPI_COMM_WORLD, &status);	
			//integrate the data from all the other ring processes into a complete vector
			for(int ij=0;ij<M;ij++)
			{
				C[ij*N+i]=C[ij*N+i]+V[ij];
			}
		}		
	}
	//print the results
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			cout<<C[i*N+j]<<"\t";
		}
		cout<<endl;
	}	
	return 0;
}
//the other ring processes
int ring()
{
	int rank; //
	int size; //
	MPI_Status status;
	
	//get the current process id
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//get the number of total processes
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	int const M(4),K(3),N(2);
	vector<int> indexes;	

	int V[M];
	
	int A[M*K]={1,2,3,4,5,6,7,8,9,8,7,6};
	int B[K*N]={5,4,3,2,1,0};
	 
	//compute which rows of matrix A should be distributed to the current process, save the row data to indexes
	//eg. 2 ring processes, A has 4 rows, then the process with id==1 get row 0 and 2
	for(int i=rank-1;i<M;i+=size-1)
	{
		indexes.push_back(i);
	} 
	//for each column of matrix B (N columns)
	for(int i=0;i<N;i++)
	{
		//initialize every cell in the vector V to be 0
		for(int j=0;j<M;j++)
		{
			V[j]=0;
		}
		if(indexes.size()>0)
		{						
			for(int j=0;j<indexes.size();j++)
			{
				for(int ij=0;ij<K;ij++)
				{
					V[indexes[j]]+=A[indexes[j]*K+ij]*B[ij*N+i];
				}
			}
			MPI_Send(V, M, MPI_INT, 0,i,MPI_COMM_WORLD);				
		}
		else
		{
			//send V with all 0 elements
			MPI_Send(V, M, MPI_INT, 0,i,MPI_COMM_WORLD);	
		}
	}	
	return 0;
}

