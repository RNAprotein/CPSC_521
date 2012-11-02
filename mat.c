#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mpi.h"

#define tag_code 3333
static int node = 4;

void cal_mat(int *data, int *raw, int *temp, int size, int rank, int gran);
void get_row(int *data, int *row, int num_row);
void get_col(int *data, int *col, int num_col);
int main(int argc,char *argv[])
{
	int size, rank, i, j, k, right, left, rounds,gran;

	rounds = atoi(argv[1]);

	char* file_name;
	//initilize the MPI environment
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	gran = atoi(argv[2]);
	file_name = argv[3];

	MPI_Status recv_status;
	int data[node * node];
	int raw[node*node];
	int temp[node * node];
	right=rank+1;
	if(right == size) right=0;
	left=rank-1;
	if(left == -1) left = size-1;



	FILE *fp;
	int *p = data;

	//Initialization Part

	if(rank == 0){
		//Store each value in data[]

		fp = fopen(file_name, "r");
		for (i = 0; i < node * node; i++)
		{
			fscanf(fp, "%d", p);
			p++;
		}
		p = data;


		MPI_Send(p, node * node, MPI_INT, right, tag_code, MPI_COMM_WORLD);
	}
	else if(rank == size -1){
		MPI_Recv(data, node * node, MPI_INT, left, tag_code, MPI_COMM_WORLD, &recv_status);
	}
	else{
		MPI_Recv(data, node * node, MPI_INT, left, tag_code, MPI_COMM_WORLD, &recv_status);
		MPI_Send(p, node * node, MPI_INT, right, tag_code, MPI_COMM_WORLD);

	}

	for(k = 0; k < node * node; k++){
		// printf("%d\t", data[k]);
		raw[k] = data[k];
	}

	//Simulation Part

	for(j = 0; j < rounds; j++){

		cal_mat(data, raw, temp, size, rank, gran);
		//printf("%d %d", rounds, gran * node);
		// for(i = 0; i < gran * node; i++){
		//
		// printf("%d ", temp[i]);
		// }
		p = temp;
		for(i = 0; i < size-1; i++){

			MPI_Send(p, node * gran, MPI_INT, left, tag_code, MPI_COMM_WORLD);
			p += node * gran;
			MPI_Recv(p, node * gran, MPI_INT, right, tag_code, MPI_COMM_WORLD, &recv_status);
		}

		if(rank == 0){
			for(i = 0; i< node * node; i++){
				data[i] = temp[i];
			}
			MPI_Send(p, node * node, MPI_INT, right, tag_code, MPI_COMM_WORLD);
		}
		else if(rank == size -1){
			MPI_Recv(data, node * node, MPI_INT, left, tag_code, MPI_COMM_WORLD, &recv_status);
		}
		else{
			MPI_Recv(data, node * node, MPI_INT, left, tag_code, MPI_COMM_WORLD, &recv_status);
			MPI_Send(p, node * node, MPI_INT, right, tag_code, MPI_COMM_WORLD);

		}
	}


	p = data;
	if(rank == 0){
		for(i = 0; i < node; i++){
			for (j = 0; j < node; j++)
				printf("%d\t", *(p + i*node + j));
			printf("\n");
		}
	}

	MPI_FINALIZE();

}


void cal_mat(int *data, int *raw, int *temp, int size, int rank, int gran){
	int row[node];
	int col[node];
	int l, u,add_row, act_col, act_row, val;
	int k;
	for(l=0; l < gran * node; l++){
		add_row = l/node;
		act_col = l%node;
		act_row = gran * rank + act_row;
		get_row(data, row, act_row);
		//for (k = 0; k<node;k++)
		// printf("%d\t", row[k]);
		get_col(raw, col, act_col);
		val = 0;
		for(u = 0; u < node; u++){
			val += data[u] * raw[u];
		}
		*(temp + l) = val;

	}
	return;
}

void get_row(int *data, int *row, int num_row){
	int i;
	for(i=0; i < node; i++){
		*(row + i) = *(data + node * num_row + i);
	}
}
void get_col(int *data, int *col, int num_col){
	int i;
	for(i=0; i < node; i++){
		*(col + i) = *(data + num_col + node * i );
	}
}
