/*
CPSC 521 assignment 1: Nbody program based on [1].

[1] Wilkinson B., Allen M., "Parallel Programming: Techniques and Applications Using Networked
Workstations and Parallel Computers", Prentice-Hall, 2005 (2nd Edition): page 126-130


  (c) Shu Yang 2012
  Email: syang11@cs.ubc.ca

*/

#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

const float deltaT = 1.0;	//delta t constant
const double G= 0.0000000000667384;		//gravitational constant

//main function
//argv: (rounds, granularity, filename)
int main(int argc, char *argv[]) {
	//initilize MPI
	int test, size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//exit if the number of processes is not satisfied
	if(size<2)
	{
		printf("Need at least 2 processes!\n");
		MPI_Abort(MPI_COMM_WORLD,99);
	}

	MPI_Status status;

	int rounds = atoi(argv[1]);
	char* filename= argv[2];

	float data[size * 3];


	int i = 0;
	int j = 0;
	int k = 0;

	//find the neighbors
	int previous = rank - 1;
	int next = rank + 1;
	if (next == size) { //upperbound
		next = 0;
	}
	if (previous == -1){	//lowerbound
		previous = size - 1;
	}

	//INITIALIZATION
	if (rank == 0) {
		//read from input file
		FILE *fp;
		fp = fopen(filename, "r");

		for (i = 0; i < size * 3; i++) {
			fscanf(fp, "%f", data+i);
		}

		MPI_Send(data + 3 , (size - rank-1) * 3, MPI_FLOAT, next, 0,
				MPI_COMM_WORLD);
	} else if (rank == size - 1) {
		MPI_Recv(data, (size - rank) * 3, MPI_FLOAT, previous, 0, MPI_COMM_WORLD,
				&status);
	} else {
		MPI_Recv(data, (size - rank) * 3, MPI_FLOAT, previous, 0,
				MPI_COMM_WORLD, &status);
		MPI_Send(data + 3, (size - rank - 1) * 3, MPI_FLOAT, next,
				0, MPI_COMM_WORLD);
	}

	//SIMULATION
	float x, y;
	float vx, vy;
	for (i = 0; i < rounds; i++) {
		int tag_code=i; //set tag_code as the current round
		float *p = data;

		for (j = 0; j < size; j++) {
			MPI_Send(p, 3, MPI_FLOAT, next, tag_code, MPI_COMM_WORLD);
			p += 3;
			MPI_Recv(p, 3, MPI_FLOAT, previous, tag_code, MPI_COMM_WORLD,
					&status);
		}



		float mass = *(p + 2);

		//compute force
		float x_new = 0, y_new = 0, r2, fc;
		for (k = 0; k < size; k++) {
			//if r2 is 0, there is no force
			r2 = pow((data[k * 3] - data[0]), 2) + pow((data[k * 3 + 1] - data[1]), 2);
			if (r2 != 0){
				fc = G * data[k * 3 + 2] * mass / r2;
			}else{
				fc = 0;
			}

			//if the change for x (or y) is 0, there is no force
			if ((data[k * 3] - data[0]) != 0 || r2 != 0){
				x_new += fc * (data[k * 3] - data[0]) / sqrt(r2);
			}

			if ((data[k * 3 + 1] - data[1]) != 0 || r2 != 0){
				y_new += fc * (data[k * 3 + 1] - data[1]) / sqrt(r2);
			}
		}
		x = x_new;
		y = y_new;
		vx += x * deltaT / mass;
		vy += y * deltaT / mass;
		data[0] += vx * deltaT;
		data[1] += vy * deltaT;

	}

	//TERMINATION
	for (i = 0; i < size - 1; i++) {
		MPI_Send(data+i * 3, 3, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
		MPI_Recv(data+(i+1) * 3, 3, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
				&status);
	}


	//head node prints the results
	if (rank == 0) {
		for (i = 0; i < size; i++) {
			printf("%f\t%f\t%f\n", data[i*3], data[i*3+1], data[i*3+2]);
		}
	}

	MPI_FINALIZE();

}

