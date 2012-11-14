/*
CPSC 521 assignment 2: Nbody with granularity.

The original version is simple N-body script without granularity.

The later version, according to assignment 2, has granularity bodies

  (c) Shu Yang 2012
  Email: syang11@cs.ubc.ca

*/

#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

const float deltaT = 1.0;	//delta t constant
const float G = 6.67384 * pow(10, -11); //G constant

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
	int granularity = atoi(argv[2]);
	char* filename= argv[3];

	float T[granularity * 3];
	int bodies = size * granularity;
	float data[bodies * 3];


	//find the neighbors
	int previous = rank - 1;
	int next = rank + 1;
	if (next == size) { //upperbound
		next = 0;
	}
	if (previous == -1){	//lowerbound
		previous = size - 1;
	}

	//initialization
	if (rank == 0) {
		//read from input file
		FILE *fp;
		fp = fopen(filename, "r");
		float *p = data;
		for (int i = 0; i < bodies * 3; i++) {
			fscanf(fp, "%f", p);
			p++;
		}
		p = data;

		MPI_Send(p + 3 * granularity, (bodies - granularity) * 3, MPI_FLOAT, next, 0,
				MPI_COMM_WORLD);
	} else if (rank == size - 1) {
		MPI_Recv(data, granularity * 3, MPI_FLOAT, previous, 0, MPI_COMM_WORLD,
				&status);
	} else {
		MPI_Recv(data, (size - rank) * granularity * 3, MPI_FLOAT, previous, 0,
				MPI_COMM_WORLD, &status);
		MPI_Send(data + 3 * granularity, (size - rank - 1) * granularity * 3, MPI_FLOAT, next,
				0, MPI_COMM_WORLD);
	}

	//simulation
	for (int i = 0; i < granularity * 3; i++) { //make a copy
		T[i] = data[i];
	}

	float x, y;
	for (int i = 0; i < rounds; i++) {
		int tag_code=i; //set tag_code as the current round
		float *p = data;

		for (int j = 0; j < size - 1; j++) {
			MPI_Send(data, granularity * 3, MPI_FLOAT, next, tag_code, MPI_COMM_WORLD);
			p += granularity * 3;
			MPI_Recv(p, granularity * 3, MPI_FLOAT, previous, tag_code, MPI_COMM_WORLD,
					&status);
		}

		float vx, vy;
		for (int j = 0; j < granularity; j++) {
			float mass = *(p + j * 3 + 2);

			//compute force
			float x_new = 0, y_new = 0, r2, fc;
			for (int k = 0; k < size * granularity; k++) {
				if (k != j) {

					//if the change for x (or y) is 0, there is no force
					if ((data[k * 3] - data[j * 3]) != 0 || r2 != 0){
						x_new += fc * (data[k * 3] - data[j * 3]) / sqrt(r2);
					}

					if ((data[k * 3 + 1] - data[j * 3 + 1]) != 0 || r2 != 0){
						y_new += fc * (data[k * 3 + 1] - data[j * 3 + 1]) / sqrt(r2);
					}

					//if r2 is 0, there is no force
					r2 = pow((data[k * 3] - data[j * 3]), 2) + pow((data[k * 3 + 1] - data[j * 3 + 1]), 2);
					if (r2 != 0){
						fc = G * data[k * 3 + 2] * mass / r2;
					}else{
						fc = 0;
					}
				}
			}
			x = x_new;
			y = y_new;
			vx += x * deltaT / mass;
			vy += y * deltaT / mass;
			T[j * 3] = data[j * 3] + vx * deltaT;
			T[j * 3+ 1 ] = data[j * 3 + 1] + vy * deltaT;
		}
		for (int j = 0; j < granularity; j++) {
			data[j * 3] = T[j * 3];
			data[j * 3 + 1] = T[j * 3 + 1];
		}
	}

	//termination
	if (rank == 0) {
		float *p = data;
		for (int i = 0; i < size - 1; i++) {
			MPI_Send(data, granularity * 3, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
			p += granularity * 3;
			MPI_Recv(p, granularity * 3, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
					&status);
		}
	} else {
		float *p = data;
		for (int i = 0; i < size - 1; i++) {
			MPI_Send(data, granularity * 3, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
			p += granularity * 3;
			MPI_Recv(p, granularity * 3, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
					&status);
		}
	}

	//head node prints the results
	if (rank == 0) {
		for (int i = 0; i < bodies; i++) {
			printf("%f\t%f\t%f\n", data[i*3], data[i*3+1], data[i*3+2]);
		}
	}

	MPI_FINALIZE();

}

