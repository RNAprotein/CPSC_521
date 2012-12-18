/*
CPSC 521 assignment 4: Visualize Nbody using MPE

Since the gravitational constant is too small that the bodies move too slow to demonstrate, use
gravitational constant * 1000000000 instead.


To build the program:
	mpicc hw4Nbody.c -o hw4Nbody -lm -lmpe -lX11

To run the program:
	mpiexec -n 10 ./hw4Nbody 50 data4Nbody.txt

  (c) Shu Yang 2012
  Email: syang11@cs.ubc.ca

*/

#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpe.h>
#include <mpe_graphics.h>

const float deltaT = 1.0;	//delta t constant
const double G= 0.0067;		//gravitational constant * 1000000000 for visualization purpose


const int WIDTH = 1000; //window width
const int HEIGHT = 1000; //window height
const int POSITION_DIMENSION = 3;	//position info (x, y, z)

MPE_XGraph window;

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

	//open a window
	MPE_Open_graphics(&window,MPI_COMM_WORLD,0,0,0,WIDTH,HEIGHT,0);

	MPI_Status status;

	int rounds = atoi(argv[1]);
	char* filename= argv[2];

	float data[size * POSITION_DIMENSION];


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

		for (i = 0; i < size * POSITION_DIMENSION; i++) {
			fscanf(fp, "%f", data+i);
		}

		MPI_Send(data + POSITION_DIMENSION , (size - rank-1) * POSITION_DIMENSION, MPI_FLOAT, next, 0,
				MPI_COMM_WORLD);
	} else if (rank == size - 1) {
		MPI_Recv(data, (size - rank) * POSITION_DIMENSION, MPI_FLOAT, previous, 0, MPI_COMM_WORLD,
				&status);
	} else {
		MPI_Recv(data, (size - rank) * POSITION_DIMENSION, MPI_FLOAT, previous, 0,
				MPI_COMM_WORLD, &status);
		MPI_Send(data + POSITION_DIMENSION, (size - rank - 1) * POSITION_DIMENSION, MPI_FLOAT, next,
				0, MPI_COMM_WORLD);
	}

	//SIMULATION
	float x, y;
	float vx, vy;
	for (i = 0; i < rounds; i++) {
		int tag_code=i; //set tag_code as the current round
		float *p = data;

		for (j = 0; j < size; j++) {
			MPI_Send(p, POSITION_DIMENSION, MPI_FLOAT, next, tag_code, MPI_COMM_WORLD);
			p += POSITION_DIMENSION;
			MPI_Recv(p, POSITION_DIMENSION, MPI_FLOAT, previous, tag_code, MPI_COMM_WORLD,
					&status);
		}

		float mass = *(p + 2);

		//compute force
		float x_new = 0, y_new = 0, r2, fc;
		for (k = 0; k < size; k++) {
			//if r2 is 0, there is no force
			r2 = pow((data[k * POSITION_DIMENSION] - data[0]), 2) + pow((data[k * POSITION_DIMENSION + 1] - data[1]), 2);
			if (r2 != 0){
				fc = G * data[k * POSITION_DIMENSION + 2] * mass / r2;
			}else{
				fc = 0;
			}

			//if the change for x (or y) is 0, there is no force
			if ((data[k * POSITION_DIMENSION] - data[0]) != 0 || r2 != 0){
				x_new += fc * (data[k * POSITION_DIMENSION] - data[0]) / sqrt(r2);
			}

			if ((data[k * POSITION_DIMENSION + 1] - data[1]) != 0 || r2 != 0){
				y_new += fc * (data[k * POSITION_DIMENSION + 1] - data[1]) / sqrt(r2);
			}
		}
		x = x_new;
		y = y_new;
		vx += x * deltaT / mass;
		vy += y * deltaT / mass;
		data[0] += vx * deltaT;
		data[1] += vy * deltaT;

		/*for (i= 0; i < size; i++)
		{
		  MPE_Draw_circle (window, data[3*i],data[3*i+1] ,4, MPE_BLUE);
		}*/
		char *string="sy";
		MPE_Draw_string(window, data[0],data[1] ,rank+1,string);
	    MPE_Update(window);
	    usleep(100000); //sleep 100ms time

	    if(i < rounds - 1) {
	    	MPE_Draw_string(window, data[0],data[1] ,MPE_WHITE,string);
	    }
	}

	//TERMINATION
	for (i = 0; i < size - 1; i++) {
		MPI_Send(data+i * POSITION_DIMENSION, POSITION_DIMENSION, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
		MPI_Recv(data+(i+1) * POSITION_DIMENSION, POSITION_DIMENSION, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
				&status);
	}


	//head node prints the results
	if (rank == 0) {
		for (i = 0; i < size; i++) {
			printf("%f, %f, %f\n", data[i*POSITION_DIMENSION], data[i*POSITION_DIMENSION+1], data[i*POSITION_DIMENSION+2]);
		}
	}

	//Since MPE does not provide any interface functions to interact with the users,
	//here I use touch-like design: exit once the user drag a square (or whatever)
	int startx, starty, endx, endy, button, dragVisual=2;		//2: line, 1: rect, 3: cycle
	if(rank == 0) {
	    MPE_Get_drag_region(window, button, dragVisual, &startx, &starty, &endx, &endy);
	}
	MPE_Close_graphics(&window);
	MPI_FINALIZE();

}
