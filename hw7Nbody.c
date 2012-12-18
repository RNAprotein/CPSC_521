/*
CPSC 521 assignment 7: enhanced Nbody program with MPI_I/O support.


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
	MPI_Offset filesize;

	char* filename= argv[2];
	int rounds = atoi(argv[1]);


	float *data;


	int i = 0;
	int j = 0;
	int k = 0;
	int bufsize, count;
	float *buf;

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

	//each process read file individually, no need for the head to pass the file for one round
	MPI_File infile;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);

	MPI_File_get_size(infile, &filesize); // in bytes
	filesize = filesize/sizeof(float);
	bufsize = filesize/size + 1; // local number to read

	buf=(float *)malloc(bufsize * sizeof(float));
	MPI_File_set_view(infile, rank * bufsize * sizeof(float), MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
	MPI_File_read(infile, buf, bufsize, MPI_FLOAT, &status);
	MPI_Get_count(&status, MPI_FLOAT, &count);

	printf("process %d read %d ints\n", rank, count);
	MPI_File_close(&infile);

	data=buf;

	//SIMULATION
	float x, y;
	float vx, vy;
	for (i = 0; i < rounds; i++) {
		int tag_code=i; //set tag_code as the current round
		float *p = data;

		for (j = 0; j < size-1; j++) {
			MPI_Send(p, 3, MPI_FLOAT, next, tag_code, MPI_COMM_WORLD);
			p += 3;
			MPI_Recv(p, 3, MPI_FLOAT, previous, tag_code, MPI_COMM_WORLD,&status);
		}

		float mass = data[2];

		//compute force
		float x_new = 0, y_new = 0, r2, fc;
		for (k = 1; k < size; k++) {
			//if r2 is 0, there is no force
			r2 = pow((data[k * 3] - data[0]), 2) + pow((data[k * 3 + 1] - data[1]), 2);
			if (r2 != 0){
				fc = G * data[k * 3 + 2] * mass / r2;
			}else{
				fc = 0;
			}

			//if the change for x (or y) is 0, there is no force
			if ((data[k * 3] - data[0]) != 0){
				x_new += fc * (data[k * 3] - data[0]) / sqrt(r2);
			}

			if ((data[k * 3 + 1] - data[1]) != 0){
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


	//each process prints its results
	MPI_File outfile;
	MPI_File_open(MPI_COMM_WORLD, "testfile", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);
	MPI_File_set_view(outfile, rank * bufsize * sizeof(float),MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
	MPI_File_write(outfile, data, bufsize, MPI_FLOAT, MPI_STATUS_IGNORE);
	MPI_File_close(&outfile);

	MPI_FINALIZE();

}
