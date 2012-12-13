/*
CPSC 521 assignment 2: Nbody with granularity.

The original version is simple N-body script without granularity.

The later version, according to assignment 2, has granularity bodies

Compile on the head node:
[syang@cyclops Nbody]$ mpicc hw2Nbody.c -o hw2Nbody -lm

Run on the assigned node:
[syang@cyclops ~]$ ique1
qsub: waiting for job 8797.cyclops to start
qsub: job 8797.cyclops ready

----------------------------------------
Begin PBS Prologue Thu Dec 13 01:27:17 PST 2012 1355390837
Job ID:         8797.cyclops
Username:       syang
Group:          students
Rsrcs:          neednodes=1:ppn=8,nodes=1:ppn=8,pmem=1500000kb,pvmem=1500000kb,walltime=03:00:00
Queue:          dque
Nodes:          node26

Note: Global temp storage set at /scratch/tmp/8797.cyclops ($PBS_TMPDIR)
      Fast local temp storage set at /tmp/8797.cyclops     ($TMPDIR)

End PBS Prologue Thu Dec 13 01:27:17 PST 2012 1355390837
----------------------------------------

[syang@node26 ~]$ cat $PBS_NODEFILE
node26
node26
node26
node26
node26
node26
node26
node26

[syang@node26 Nbody]$ mpiexec -machinefile $PBS_NODEFILE -n 8 ./hw2Nbody 10 1000 data4Nbody.txt >> 8cores1K.txt





[syang@cyclops ~]$ ique2
qsub: waiting for job 8799.cyclops to start
qsub: job 8799.cyclops ready

----------------------------------------
Begin PBS Prologue Thu Dec 13 01:58:33 PST 2012 1355392713
Job ID:         8799.cyclops
Username:       syang
Group:          students
Rsrcs:          neednodes=2:ppn=8,nodes=2:ppn=8,pmem=1500000kb,pvmem=1500000kb,walltime=03:00:00
Queue:          dque
Nodes:          node21 node26

Note: Global temp storage set at /scratch/tmp/8799.cyclops ($PBS_TMPDIR)
      Fast local temp storage set at /tmp/8799.cyclops     ($TMPDIR)

End PBS Prologue Thu Dec 13 01:58:33 PST 2012 1355392713
----------------------------------------

[syang@node26 ~]$ cat $PBS_NODEFILE >> mymachinefile
[syang@node26 ~]$ vi mymachinefile 
[syang@node26 ~]$ more mymachinefile 
node26:8
node26:8
node26:8
node26:8
node26:8
node26:8
node26:8
node26:8
node21:8
node21:8
node21:8
node21:8
node21:8
node21:8
node21:8
node21:8
[syang@node26 ~]$ mpiexec -machinefile mymachinefile -n 16 521/a2/Nbody/hw2Nbody 10 1000 521/a2/Nbody/data4Nbody.txt >>521/a2/Nbody/16cores1K.txt





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
const int POSITION_DIMENSION = 3;	//position info (x, y, z)

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
	char* filename= argv[POSITION_DIMENSION];

	float T[granularity * POSITION_DIMENSION];
	int bodies = size * granularity;
	float data[bodies * POSITION_DIMENSION];


	int i = 0;
	int j = 0;
	int k = 0;
	double start, end, elapsed;

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

		for (i = 0; i < bodies * POSITION_DIMENSION; i++) {
			fscanf(fp, "%f", data+i);
		}
		
		start=MPI_Wtime();
		printf("....Start to distribute data at= %lf....\n", start);
		MPI_Send(data + POSITION_DIMENSION * granularity, (bodies - granularity) * POSITION_DIMENSION, MPI_FLOAT, next, 0,
				MPI_COMM_WORLD);
	} else if (rank == size - 1) {
		MPI_Recv(data, granularity * POSITION_DIMENSION, MPI_FLOAT, previous, 0, MPI_COMM_WORLD,
				&status);
		end==MPI_Wtime();
		elapsed=end-start;
		printf("''''Finish distribute data at= %lf\tElapsed time= %lf''''\n", end, elapsed);
	} else {
		MPI_Recv(data, (size - rank) * granularity * POSITION_DIMENSION, MPI_FLOAT, previous, 0,
				MPI_COMM_WORLD, &status);
		MPI_Send(data + POSITION_DIMENSION * granularity, (size - rank - 1) * granularity * POSITION_DIMENSION, MPI_FLOAT, next,
				0, MPI_COMM_WORLD);
	}

	//SIMULATION
	for (i = 0; i < granularity * POSITION_DIMENSION; i++) { //make a copy
		T[i] = data[i];
	}

	float x, y;
	float vx, vy;
	for (i = 0; i < rounds; i++) {
		int tag_code=i; //set tag_code as the current round
		float *p = data;

		for (j = 0; j < size - 1; j++) {
			MPI_Send(p, granularity * POSITION_DIMENSION, MPI_FLOAT, next, tag_code, MPI_COMM_WORLD);
			p += granularity * POSITION_DIMENSION;
			MPI_Recv(p, granularity * POSITION_DIMENSION, MPI_FLOAT, previous, tag_code, MPI_COMM_WORLD,
					&status);
		}


		for (j = 0; j < granularity; j++) {
			float mass = *(p + j * POSITION_DIMENSION + 2);

			//compute force
			float x_new = 0, y_new = 0, r2, fc;
			for (k = 0; k < size * granularity; k++) {
				if (k != j) {

					//if r2 is 0, there is no force
					r2 = pow((data[k * POSITION_DIMENSION] - data[j * POSITION_DIMENSION]), 2) + pow((data[k * POSITION_DIMENSION + 1] - data[j * POSITION_DIMENSION + 1]), 2);
					if (r2 != 0){
						fc = (G) * data[k * POSITION_DIMENSION + 2] * mass / r2;
					}else{
						fc = 0;
					}

					//if the change for x (or y) is 0, there is no force
					if ((data[k * POSITION_DIMENSION] - data[j * POSITION_DIMENSION]) != 0 || r2 != 0){
						x_new += fc * (data[k * POSITION_DIMENSION] - data[j * POSITION_DIMENSION]) / sqrt(r2);
					}

					if ((data[k * POSITION_DIMENSION + 1] - data[j * POSITION_DIMENSION + 1]) != 0 || r2 != 0){
						y_new += fc * (data[k * POSITION_DIMENSION + 1] - data[j * POSITION_DIMENSION + 1]) / sqrt(r2);
					}

				}
			}
			x = x_new;
			y = y_new;
			vx += x * deltaT / mass;
			vy += y * deltaT / mass;
			T[j * POSITION_DIMENSION] = data[j * POSITION_DIMENSION] + vx * deltaT;
			T[j * POSITION_DIMENSION+ 1 ] = data[j * POSITION_DIMENSION + 1] + vy * deltaT;
		}
		for (j = 0; j < granularity; j++) {
			data[j * POSITION_DIMENSION] = T[j * POSITION_DIMENSION];
			data[j * POSITION_DIMENSION + 1] = T[j * POSITION_DIMENSION + 1];
		}
	}

	//TERMINATION
	if (rank == 0) {
		start=MPI_Wtime();
		printf("....Start to gather data at: %lf....\n", start);
		for (i = 0; i < size - 1; i++) {
			MPI_Send(data+i*granularity * POSITION_DIMENSION, granularity * POSITION_DIMENSION, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
			MPI_Recv(data+(i+1)*granularity * POSITION_DIMENSION, granularity * POSITION_DIMENSION, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
					&status);
		}
		end==MPI_Wtime();
		elapsed=end-start;
		printf("''''Finish distribute data at= %lf\tElapsed time= %lf''''\n", end, elapsed);
	} else {
		for (i = 0; i < size - 1; i++) {
			MPI_Send(data+i*granularity * POSITION_DIMENSION, granularity * POSITION_DIMENSION, MPI_FLOAT, previous, 1, MPI_COMM_WORLD);
			MPI_Recv(data+(i+1)*granularity * POSITION_DIMENSION, granularity * POSITION_DIMENSION, MPI_FLOAT, next, 1, MPI_COMM_WORLD,
					&status);
		}
	}

	//head node prints the results
	if (rank == 0) {
		for (i = 0; i < bodies; i++) {
			printf("%f\t%f\t%f\n", data[i*POSITION_DIMENSION], data[i*POSITION_DIMENSION+1], data[i*POSITION_DIMENSION+2]);
		}
	}

	MPI_FINALIZE();

}

