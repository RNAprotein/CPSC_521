#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpe.h>
#include <mpe_graphics.h>

const double GRAVITY = 0.00000667384;
const double dt = 1000.0;
const int MAX_LINE_LENGTH = 40; //max line length

#define RADIUS 4
#define WINDOW_X 1000
#define WINDOW_Y 1000
MPE_XGraph window;

typedef struct {
  float x;
  float y;
  float mass;
} Particle;

Particle *particles;
int rank, size, granularity;

void simulation(int rounds);

void initialize();

void termination();

void recalculateForce(Particle *other, double *forceV);

void updatePositions(double *forceV);

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rounds;
  char *fileName;
  if (argc == 4) {
      rounds = atoi(argv[1]);
      granularity = atoi(argv[2]);
      fileName = argv[3];
    } else {
      if (rank ==0)
	printf("USAGE ./assignment1 ROUNDS GRANULARITY FILENAME\n");
      MPI_Finalize();
      return 1;
  }
  MPE_Open_graphics( &window, MPI_COMM_WORLD, 0, -1, -1, WINDOW_X,
		     WINDOW_Y, 0 );
  particles = (Particle *)malloc(sizeof(Particle) * granularity);
  double ib = MPI_Wtime();
  initialize(fileName);
  double ie = MPI_Wtime();
  double sb = MPI_Wtime();
  simulation(rounds);
  double se = MPI_Wtime();
  double gb = MPI_Wtime();
  termination();
  double ge = MPI_Wtime();
  if (rank == 0)
    printf("%.10f, %.10f, %.10f\n", ie - ib, se - sb, ge - gb);
  MPE_Close_graphics(&window);
  MPI_Finalize();
}

//Simulates the ring communication and recalculates new value of force at each step
//updates positions for each round
void simulation(int rounds) {
  int i, j;
  int partSize = sizeof(Particle) * granularity;
  Particle *buff = (Particle*)malloc(sizeof(Particle) * granularity);
  for (i = 0; i < rounds; i++) {
    double forceV[2 * granularity];  //force vectors
    if (!rank) {
      MPE_Fill_rectangle (window, 0, 0, WINDOW_X, WINDOW_Y, MPE_WHITE);
    }
    //if (!rank) {
      int index;
      for (index = 0; index < granularity; index++)
	{
	  MPE_Draw_circle (window, (float)(particles[index].x + WINDOW_X / 2), (float)(particles[index].y + WINDOW_Y / 2), RADIUS, MPE_BLUE);
	}
      //}
    for (j = 0; j < 2 * granularity; j++) //reset force vectors at every round
      forceV[j] = 0.0;
    recalculateForce(particles, forceV);  //calculate force locally
    if (size > 1) {
      int n = 0;
      if (rank == 0){
	MPI_Recv(buff, partSize, MPI_BYTE, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(particles, partSize, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
      } else {
	MPI_Send(particles, partSize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
	MPI_Recv(buff, partSize, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      while (n < size - 1) {
	recalculateForce(buff, forceV);
	MPI_Send(buff, partSize, MPI_BYTE, (rank + 1) % size, 0, MPI_COMM_WORLD);
	MPI_Recv(buff, partSize, MPI_BYTE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	n++;
      }
  }
    updatePositions(forceV);
    MPE_Update(window);
    usleep(100000);
  }
  free(buff);
}

//Recalculates new value of force
void recalculateForce(Particle *other, double* forceV) {
  int i, j;
  for (i = 0; i < granularity; i++) {  //calculate force locally
      for (j = 0; j < granularity; j++) {
	  double rsq = pow(particles[i].x - other[j].x, 2) + pow(particles[i].y - other[j].y, 2);
	  if (rsq > 0) {
	    double r = sqrt(rsq);
	    double force = GRAVITY * particles[i].mass * other[j].mass / rsq;
	    forceV[i * 2] += force * (-particles[i].x + other[j].x) / r;
	    forceV[i * 2 + 1] += force * (-particles[i].y + other[j].y) / r;
	  }
	}
      }
}

//Updates value of x and y positions according to total force from other particles
void updatePositions(double* forceV) {
  int i;
  for (i = 0; i < granularity; i++) {
    double vx = forceV[i * 2] * dt / particles[i].mass;
    double vy = forceV[i * 2 + 1] * dt / particles[i].mass;
    particles[i].x += vx * dt;
    particles[i].y += vy * dt;
  }
}

//Saves values in a struct and forwards the rest to the next processor
void saveParticle(char *buffer) {
  float mass, x, y;
  int j;
  for (j = 0; j < granularity; j++)
    {
      sscanf(buffer + MAX_LINE_LENGTH * j, "%f %f %f", &x, &y, &mass);
      particles[j].x = x;
      particles[j].y = y;
      particles[j].mass = mass;
    }
}

//Root processor read from the file into a string buffer and forwards to the next process
void initialize(char *fileName) {
  char *buffer = (char*)malloc(MAX_LINE_LENGTH * granularity * sizeof(char));
  FILE *file;
  if (rank == 0)
  {
    file = fopen(fileName, "rt");
    if (file == NULL) {
      printf("File doesn't exist");
      exit(1);
    }
  }
  int i = 0;
  int n = rank;
  while (n < size) {
    if (rank != 0) {
       MPI_Recv(buffer, MAX_LINE_LENGTH * granularity, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      int index;
      for (index = 0; index < granularity; index++) {//read a chunk
          fgets(buffer + index * MAX_LINE_LENGTH, MAX_LINE_LENGTH, file);
      }
    }
      if (i == 0) {
	saveParticle(buffer);
	i = 1;
      } else {
	MPI_Send(buffer, MAX_LINE_LENGTH * granularity, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
      }
      n++;
  }
  if (rank == 0)
    fclose(file);
  free(buffer);
}

//forwards particle information down and root proc prints them
void termination() {
  int n = rank;
  int sent = 0;
  while (n < size) { //reuse particle to forward
      if (sent == 1) {
	MPI_Recv(particles, sizeof(Particle) * granularity, MPI_BYTE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if (rank == 0) {
	int i;
	//	for (i = 0; i < granularity; i++)
	  //	  printf("%f %f %.2f\n", particles[i].x, particles[i].y, particles[i].mass);
      } else {
       	MPI_Send(particles, sizeof(Particle) * granularity, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD);
      }
      sent = 1;
      n++;
  }
  free(particles);
}
