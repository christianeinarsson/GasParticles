#include <stdlib.h>
#include <math.h>
#include "definitions.h"
#include "debug.h"

int msglevel = 0;

int main(int argc, char *argv[])
{
	pemsg(30, "Initiate MPI\n");
	int np;
	int me;
	{
		int ierr = MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &np);
		MPI_Comm_rank(MPI_COMM_WORLD, &me);

		// MPI commiting type particle
		{
			particle_t item; // Element of the type
			MPI_Datatype MPI_PARTICLE; // MPI type to commit
			int block_lengths[] = { sizeof(pcord_t), sizeof(int) }; // Set lengths of type elements
			MPI_Datatype block_types[] = { MPI_UNSIGNED_CHAR (FIX), MPI_INT }; // Set types
			MPI_Aint start, displ[2];
			MPI_Address(&item, &start);
			MPI_Address(&item.r, &displ[0]);
			MPI_Address(&item.g, &displ[1]);
			displ[0] -= start; // Displacement relative to address of start
			displ[1] -= start; // Displacement relative to address of start
			MPI_Type_struct(2, block_lengths, displ, block_types, &MPI_PARTICLE);
			MPI_Type_commit(&MPI_PARTICLE);
		}
	}
	pemsg(35, "Initiate MPI: %d/%d\n", me, np);

	int master = (me == 0);

	if(master) pemsg("Divide rectangle into horizontal slices\n");
	float xstart = 0;
	float xstop = BOX_HORIZ_SIZE;
	float ystart = me * (BOX_VERT_SIZE/np);
	float ystop = (me + 1) * (BOX_VERT_SIZE/np);

	// Initiate particles
	particle_t * particles;
	{
		particles = (particle_t *) malloc(sizeof(particle_t) * INIT_NO_PARTICLES * np);

		for(int i = (me * INIT_NO_PARTICLES); i < ((me + 1) * INIT_NO_PARTICLES); i++){
			particles[i].pcord->x = xstart + (rand() * (xstop - xstart));
			particles[i].pcord->y = ystart + (rand() * (ystop - ystart));
			float r = rand() * MAX_INITIAL_VELOCITY;
			float theta = rand() * 2 * PI;
			particles[i].pcord->vx = r * cos(theta);
			particles[i].pcord->vy = r * sin(theta);
		}
	}

	if(master) pemsg(30, "Main loop: for each time-step do\n");
	// for
	{
		if(master) pemsg(25, "for all particles do\n");
		// for
		{
			if(master) pemsg(20, "Check for collisions\n");
			if(master) pemsg(20, "Move particles that has not collided with another\n");
			if(master) pemsg(20, "Check for wall interaction and add the momentum\n");
			if(master) pemsg(20, "Communicate if needed\n");
		}
	}

	if(master){
		pemsg(30, "Calculate pressure\n");
		{

		}
	}
}
