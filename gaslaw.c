#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "debug.h"

int msglevel = 94;

int main(int argc, char *argv[])
{
	int np;
	int me;
	{
		int ierr = MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &np);
		MPI_Comm_rank(MPI_COMM_WORLD, &me);
	}
	pmesg(75, "Initiate MPI: %d/%d\n", me, np);

	int isMaster = (me == 0);
	int sliceAbove = (me-1);
	int sliceBelow = (me+1);
	int isFirstSlice = (me == 0 ? 1 : 0);
	int isLastSlice = (me == (np-1) ? 1 : 0);

	// Check preconditions
	if(isMaster){
		if(np < 3){
			fprintf(stderr, "Just %d MPI nodes is too few - need at least three.\n", np);
			MPI_Finalize();
			exit(1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Datatype MPI_COORDINATE; // MPI type to commit
	MPI_Datatype MPI_SLICE_DATA; // MPI type to commit
	{
		// MPI commiting type particle
		{
			pcord_t item; // Element of the type
			int block_lengths[] = { 1, 1, 1, 1 }; // Set lengths of type elements
			MPI_Datatype block_types[] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT }; // Set types
			MPI_Aint start, displ[4];
			MPI_Address(& item, & start);
			MPI_Address(& item.x, & displ[0]);
			MPI_Address(& item.y, & displ[1]);
			MPI_Address(& item.vx, & displ[2]);
			MPI_Address(& item.vy, & displ[3]);
			displ[0] -= start; // Displacement relative to address of start
			displ[1] -= start; // Displacement relative to address of start
			displ[2] -= start; // Displacement relative to address of start
			displ[3] -= start; // Displacement relative to address of start
			MPI_Type_struct(4, block_lengths, displ, block_types, & MPI_COORDINATE);
			MPI_Type_commit(& MPI_COORDINATE);
		}

		// MPI commiting type slice_send_size_t
		{
			slice_data_t item; // Element of the type
			int block_lengths[] = { 1, 1, 1, 1, 1, 1 }; // Set lengths of type elements
			MPI_Datatype block_types[] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_FLOAT }; // Set types
			MPI_Aint start, displ[6];
			MPI_Address(& item, & start);
			MPI_Address(& item.current, & displ[0]);
			MPI_Address(& item.keeping, & displ[1]);
			MPI_Address(& item.sendingUp, & displ[2]);
			MPI_Address(& item.sendingDown, & displ[3]);
			MPI_Address(& item.wallCollisions, & displ[4]);
			MPI_Address(& item.pressure, & displ[5]);
			displ[0] -= start; // Displacement relative to address of start
			displ[1] -= start; // Displacement relative to address of start
			displ[2] -= start; // Displacement relative to address of start
			displ[3] -= start; // Displacement relative to address of start
			displ[4] -= start; // Displacement relative to address of start
			displ[5] -= start; // Displacement relative to address of start
			MPI_Type_struct(6, block_lengths, displ, block_types, & MPI_SLICE_DATA);
			MPI_Type_commit(& MPI_SLICE_DATA);
		}
	}

	if(isMaster) pmesg(70, "Divide rectangle into horizontal slices\n");
	cord_t slice;
	{
		slice.x0 = 0;
		slice.x1 = BOX_HORIZ_SIZE;
		slice.y0 = me * (BOX_VERT_SIZE / np);
		slice.y1 = (me + 1) * (BOX_VERT_SIZE / np);
	}

	if(isMaster) pmesg(70, "Initiate walls\n");
	cord_t walls;
	{
		walls.x0 = 0;
		walls.x1 = BOX_HORIZ_SIZE;
		walls.y0 = 0;
		walls.y1 = BOX_VERT_SIZE;
	}

	if(isMaster) pmesg(70, "Initiate particles\n");
	pcord_t * particles;
	pcord_t * particlesNext;
	pcord_t * particlesSendUp;
	pcord_t * particlesSendDown;
	slice_data_t * sliceData;
	slice_data_t * sd;
	{
		particles = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesNext = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesSendUp = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesSendDown = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);

		sliceData = (slice_data_t *) malloc(sizeof(slice_data_t) * np);
		sd = & sliceData[me];

		sd->current = INIT_NO_PARTICLES;
		sd->keeping = 0;
		sd->sendingUp = 0;
		sd->sendingDown = 0;
		sd->wallCollisions = 0;
		sd->pressure = 0.0;

		// The functions start at the same time
		// Multiplier as spacing to make sure time doesn't
		// tick and gives me and me+1 the same value
		srand(time(NULL) + (me * 1024));
		for(int p = 0; p < sd->current; p++)
		{
			particles[p].x = slice.x0 + (rand01() * (slice.x1 - slice.x0));
			particles[p].y = slice.y0 + (rand01() * (slice.y1 - slice.y0));
			float r = rand01() * MAX_INITIAL_VELOCITY;
			float theta = rand01() * 2 * PI;
			particles[p].vx = r * cos(theta);
			particles[p].vy = r * sin(theta);
		}
	}

	// Print stats
	MPI_Barrier(MPI_COMM_WORLD);
	if(isMaster)
	{
		pmesg(70, "Slice data at startup in node\n");
		pmesg(70, "ID\tCURR\tKEEP\tUP  \tDOWN\tPRESSURE\tWCOL\n");
	}
	pmesg(70, "%3d\t%4d\t%4d\t%4d\t%4d\t%4.04f\t\t%4d\n", me, sliceData[me].current, sliceData[me].keeping, sliceData[me].sendingUp, sliceData[me].sendingDown, sliceData[me].pressure, sliceData[me].wallCollisions);
	MPI_Barrier(MPI_COMM_WORLD);

	if(isMaster) pmesg(70, "Main loop: for each time-step do\n");
	for(int t = 0; t < TIME_STEPS; t += TIME_STEP)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(isMaster) pmesg(75, "Time is %d\n", t);

		sd->keeping = 0;
		sd->sendingUp = 0;
		sd->sendingDown = 0;

		// Outer particle loop
		if(isMaster) pmesg(75, "for all particles do\n");
		for(int po = 0; po < sd->current; po++)
		{
			float pressure;

			// Inner particle loop
			// Checks only the lower part of the
			// (sd->current)^2 matrix (from po and down)
			for(int pi = po+1; pi < sd->current; pi++)
			{
				if(isMaster) pmesg(96, "Check for collisions between %d and %d\n", po, pi);
				float collision = collide(&particles[po], &particles[pi]);
				if(collision != -1)
				{
					if(isMaster) pmesg(80, "Collision between %d and %d, returned %f\n", po, pi, collision);
					interact(&particles[po], &particles[pi], collision);
					// Don't do further collision checks to cut calculations
					break;
				}else{
					if(isMaster) pmesg(95, "No collision between %d and %d\n", po, pi);
					feuler(&particles[po], TIME_STEP);
				}
			}

			if(isMaster) pmesg(80, "Check for wall interaction and add the momentum\n");
			// Check left and right wall collisions
			// Done first, they may bounce out of slice bounds
			pressure = wall_collide_leftright(&particles[po], walls);
			if(pressure !=  0.0)
			{
				sd->pressure += pressure;
				sd->wallCollisions++;
			}

			// Check top wall collisions only for the first slice
			// otherwise check for escaping particles
			int escaped = 0;
			if(isFirstSlice)
			{
				pressure = wall_collide_top(&particles[po], walls);
				if(pressure !=  0.0)
				{
					sd->pressure += pressure;
					sd->wallCollisions++;
				}
			}
			else if(particle_escape_top(&particles[po], slice))
			{
				particlesSendUp[sd->sendingUp] = particles[po];
				sd->sendingUp++;
				escaped = 1;
			}

			// Check bottom wall collisions only for the last slice
			// otherwise check for escaping particles
			if(isLastSlice)
			{
				pressure = wall_collide_bottom(&particles[po], walls);
				if(pressure !=  0.0)
				{
					sd->pressure += pressure;
					sd->wallCollisions++;
				}
			}
			else if(particle_escape_bottom(&particles[po], slice))
			{
				particlesSendDown[sd->sendingDown] = particles[po];
				sd->sendingDown++;
				escaped = 1;
			}

			if(!escaped)
			{
				particlesNext[sd->keeping] = particles[po];
				sd->keeping++;
			}
		}

		if(isMaster) pmesg(80, "Communicate if needed\n");

		// All processors knows all the other data
		// Might be overkill, could just send to sliceAbove, sliceBelow and perhaps root
		MPI_Allgather(sd, 1, MPI_SLICE_DATA, sliceData, 1, MPI_SLICE_DATA, MPI_COMM_WORLD);

		// Print stats
		MPI_Barrier(MPI_COMM_WORLD);
		if(isMaster)
		{
			pmesg(70, "Slice data for iteration %d\n", t);
			pmesg(70, "ID\tCURR\tKEEP\tUP  \tDOWN\tPRESSURE\tWCOL\n");
			for(int i = 0; i < np; i++){
				pmesg(70, "%3d\t%4d\t%4d\t%4d\t%4d\t%4.04f\t\t%4d\n", i, sliceData[i].current, sliceData[i].keeping, sliceData[i].sendingUp, sliceData[i].sendingDown, sliceData[i].pressure, sliceData[i].wallCollisions);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Send up, but only if there's something to send
		if(!isFirstSlice && sd->sendingUp != 0)
		{
			pmesg(70, "Sending %d->%d\n", me, sliceAbove);
			MPI_Send(&particlesSendUp, sd->sendingUp, MPI_COORDINATE, sliceAbove, TAG_SEND_UP, MPI_COMM_WORLD);
		}

		// Receive from below, but only if there's something to receive
		if(!isLastSlice && sliceData[sliceBelow].sendingUp != 0)
		{
			pmesg(70, "Receiving %d<-%d\n", me, sliceBelow);
			MPI_Status status;
			// TODO: Add MPI_Probe first to get buffer size of incoming data
			// http://www.mpi-forum.org/docs/mpi-11-html/node50.html#Node50
			MPI_Recv(&particlesSendUp, MAX_NO_PARTICLES, MPI_COORDINATE, sliceBelow, TAG_SEND_UP, MPI_COMM_WORLD, &status);
			pmesg(70, "Received %d<-%d\n", me, sliceBelow);

			// Insert sorted according to x position if possible?
			// Would increase the chances of a collision with a neighbor
			// and since collisions cut calculations, it would be efficient
			for(int pr = 0; pr < sliceData[sliceBelow].sendingUp; pr++)
			{
				pmesg(70, "Adding %d[%d]<-%d[%d]\n", me, sd->keeping, sliceBelow, pr);
				// TODO: THIS IS WHERE IT BREAKS *************************************************************************
				// Segfault if a particle has escaped to here, need to allocate memory by MPI_Probe above
				particlesNext[sd->keeping] = particlesSendUp[pr];
				//memcpy(&particlesNext[sd->keeping], &particlesSendUp[pr], sizeof(pcord_t));
				sd->keeping++;
				pmesg(70, "Added %d[%d]<-%d[%d]\n", me, sd->keeping, sliceBelow, pr);
			}
		}
		if(isMaster)
		{
			pmesg(70, "Done sending/receiving up %d\n", me);
		}
		else{
			pmesg(95, "Done sending/receiving up %d\n", me);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(isMaster) pmesg(70, "Done sending/receiving up iteration %d\n", t);

		// Send down, but only if there's something to send
		if(!isLastSlice && sd->sendingDown != 0)
		{
			pmesg(70, "Sending %d->%d\n", me, sliceBelow);
			MPI_Send(&particlesSendDown, sd->sendingDown, MPI_COORDINATE, sliceBelow, TAG_SEND_DOWN, MPI_COMM_WORLD);
		}

		// Receive from above, but only if there's something to receive
		if(!isFirstSlice && sliceData[sliceAbove].sendingDown != 0)
		{
			pmesg(70, "Receiving %d<-%d\n", me, sliceAbove);
			// TODO: Add MPI_Probe first to get buffer size of incoming data
			// http://www.mpi-forum.org/docs/mpi-11-html/node50.html#Node50
			MPI_Status status;
			MPI_Recv(&particlesSendDown, MAX_NO_PARTICLES, MPI_COORDINATE, sliceAbove, TAG_SEND_DOWN, MPI_COMM_WORLD, &status);

			// Insert sorted according to x position if possible?
			// Would increase the chances of a collision with a neighbor
			// and since collisions cut calculations, it would be efficient
			for(int pr = 0; pr < sliceData[sliceAbove].sendingDown; pr++)
			{
				// TODO: THIS IS WHERE IT BREAKS *************************************************************************
				// Segfault if a particle has escaped to here, need to allocate memory by MPI_Probe above
				particlesNext[sd->keeping] = particlesSendDown[pr];
				sd->keeping++;
			}
		}
		if(isMaster)
		{
			pmesg(70, "Done sending/receiving down %d\n", me);
		}
		else
		{
			pmesg(95, "Done sending/receiving down %d\n", me);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(isMaster) pmesg(70, "Done sending/receiving down iteration %d\n", t);

		// Use the kept particles for the next iteration
		sd->current = sd->keeping;
		swap(pcord_t *, particles, particlesNext);
	}

	// One last gather to calculate total pressure
	MPI_Gather(sd, 1, MPI_SLICE_DATA, sliceData, 1, MPI_SLICE_DATA, 0, MPI_COMM_WORLD);

	if(isMaster)
	{
		pmesg(70, "Slice data when done\n");
		pmesg(70, "ID\tCURR\tKEEP\tUP  \tDOWN\tPRESSURE\tWCOL\n");
		for(int i = 0; i < np; i++){
			pmesg(70, "%3d\t%4d\t%4d\t%4d\t%4d\t%4.04f\t\t%4d\n", i, sliceData[i].current, sliceData[i].keeping, sliceData[i].sendingUp, sliceData[i].sendingDown, sliceData[i].pressure, sliceData[i].wallCollisions);
		}

		pmesg(70, "Calculate pressure\n");
		{
			float totalPressure = 0;
			float averagePressure;
			for(int i = 0; i < np; i++){
				totalPressure += sliceData[i].pressure;
			}
			averagePressure = (totalPressure / (WALL_LENGTH * TIME_STEPS));
			pmesg(70, "Calculated pressure to be %f\n", averagePressure);
		}
	}
	MPI_Finalize();
}
