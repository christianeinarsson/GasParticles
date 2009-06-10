#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "debug.h"

int msglevel = 60;

int main(int argc, char *argv[])
{
	int np;
	int me;
	{
		int ierr = MPI_Init(& argc, & argv);
		MPI_Comm_size(MPI_COMM_WORLD, & np);
		MPI_Comm_rank(MPI_COMM_WORLD, & me);
	}
	pmesg(99, "Hello from node %d\n", me);

	int isMaster = (me == 0 ? 1 : 0);
	int sliceAbove = (me - 1);
	int sliceBelow = (me + 1);
	int isFirstSlice = (me == 0 ? 1 : 0);
	int isLastSlice = (me == (np - 1) ? 1 : 0);

	// Check preconditions
	if(isMaster)
	{
		if(np < 3)
		{
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
			pcord_t pcord_t_item; // Element of the type
			int pcord_t_block_lengths[] = { 1, 1, 1, 1 }; // Set lengths of type elements
			MPI_Datatype pcord_t_block_types[] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT }; // Set types
			MPI_Aint pcord_t_start, pcord_t_displ[4];
			MPI_Address(& pcord_t_item, & pcord_t_start);
			MPI_Address(& pcord_t_item.x, & pcord_t_displ[0]);
			MPI_Address(& pcord_t_item.y, & pcord_t_displ[1]);
			MPI_Address(& pcord_t_item.vx, & pcord_t_displ[2]);
			MPI_Address(& pcord_t_item.vy, & pcord_t_displ[3]);
			pcord_t_displ[0] -= pcord_t_start; // Displacement relative to address of start
			pcord_t_displ[1] -= pcord_t_start; // Displacement relative to address of start
			pcord_t_displ[2] -= pcord_t_start; // Displacement relative to address of start
			pcord_t_displ[3] -= pcord_t_start; // Displacement relative to address of start
			MPI_Type_create_struct(4, pcord_t_block_lengths, pcord_t_displ, pcord_t_block_types, & MPI_COORDINATE);
			MPI_Type_commit(& MPI_COORDINATE);
		}

		// MPI commiting type slice_data_t
		{
			slice_data_t slice_data_t_item; // Element of the type
			int slice_data_t_block_lengths[] = { 1, 1, 1, 1, 1, 1 }; // Set lengths of type elements
			MPI_Datatype slice_data_t_block_types[] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_FLOAT }; // Set types
			MPI_Aint slice_data_t_start, slice_data_t_displ[6];
			MPI_Address(& slice_data_t_item, & slice_data_t_start);
			MPI_Address(& slice_data_t_item.current, & slice_data_t_displ[0]);
			MPI_Address(& slice_data_t_item.keeping, & slice_data_t_displ[1]);
			MPI_Address(& slice_data_t_item.sendingUp, & slice_data_t_displ[2]);
			MPI_Address(& slice_data_t_item.sendingDown, & slice_data_t_displ[3]);
			MPI_Address(& slice_data_t_item.wallCollisions, & slice_data_t_displ[4]);
			MPI_Address(& slice_data_t_item.pressure, & slice_data_t_displ[5]);
			slice_data_t_displ[0] -= slice_data_t_start; // Displacement relative to address of start
			slice_data_t_displ[1] -= slice_data_t_start; // Displacement relative to address of start
			slice_data_t_displ[2] -= slice_data_t_start; // Displacement relative to address of start
			slice_data_t_displ[3] -= slice_data_t_start; // Displacement relative to address of start
			slice_data_t_displ[4] -= slice_data_t_start; // Displacement relative to address of start
			slice_data_t_displ[5] -= slice_data_t_start; // Displacement relative to address of start
			MPI_Type_create_struct(6, slice_data_t_block_lengths, slice_data_t_displ, slice_data_t_block_types, & MPI_SLICE_DATA);
			MPI_Type_commit(& MPI_SLICE_DATA);
		}
	}

	if(isMaster) pmesg(97, "Divide rectangle into horizontal slices\n");
	cord_t slice;
	{
		slice.x0 = 0;
		slice.x1 = BOX_HORIZ_SIZE;
		slice.y0 = me * (BOX_VERT_SIZE / np);
		slice.y1 = (me + 1) * (BOX_VERT_SIZE / np);
	}

	if(isMaster) pmesg(97, "Initiate walls\n");
	cord_t walls;
	{
		walls.x0 = 0;
		walls.x1 = BOX_HORIZ_SIZE;
		walls.y0 = 0;
		walls.y1 = BOX_VERT_SIZE;
	}

	if(isMaster) pmesg(97, "Initiate particles\n");
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
		printParticles(90, sd->current, particles, me, "particles");

		pmesg(90, "Slice data at startup in node\n");
		pmesg(90, "ID\tCURR\tKEEP\tUP  \tDOWN\tPRESSURE\tWCOL\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	pmesg(90, "%3d\t%4d\t%4d\t%4d\t%4d\t%4.04f\t\t%4d\n", me, sliceData[me].current, sliceData[me].keeping, sliceData[me].sendingUp, sliceData[me].sendingDown, sliceData[me].pressure, sliceData[me].wallCollisions);
	MPI_Barrier(MPI_COMM_WORLD);

	if(isMaster) pmesg(97, "Main loop: for each time-step do\n");
	for(int t = 0; t < TIME_STEPS; t += TIME_STEP)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(isMaster) pmesg(75, "Time is %d\n", t);

		sd->keeping = 0;
		sd->sendingUp = 0;
		sd->sendingDown = 0;

		// Outer particle loop
		if(isMaster) pmesg(97, "for all particles do\n");
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

			if(isMaster) pmesg(97, "Check for wall interaction and add the momentum\n");
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

			if(escaped == 0)
			{
				particlesNext[sd->keeping] = particles[po];
				sd->keeping++;
			}
		}

		if(isMaster) pmesg(97, "Communicate if needed\n");

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
		if(isMaster){
			printParticles(85, sd->current, particles, me, "particles");
			printParticles(85, sd->keeping, particlesNext, me, "particlesNext");
			printParticles(85, sd->sendingUp, particlesSendUp, me, "particlesSendUp");
			printParticles(85, sd->sendingDown, particlesSendDown, me, "particlesSendDown");
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Send up, but only if there's something to send
		if(!isFirstSlice && sd->sendingUp != 0)
		{
			pmesg(70, "Sending %d->%d\n", me, sliceAbove);
			MPI_Send(particlesSendUp, sd->sendingUp, MPI_COORDINATE, sliceAbove, TAG_SEND_UP, MPI_COMM_WORLD);
		}

		// Receive from below, but only if there's something to receive
		if(!isLastSlice && sliceData[sliceBelow].sendingUp != 0)
		{
			pmesg(70, "Receiving %d<-%d (%d coordinates)\n", me, sliceBelow, sliceData[sliceBelow].sendingUp);
			pcord_t * particlesReceive = (pcord_t *) malloc(sizeof(pcord_t) * sliceData[sliceBelow].sendingUp);

			MPI_Status status;
			// TODO: Add MPI_Probe first to get buffer size of incoming data
			// http://www.mpi-forum.org/docs/mpi-11-html/node50.html#Node50
			MPI_Recv(particlesReceive, sliceData[sliceBelow].sendingUp, MPI_COORDINATE, sliceBelow, TAG_SEND_UP, MPI_COMM_WORLD, &status);
			//pmesg(70, "Received %d<-%d with status %d and size %d\n", me, sliceBelow, status->ERROR, status->size);
			//pmesg(70, "Received %d<-%d with status %d and size %s\n", me, sliceBelow, status);

			printParticles(85, sd->keeping, particlesNext, me, "particlesNext");

			// Insert sorted according to x position if possible?
			// Would increase the chances of a collision with a neighbor
			// and since collisions cut calculations, it would be efficient
			for(int pr = 0; pr < sliceData[sliceBelow].sendingUp; pr++)
			{
				pmesg(70, "Adding %d[%d]<-%d[%d]\n", me, sd->keeping, sliceBelow, pr);
				// TODO: THIS IS WHERE IT BREAKS *************************************************************************
				// Segfault if a particle has escaped to here, need to allocate memory by MPI_Probe above
				//particlesNext[sd->keeping] = particlesSendUp[pr];
				pmesg(70, "%d says %d[%d] = %d\t%d\t%d\t%d\n", me, sliceBelow, pr, particlesReceive[pr].x, particlesReceive[pr].y, particlesReceive[pr].vx, particlesReceive[pr].vy);
				particlesNext[sd->keeping].x = particlesReceive[pr].x;
				particlesNext[sd->keeping].y = particlesReceive[pr].y;
				particlesNext[sd->keeping].vx = particlesReceive[pr].vx;
				particlesNext[sd->keeping].vy = particlesReceive[pr].vy;
				//memcpy(&particlesNext[sd->keeping], &particlesSendUp[pr], sizeof(pcord_t));
				sd->keeping++;
				pmesg(70, "Added %d[%d]<-%d[%d]\n", me, sd->keeping, sliceBelow, pr);
			}
			free(particlesReceive);
		}
		if(isMaster)
		{
			pmesg(95, "Done sending/receiving up %d\n", me);
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
			MPI_Send(particlesSendDown, sd->sendingDown, MPI_COORDINATE, sliceBelow, TAG_SEND_DOWN, MPI_COMM_WORLD);
		}

		// Receive from above, but only if there's something to receive
		if(!isFirstSlice && sliceData[sliceAbove].sendingDown != 0)
		{
			pmesg(70, "Receiving %d<-%d (%d coordinates)\n", me, sliceAbove, sliceData[sliceAbove].sendingDown);
			//pcord_t * particlesReceive = (pcord_t *) malloc(sizeof(pcord_t) * sliceData[sliceAbove].sendingDown);

			// TODO: Add MPI_Probe first to get buffer size of incoming data
			// http://www.mpi-forum.org/docs/mpi-11-html/node50.html#Node50
			MPI_Status status;
			MPI_Recv(particlesSendDown, sliceData[sliceAbove].sendingDown, MPI_COORDINATE, sliceAbove, TAG_SEND_DOWN, MPI_COMM_WORLD, &status);

			printParticles(85, sd->keeping, particlesNext, me, "particlesNext");

			// Insert sorted according to x position if possible?
			// Would increase the chances of a collision with a neighbor
			// and since collisions cut calculations, it would be efficient
			for(int pr = 0; pr < sliceData[sliceAbove].sendingDown; pr++)
			{
				// TODO: THIS IS WHERE IT BREAKS *************************************************************************
				// Segfault if a particle has escaped to here, need to allocate memory by MPI_Probe above
				//particlesNext[sd->keeping] = particlesSendDown[pr];
				// Try to copy values one by one
				particlesNext[sd->keeping].x = particlesSendDown[pr].x;
				particlesNext[sd->keeping].y = particlesSendDown[pr].y;
				particlesNext[sd->keeping].vx = particlesSendDown[pr].vx;
				particlesNext[sd->keeping].vy = particlesSendDown[pr].vy;
				sd->keeping++;
			}
			//free(particlesReceive);
		}
		if(isMaster)
		{
			pmesg(95, "Done sending/receiving down %d\n", me);
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

	MPI_Barrier(MPI_COMM_WORLD);
	if(isMaster)
	{
		pmesg(40, "Slice data when done\n");
		pmesg(40, "ID\tCURR\tKEEP\tUP  \tDOWN\tPRESSURE\tWCOL\n");
		for(int i = 0; i < np; i++){
			pmesg(40, "%3d\t%4d\t%4d\t%4d\t%4d\t%4.04f\t\t%4d\n", i, sliceData[i].current, sliceData[i].keeping, sliceData[i].sendingUp, sliceData[i].sendingDown, sliceData[i].pressure, sliceData[i].wallCollisions);
		}

		pmesg(70, "Calculate pressure\n");
		{
			float totalPressure = 0;
			float averagePressure;
			for(int i = 0; i < np; i++){
				totalPressure += sliceData[i].pressure;
			}
			averagePressure = (totalPressure / (WALL_LENGTH * TIME_STEPS));
			pmesg(0, "Calculated pressure to be %.5f in a box with size %.1f by %.1f over %d time steps with %d particles.\n", averagePressure, BOX_HORIZ_SIZE, BOX_VERT_SIZE, TIME_STEPS, INIT_NO_PARTICLES*np);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
