#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <VT.h>
#include <limits.h>
#include <float.h>
#include "definitions.h"
#include "debug.h"

int msglevel = 60;

int main(int argc, char *argv[])
{
	VT_initialize(&argc, &argv);
	int symdef_main;
	VT_funcdef("Main", VT_NOCLASS, &symdef_main);
	VT_enter(symdef_main, VT_NOSCL);

	int symdef_init;
	VT_funcdef("Init", VT_NOCLASS, &symdef_init);
	VT_enter(symdef_init, VT_NOSCL);

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
			VT_finalize();
			exit(1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// MPI types to commit
	MPI_Datatype MPI_COORDINATE;
	MPI_Datatype MPI_SLICE_DATA;
	{
		// MPI commiting type particle
		{
			// Element of the type
			pcord_t item;
			// Set lengths of type elements
			int block_lengths[] = { 1, 1, 1, 1 };
			// Set types
			MPI_Datatype block_types[] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT };
			MPI_Aint start, displ[4];
			MPI_Address(& item, & start);
			MPI_Address(& item.x, & displ[0]);
			MPI_Address(& item.y, & displ[1]);
			MPI_Address(& item.vx, & displ[2]);
			MPI_Address(& item.vy, & displ[3]);
			// Displacement relative to address of start
			displ[0] -= start;
			displ[1] -= start;
			displ[2] -= start;
			displ[3] -= start;
			MPI_Type_create_struct(4, block_lengths, displ, block_types, & MPI_COORDINATE);
			MPI_Type_commit(& MPI_COORDINATE);
		}

		// MPI commiting type slice_data_t
		{
			// Element of the type
			slice_data_t item;
			// Set lengths of type elements
			int block_lengths[] = { 1, 1, 1, 1, 1, 1, 1 };
			// Set types
			MPI_Datatype block_types[] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_FLOAT, MPI_DOUBLE };
			MPI_Aint start, displ[7];
			MPI_Address(& item, & start);
			MPI_Address(& item.current, & displ[0]);
			MPI_Address(& item.keeping, & displ[1]);
			MPI_Address(& item.sendingUp, & displ[2]);
			MPI_Address(& item.sendingDown, & displ[3]);
			MPI_Address(& item.wallCollisions, & displ[4]);
			MPI_Address(& item.pressure, & displ[5]);
			MPI_Address(& item.timing, & displ[6]);
			// Displacement relative to address of start
			displ[0] -= start;
			displ[1] -= start;
			displ[2] -= start;
			displ[3] -= start;
			displ[4] -= start;
			displ[5] -= start;
			displ[6] -= start;
			MPI_Type_create_struct(7, block_lengths, displ, block_types, & MPI_SLICE_DATA);
			MPI_Type_commit(& MPI_SLICE_DATA);
		}
	}

	if(isMaster) pmesg(97, "Initiate walls\n");
	cord_t walls;
	{
		walls.x0 = 0.0;
		walls.x1 = BOX_HORIZ_SIZE;
		walls.y0 = 0.0;
		walls.y1 = BOX_VERT_SIZE;
	}

	if(isMaster) pmesg(97, "Divide rectangle into horizontal slices\n");
	cord_t slice;
	{
		slice.x0 = 0.0;
		slice.x1 = BOX_HORIZ_SIZE;
		slice.y0 = me * (BOX_VERT_SIZE / np);
		slice.y1 = (me + 1) * (BOX_VERT_SIZE / np);
	}

	if(isMaster) pmesg(97, "Initiate particles\n");
	pcord_t * particles;
	pcord_t * particlesKeep;
	pcord_t * particlesSendUp;
	pcord_t * particlesSendDown;
	slice_data_t * sliceData;
	slice_data_t * sd;
	{
		// Assuming that there is a statistical chance that ALL particles
		// either go up, down or stay
		particles = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesKeep = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesSendUp = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);
		particlesSendDown = (pcord_t *) malloc(sizeof(pcord_t) * MAX_NO_PARTICLES);

		sliceData = (slice_data_t *) malloc(sizeof(slice_data_t) * np);
		sd = (slice_data_t *) malloc(sizeof(slice_data_t));

		sd->current = INIT_NO_PARTICLES;
		sd->keeping = 0;
		sd->sendingUp = 0;
		sd->sendingDown = 0;
		sd->wallCollisions = 0;
		sd->pressure = 0.0;
		sd->timing = 0.0;

		// The functions start at the same time
		// Multiplier as spacing to make sure time doesn't
		// tick and gives me and me+1 the same value
		srand(time(NULL) + (me * 3571));
		for(int p = 0; p < sd->current; p++)
		{
			particles[p].x = slice.x0 + (rand01() * (slice.x1 - slice.x0));
			particles[p].y = slice.y0 + (rand01() * (slice.y1 - slice.y0));
			float r = rand01() * MAX_INITIAL_VELOCITY;
			float theta = rand01() * 2 * PI;
			particles[p].vx = r * cos(theta);
			particles[p].vy = r * sin(theta);

			// TODO: Check if the particle is too far away from a wall/slice
			// limit to ever reach it at it's speed? This would eliminate
			// a lot of calculations later.

			// TODO: Hey, while we're at it - particles that are more than
			// MAX_INITIAL_VELOCITY * TIME_STEPS
			// units away from any wall can't affect the final
			// calculated pressure at all!
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

	// Vampir Trace code grouping
	int symdef_step;
	int symdef_step_outer;
	int symdef_step_particle_inner;
	int symdef_particle_walls;
	int symdef_step_send;
	int symdef_step_send_up;
	int symdef_step_send_down;
	int symdef_step_recv;
	int symdef_step_recv_up;
	int symdef_step_recv_down;
	int * counters;
	long long lnglngCounterHit = 0;
	{
		VT_funcdef("Step", VT_NOCLASS, &symdef_step);
		VT_funcdef("Step:Particle", VT_NOCLASS, &symdef_step_outer);
		VT_funcdef("Step:Particle:Inner", VT_NOCLASS, &symdef_step_particle_inner);
		VT_funcdef("Step:Particle:Walls", VT_NOCLASS, &symdef_particle_walls);
		VT_funcdef("Step:Send", VT_NOCLASS, &symdef_step_send);
		VT_funcdef("Step:Send:Up", VT_NOCLASS, &symdef_step_send_up);
		VT_funcdef("Step:Send:Down", VT_NOCLASS, &symdef_step_send_down);
		VT_funcdef("Step:Recv", VT_NOCLASS, &symdef_step_recv);
		VT_funcdef("Step:Recv:Up", VT_NOCLASS, &symdef_step_recv_up);
		VT_funcdef("Step:Recv:Down", VT_NOCLASS, &symdef_step_recv_down);

		counters = (int *) malloc(sizeof(int) * 6);
		int * uintBounds = (int *) malloc(sizeof(int) * 2);
		uintBounds[0] = 0;
		uintBounds[1] = INT_MAX;
		float * floatBounds = (float *) malloc(sizeof(float) * 2);
		floatBounds[0] = 0.0;
		floatBounds[1] = FLT_MAX;

		int uintDataType = VT_COUNT_INTEGER64 | VT_COUNT_RATE | VT_COUNT_VALID_SAMPLE;
		int floatDataType = VT_COUNT_FLOAT | VT_COUNT_ABSVAL | VT_COUNT_VALID_SAMPLE;

		VT_countdef("Wall collisions up", VT_NOCLASS, uintDataType, VT_ME, uintBounds, "hits", & counters[0]);
		VT_countdef("Wall collisions down", VT_NOCLASS, uintDataType, VT_ME, uintBounds, "hits", & counters[1]);
		VT_countdef("Wall collisions left", VT_NOCLASS, uintDataType, VT_ME, uintBounds, "hits", & counters[2]);
		VT_countdef("Wall collisions right", VT_NOCLASS, uintDataType, VT_ME, uintBounds, "hits", & counters[3]);
		VT_countdef("Particle collisions", VT_NOCLASS, uintDataType, VT_ME, uintBounds, "hits", & counters[4]);
		VT_countdef("Pressure", VT_NOCLASS, floatDataType, VT_ME, floatBounds, "unit", & counters[5]);
	}
	VT_end(0); // Init

	// Basic timing
	double starttime, endtime;
	starttime = MPI_Wtime();

	if(isMaster) pmesg(97, "Main loop: for each time-step do\n");
	for(int t = 0; t < TIME_STEPS; t += TIME_STEP)
	{
		VT_enter(symdef_step, VT_NOSCL);

		if(isMaster) pmesg(75, "Time is %d\n", t);

		sd->keeping = 0;
		sd->sendingUp = 0;
		sd->sendingDown = 0;

		// Outer particle loop
		VT_enter(symdef_step_outer, VT_NOSCL);
		for(int po = 0; po < sd->current; po++)
		{
			int collision = 0;

			// Inner particle loop
			VT_enter(symdef_step_particle_inner, VT_NOSCL);
			// Checks only the particles from po and down the line
			// to reduce number of checks.

			// TODO: If particles[] was sorted according to, for example, y position
			// it would mean a lot better peformanceas checks
			for(int pi = po+1; pi < sd->current; pi++)
			{
				float collision = collide(&particles[po], &particles[pi]);
				if(collision != -1)
				{
					interact(&particles[po], &particles[pi], collision);

					lnglngCounterHit++;
					VT_countval (1, & counters[4], & lnglngCounterHit);

					// Don't do further collision checks to cut calculations
					collision = 1;
					break;
				}
			}
			VT_end(0); // Inner

			if(!collision){
				feuler(&particles[po], TIME_STEP);
			}

			// Collisions
			VT_enter(symdef_particle_walls, VT_NOSCL);
			// Check left and right wall collisions
			// Done first, they may bounce out of slice bounds
			sd->pressure += wall_collide_leftright(&particles[po], walls, &sd->wallCollisions, counters);

			// Check top wall collisions only for the first slice
			// otherwise check for escaping particles
			int escaped = 0;
			if(isFirstSlice)
			{
				sd->pressure += wall_collide_top(&particles[po], walls, &sd->wallCollisions, counters);
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
				sd->pressure += wall_collide_bottom(&particles[po], walls, &sd->wallCollisions, counters);
			}
			else if(particle_escape_bottom(&particles[po], slice))
			{
				particlesSendDown[sd->sendingDown] = particles[po];
				sd->sendingDown++;
				escaped = 1;
			}

			if(escaped == 0)
			{
				particlesKeep[sd->keeping] = particles[po];
				sd->keeping++;
			}

			double dblPressure = (double) sd->pressure;
			VT_countval(1, & counters[5], & dblPressure);
			VT_end(0); // Walls
		}
		VT_end(0); // Outer

		// If this is the last time step, sending and receiving isn't really necessary,
		// it is done to keep the stats correct at the summary

		VT_enter(symdef_step_send, VT_NOSCL);
		VT_enter(symdef_step_send_up, VT_NOSCL);
		// Send up, but only if there's something to send
		if(!isFirstSlice && sd->sendingUp != 0)
		{
			MPI_Request request;
			MPI_Isend(particlesSendUp, sd->sendingUp, MPI_COORDINATE, sliceAbove, TAG_SEND_UP, MPI_COMM_WORLD, & request);
		}
		VT_end(0); // Send up

		VT_enter(symdef_step_send_down, VT_NOSCL);
		// Send down, but only if there's something to send
		if(!isLastSlice && sd->sendingDown != 0)
		{
			MPI_Request request;
			MPI_Isend(particlesSendDown, sd->sendingDown, MPI_COORDINATE, sliceBelow, TAG_SEND_DOWN, MPI_COMM_WORLD, & request);
		}
		VT_end(0); // Send down
		VT_end(0); // Send

		// All sends done, then check the receive buffers
		MPI_Barrier(MPI_COMM_WORLD);

		VT_enter(symdef_step_recv, VT_NOSCL);
		VT_enter(symdef_step_recv_up, VT_NOSCL);
		// Receive from below, but only if there's something to receive
		if(!isLastSlice){
			int flag;
			MPI_Status status;
			MPI_Iprobe(sliceBelow, TAG_SEND_UP, MPI_COMM_WORLD, & flag, & status);
			if(flag)
			{
				int count;
				MPI_Get_count(& status, MPI_COORDINATE, & count);
				MPI_Recv(particlesSendUp, count, MPI_COORDINATE, sliceBelow, TAG_SEND_UP, MPI_COMM_WORLD, &status);

				// TODO: Insert sorted according to x position if possible?
				// Would increase the chances of a collision with a neighbor
				// and since collisions cut calculations, it would be efficient
				for(int pr = 0; pr < count; pr++)
				{
					particlesKeep[sd->keeping] = particlesSendUp[pr];
					sd->keeping++;
				}
			}
		}
		VT_end(0); // Recv up

		VT_enter(symdef_step_recv_down, VT_NOSCL);
		// Receive from above, but only if there's something to receive
		if(!isFirstSlice){
			int flag;
			MPI_Status status;
			MPI_Iprobe(sliceAbove, TAG_SEND_DOWN, MPI_COMM_WORLD, & flag, & status);
			if(flag)
			{
				int count;
				MPI_Get_count(& status, MPI_COORDINATE, & count);
				MPI_Recv(particlesSendDown, count, MPI_COORDINATE, sliceAbove, TAG_SEND_DOWN, MPI_COMM_WORLD, & status);

				// TODO: Insert sorted according to x position if possible?
				// Would increase the chances of a collision with a neighbor
				// and since collisions cut calculations, it would be efficient
				for(int pr = 0; pr < count; pr++)
				{
					particlesKeep[sd->keeping] = particlesSendDown[pr];
					sd->keeping++;
				}
			}
		}
		VT_end(0); // Recv down
		VT_end(0); // Recv

		// Use the kept particles for the next iteration
		swap(int, sd->current, sd->keeping);
		swap(pcord_t *, particles, particlesKeep);
		VT_end(0); // Step
	}
	endtime = MPI_Wtime();
	sd->timing = endtime - starttime;

	// Swap back so that numbers look correct for the stats
	swap(int, sd->current, sd->keeping);

	// One last gather to calculate total pressure
	MPI_Gather(sd, 1, MPI_SLICE_DATA, sliceData, 1, MPI_SLICE_DATA, 0, MPI_COMM_WORLD);

	if(isMaster)
	{
		slice_data_t sums;
		sums.current = 0;
		sums.keeping = 0;
		sums.sendingUp = 0;
		sums.sendingDown = 0;
		sums.pressure = 0;
		sums.wallCollisions = 0;
		sums.timing = 0;

		pmesg(40, "Slice data when done\n");
		printf("ID \tCURR  \tKEEP  \tUP    \tDOWN  \tPRESSURE    \tWCOLLI\tTIMING\n");
		for(int i = 0; i < np; i++){
			printf("%3d\t%6d\t%6d\t%6d\t%6d\t%12.04f\t%6d\t%5.03f\n", i, sliceData[i].current, sliceData[i].keeping, sliceData[i].sendingUp, sliceData[i].sendingDown, sliceData[i].pressure, sliceData[i].wallCollisions, sliceData[i].timing);
			sums.current += sliceData[i].current;
			sums.keeping += sliceData[i].keeping;
			sums.sendingUp += sliceData[i].sendingUp;
			sums.sendingDown += sliceData[i].sendingDown;
			sums.pressure += sliceData[i].pressure;
			sums.wallCollisions += sliceData[i].wallCollisions;
			sums.timing += sliceData[i].timing;
		}
		printf("---\t------\t------\t------\t------\t------------\t------\t------\n");
		printf("SUM\t%6d\t%6d\t%6d\t%6d\t%12.04f\t%6d\t%5.03f\n", sums.current, sums.keeping, sums.sendingUp, sums.sendingDown, sums.pressure, sums.wallCollisions, sums.timing);
		printf("AVG\t%6d\t%6d\t%6d\t%6d\t%12.04f\t%6d\t%5.03f\n", sums.current/np, sums.keeping/np, sums.sendingUp/np, sums.sendingDown/np, sums.pressure/np, sums.wallCollisions/np, sums.timing/np);

		pmesg(70, "Calculate pressure\n");
		{
			float averagePressure = (sums.pressure / (WALL_LENGTH * TIME_STEPS));
			printf("Calculated pressure to be %.5f (%e) in a box with size %.1f by %.1f over %d time steps with %d particles.\n", averagePressure, averagePressure, BOX_HORIZ_SIZE, BOX_VERT_SIZE, TIME_STEPS, INIT_NO_PARTICLES*np);
		}
	}
	MPI_Finalize();
	VT_end(0); // Main
	//VT_leave(VT_NOSCL);
	VT_finalize();
}
