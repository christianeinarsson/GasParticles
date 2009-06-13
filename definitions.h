#include<stdlib.h>
#include<math.h>

#include "coordinate.h"
#include "physics.h"

#ifndef _definitions_h
#define _definitions_h

#define PI 3.141592653

// Defined externally
//#define MAX_NO_PARTICLES  15000  /* Maximum number of particles/processor */
//#define INIT_NO_PARTICLES 500    /* Initial number of particles/processor */

#define MAX_INITIAL_VELOCITY 50

#define BOX_HORIZ_SIZE 10000.0
#define BOX_VERT_SIZE 10000.0
#define WALL_LENGTH (2.0*BOX_HORIZ_SIZE+2.0*BOX_VERT_SIZE)

#define PARTICLE_BUFFER_SIZE MAX_NO_PARTICLES/5
#define COMM_BUFFER_SIZE  5*PARTICLE_BUFFER_SIZE

struct particle {
	pcord_t  pcord;
	int ptype;        /* Used to simulate mixing of gases */
};

typedef struct particle particle_t;

#define TIME_STEPS 10
#define TIME_STEP 1

struct slice_data {
	unsigned int current;
	unsigned int keeping;
	unsigned int sendingUp;
	unsigned int sendingDown;
	unsigned int wallCollisions;
	float pressure;
};

typedef struct slice_data slice_data_t;

#define TAG_SEND_UP 666
#define TAG_SEND_DOWN 999

#define swap(type, i, j) {type t = i; i = j; j = t;}

#endif
