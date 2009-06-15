#include<stdlib.h>
#include<math.h>

#include "coordinate.h"
#include "physics.h"

#ifndef _definitions_h
#define _definitions_h

#define PI 3.141592653

// Defined externally at compilation
//#define MAX_NO_PARTICLES  15000  /* Maximum number of particles/processor */
//#define INIT_NO_PARTICLES 10000    /* Initial number of particles/processor */

#define MAX_INITIAL_VELOCITY 50

#define BOX_HORIZ_SIZE 10000.0
#define BOX_VERT_SIZE 10000.0
#define WALL_LENGTH (2.0*BOX_HORIZ_SIZE+2.0*BOX_VERT_SIZE)

#define TIME_STEPS 10
#define TIME_STEP 1

struct slice_data {
	unsigned int current;
	unsigned int keeping;
	unsigned int sendingUp;
	unsigned int sendingDown;
	unsigned int wallCollisions;
	float pressure;
	double timing;
};

typedef struct slice_data slice_data_t;

#define TAG_SEND_UP 666
#define TAG_SEND_DOWN 999

#define swap(type, i, j) {type t = i; i = j; j = t;}

#endif
