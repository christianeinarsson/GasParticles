#ifndef _physics_h
#define _physics_h

#include "coordinate.h"

#define STEP_SIZE 1.0 /* the step size use in the integration */

int feuler(pcord_t *a,
		   float time_step);

float wall_collide(pcord_t *p, cord_t wall);

float wall_collide_top(pcord_t *p, cord_t wall);
float wall_collide_bottom(pcord_t *p, cord_t wall);
float wall_collide_leftright(pcord_t *p, cord_t wall);
int particle_escape_top(pcord_t *p, cord_t limits);
int particle_escape_bottom(pcord_t *p, cord_t limits);


float collide(pcord_t *p1, pcord_t *p2);

void interact(pcord_t *p1, pcord_t *p2, float t);

double rand01();

void printParticles(int msglevel, int count, pcord_t particles[], int rank, const char * name);

#endif
