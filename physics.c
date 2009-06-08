#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "physics.h"


#ifndef sqr
#define sqr(a) ((a)*(a))
#endif

#ifndef sign
#define sign(a) ((a) > 0 ? 1 : -1)
#endif

int feuler(pcord_t *a, float time_step){
    a->x = a->x + time_step* a->vx ;
    a->y = a->y + time_step* a->vy ;	    
    return 0 ;
}

float wall_collide(pcord_t *p, cord_t wall){
    float gPreassure = 0.0 ;
    
    if(p->x < wall.x0){
	p->vx = -p->vx ;
	p->x  = wall.x0 + (wall.x0-p->x);
	gPreassure += 2.0*fabs(p->vx) ;
    }
    if(p->x > wall.x1){
	p->vx = -p->vx ;
	p->x  = wall.x1 - (p->x-wall.x1);
	gPreassure += 2.0*fabs(p->vx) ;
    }
    if(p->y < wall.y0){
	p->vy = -p->vy ;
	p->y  = wall.y0 + (wall.y0-p->y);	
	gPreassure += 2.0*fabs(p->vy) ;
    }
    if(p->y > wall.y1){
	p->vy = -p->vy ;
	p->y  = wall.y1 - (p->y-wall.y1);
	gPreassure += 2.0*fabs(p->vy) ;
    }
    return gPreassure ;
}



float collide(pcord_t *p1, pcord_t *p2){
    double a,b,c;
    double temp,t1,t2;

    a=sqr(p1->vx-p2->vx)+sqr(p1->vy-p2->vy);
    b=2*((p1->x - p2->x)*(p1->vx - p2->vx)+(p1->y - p2->y)*(p1->vy - p2->vy));
    c=sqr(p1->x-p2->x)+sqr(p1->y-p2->y)-4*1*1;
    
    if (a!=0.0){	
	temp=sqr(b)-4*a*c;
	if (temp>=0){
	    temp=sqrt(temp);
	    t1=(-b+temp)/(2*a);
	    t2=(-b-temp)/(2*a);
      
	    if (t1>t2){
		temp=t1;
		t1=t2;
		t2=temp;
	    }
	    if ((t1>=0)&(t1<=1))
		return t1;
	    else if ((t2>=0)&(t2<=1))
		return t2;    
	}
    }	
    return -1;
}



void interact(pcord_t *p1,pcord_t *p2, float t){
    float c,s,a,b,tao;
    pcord_t p1temp,p2temp;
  
    if (t>=0){

	/* Move to impact point */
	(void)feuler(p1,t);
	(void)feuler(p2,t);
    
	/* Rotate the coordinate system around p1*/
	p2temp.x=p2->x-p1->x;
	p2temp.y=p2->y-p1->y;
    
	/* Givens plane rotation, Golub, van Loan p. 216 */
	a=p2temp.x;
	b=p2temp.y;
	if (p2->y==0){
	    c=1;s=0;
	}
	else{
	    if (fabs(b)>fabs(a)){
		tao=-a/b;
		s=1/(sqrt(1+sqr(tao)));
		c=s*tao;
	    }
	    else{
		tao=-b/a;
		c=1/(sqrt(1+sqr(tao)));
		s=c*tao;
	    }
	}
    
	p2temp.x=c * p2temp.x+s * p2temp.y; /* This should be equal to 2r */
	p2temp.y=0.0;
    
	p2temp.vx= c* p2->vx + s* p2->vy;
	p2temp.vy=-s* p2->vx + c* p2->vy;
	p1temp.vx= c* p1->vx + s* p1->vy;
	p1temp.vy=-s* p1->vx + c* p1->vy;
    
	/* Assume the balls has the same mass... */
	p1temp.vx=-p1temp.vx;
	p2temp.vx=-p2temp.vx;
    
	p1->vx = c * p1temp.vx - s * p1temp.vy;
	p1->vy = s * p1temp.vx + c * p1temp.vy;
	p2->vx = c * p2temp.vx - s * p2temp.vy;
	p2->vy = s * p2temp.vx + c * p2temp.vy;

	/* Move the balls the remaining time. */
	c=1.0-t;
	(void)feuler(p1,c);
	(void)feuler(p2,c);
    }

}









