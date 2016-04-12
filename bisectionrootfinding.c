#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"bisection_root_finding.h"


void print_guess(float x, float y)
{
  /* print the guess to screen */
  /*printf("x = % e\ty = % e\n",x,y);*/
  printf("% e\t% e\n",x,y);
}
void found_root(float x, float y)
{
  /* print the root to screen */
  print_guess(x,y);
}

float bisection_root_finding(float (*f)(float x), float x_min_start, float x_max_start)
{
  /* function that uses bisection search to find a root */

  float x_min = x_min_start;  /* minimum x in bracket */
  float x_mid;		      /* 1/2 between x_min and x_max */
  float x_max = x_max_start;  /* maximum x in bracket */
  float y_min;		      /* f at x_min */
  float y_mid;		      /* f at x_mid */
  float y_max;		      /* f at x_max */

  int flag    = 0;	      /* a flag for triggering */
  int iter    = 0;	      /* current number of iterations */
  int imax    = 1000;	      /* maximum number of iterations */

  /* check our initial bracket to make sure it's valid */
  flag = check_initial_values(f,x_min, x_max);
  if(flag)
  {
    /* we made a lucky guess and we're done */
    if(flag==1)	    /*x_min is the root*/
      return x_min;
    if(flag==2)     /*x_max is the root*/
      return x_max;
  }

  /*instead, if flag==0, we continue on our way to find the root*/

  /* loop until we find the root or we reach the maximum number of iterations*/
  while(!flag)
  {
    /* find value midway between bracket values */
    x_mid = 0.5*(x_min + x_max);

    /* determine the function value at this point */
    y_mid = f(x_mid);

    /* if the absolute value of f(x) is < tolerance, then exit */
    if(fabs(y_mid) < TOLERANCE)
    {
      /* trigger an exit */
      flag = 1;

    }else{

      /*x_mid is not a root*/
  
      y_min = f(x_min);
      y_max = f(x_max);

      /* print guesses */
      print_guess(x_min,y_min);
      print_guess(x_mid,y_mid);
      print_guess(x_max,y_max);

      /* if the product of the function at the midpoint and at one of the endpoints */
      /* is greater than zero, replace this end point */
      if( f(x_min)*f(x_mid) > 0)
      {
	/*replace x_min with x_mid */
	x_min = x_mid;
      }else{
	/*replace x_max with x_mid */
	x_max = x_mid;
      }
      

      /* count the iteration */
      iter++;

      /* if we have exceeded the max number of iterations, exit*/
      if(iter>=imax)
      {
	printf("Exceeded %d iterations, exiting (x=%f, f=%f)...\n",imax,x_mid,y_mid);
	fflush(stdout);	
	exit(-1);
      }
    }
  }

  /* return x_mid == root */
  return x_mid;
  
}


int check_initial_values(float (*f)(float x), float x_min, float x_max)
{

  /*check to be sure our initial guesses are valid*/

  float y_min = f(x_min);
  float y_max = f(x_max);

  /*check for zero crossing*/
  if( y_min*y_max >= 0.0 )
  {
    printf("No zero crossing found in range [%f, %f].\n",x_min,x_max);
    printf("f(%f) = %f, f(%f) = %f\n",x_min,y_min,x_max,y_max);
    printf("Change initial values.\n");
    fflush(stdout);
    exit(-1);
  }

  /*If x_min is a root, then return a flag == 1*/
  if( fabs(y_min) < TOLERANCE)
    return 1;

  /*If x_max is a root, then return a flag == 2*/
  if( fabs(y_max) < TOLERANCE)
    return 2;

  /*if we reach this point, the bracket is valid.*/
  return 0;
}
