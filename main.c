#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"bisection_root_finding.h"

/*function we want to find the roots for*/
float function_for_roots(float x);


int main(int argc, char **argv)
{
  float x_min;  /*minimum of bracket*/
  float x_max;  /*maximim of bracket*/

  float x_root; /*actual found root*/
  float y_root; /*actual function value at root == ~zero*/

  /* check number of command line arguments */
  if(argc!=3)
  {
    printf("Usage:\n");
    printf("./bisection x_min x_max\n");
    fflush(stdout);
    exit(-1);
  }

  /* read in x_min and x_max */
  x_min = atof(argv[1]);
  x_max = atof(argv[2]);


  /* perform bisection searching */
  x_root = bisection_root_finding( function_for_roots, x_min, x_max);

  /* get the actual y value == should be ~zero*/
  y_root = function_for_roots(x_root);

  /* output the result*/
  found_root(x_root, y_root);

  /*we're done*/
  return 0;
}


/*function we want to find the root for*/
float function_for_roots(float x)
{

  /*float a =  1.0;
  float b = -3.0;
  float c =  2.0;*/
  float a =  1.01;
  float b = -3.04;
  float c =  2.07;

  return a*x*x + b*x + c;
}
