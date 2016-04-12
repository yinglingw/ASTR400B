#ifndef BISECTION_H
#define BISECTION_H
#define TOLERANCE 1.0e-6  /*f(x) tolerance for root finding*/

/*check to be sure our initial guesses are valid*/
int   check_initial_values(float (*f)(float x), float x_min, float x_max);

/* function that uses bisection search to find a root */
float bisection_root_finding(float (*f)(float x), float x_min_start, float x_max_start);

/* print the root to screen */
void  found_root(float x, float y);

/* print the guess to screen */
void  print_guess(float x, float y);

#endif /*BISECTION_H*/
