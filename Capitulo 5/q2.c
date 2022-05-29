/* 
 *
 * Compile: gcc -g -Wall -fopenmp -o q2 q2.c
 * Usage:   ./q2 <number of threads>

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);
double f(double x);	/* Function we're integrating */
void Trap(double a, double b, int n, double* global_result_p);
void PTrap(double a, double b, int n, double* global_result_p);

int main(int argc, char* argv[]) {
   double  global_result = 0.0;  /* Store result in global_result */
   double global_result2 = 0.0;
   double  a, b;             	/* Left and right endpoints  	*/
   int 	n;                	/* Total number of trapezoids	*/
   int 	thread_count;

   if (argc != 2) Usage(argv[0]);
   thread_count = strtol(argv[1], NULL, 10);
   printf("Enter a, b, and n\n");
   a=0;
   b=100;
   for(n = 10; n<=1000; n = n+10){
    #  pragma omp parallel num_threads(thread_count)
  	 Trap(a, b, n, &global_result);
    PTrap(a, b, n, &global_result2);
    if(global_result != global_result2){
   		printf("With n = %d trapezoids, we got different results\n", n);
  			printf("  %.20e  - %.20e = %.20e \n", global_result, global_result2, global_result-global_result2);    
    }
   }

   return 0;
}  /* main */

/*--------------------------------------------------------------------
 * Function:	Usage
 * Purpose: 	Print command line for function and terminate
 * In arg:  	prog_name
 */
void Usage(char* prog_name) {

   fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
   fprintf(stderr, "   number of trapezoids must be evenly divisible by\n");
   fprintf(stderr, "   number of threads\n");
   exit(0);
}  /* Usage */

/*------------------------------------------------------------------
 * Function:	f
 * Purpose: 	Compute value of function to be integrated
 * Input arg:   x
 * Return val:  f(x)
 */
double f(double x) {
   double return_val;

   return_val = x*x;
   return return_val;
}  /* f */

/*------------------------------------------------------------------
 * Function:	Trap
 * Purpose: 	Use trapezoidal rule to estimate definite integral
 * Input args:  
 *	a: left endpoint
 *	b: right endpoint
 *	n: number of trapezoids
 * Output arg:
 *	integral:  estimate of integral from a to b of f(x)
 */
void PTrap(double a, double b, int n, double* global_result_p) {
   double  h, x, my_result;
   double  local_a, local_b;
   int  i, local_n;
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();

   h = (b-a)/n;
   local_n = n/thread_count;  
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   my_result = (f(local_a) + f(local_b))/2.0;
   for (i = 1; i <= local_n-1; i++) {
 	x = local_a + i*h;
 	my_result += f(x);
   }
   my_result = my_result*h;

#  pragma omp critical
   *global_result_p += my_result;
}  /* Trap */
void Trap(double a, double b, int n, double* global_result_p) {
   double  h, x, my_result;
   double  local_a, local_b;
   int  i, local_n;
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();

   h = (b-a)/n;
   local_n = n/thread_count;  
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   my_result = (f(local_a) + f(local_b))/2.0;
   for (i = 1; i <= local_n-1; i++) {
 	x = local_a + i*h;
 	my_result += f(x);
   }
   my_result = my_result*h;

//#  pragma omp critical
   *global_result_p += my_result;
}  /* Trap */
