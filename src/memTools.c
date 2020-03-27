/*
 * memTools.c
 *
 * Created: Fri 27 Mar 2020 02:38:03 PM CET
 * Author : hhh
 */

#include <stdlib.h>

/*
 * allocate a dynamic 2D array of integers
 */
long **dyn2DintArray(long I, int J)
{
	long **arr = malloc(sizeof(long *) * I + sizeof(long) * I * J);
	long *ptr = (long *)(arr + I);
	for (int i = 0; i < I; ++i) {
		arr[i] = (ptr + J * i);
	}
	return arr;
}

/*
 * allocate a dynamic 2D array of doubles
 */
double **dyn2DdblArray(long I, int J)
{
	double **arr = malloc(sizeof(double *) * I + sizeof(double) * I * J);
	double *ptr = (double *)(arr + I);
	for (int i = 0; i < I; ++i) {
		arr[i] = (ptr + J * i);
	}
	return arr;
}
