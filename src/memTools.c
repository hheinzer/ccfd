/*
 * memTools.c
 *
 * Created: Fri 27 Mar 2020 02:38:03 PM CET
 * Author : hhh
 *
 * Manual Memory Management:
 *
 * " The manual type involves malloc and free, and is where most of your segfaults
 * happen. This memory model is why Jesus weeps when he has to code in C. "
 *	- Ben Klemens
 */

#include <stdlib.h>
#include <string.h>

/*
 * allocate a dynamic 2D array of integers
 */
long **dyn2DintArray(long I, long J)
{
	long **arr = malloc(sizeof(long *) * I + sizeof(long) * I * J);
	long *ptr = (long *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = 0;
		}
	}
	return arr;
}

/*
 * allocate a dynamic 2D array of doubles
 */
double **dyn2DdblArray(long I, long J)
{
	double **arr = malloc(sizeof(double *) * I + sizeof(double) * I * J);
	double *ptr = (double *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = 0.0;
		}
	}
	return arr;
}

/*
 * allocate a dynamic 3D array of integers
 */
long ***dyn3DintArray(long I, long J, long K)
{
	long ***arr = malloc(sizeof(long *) * I + sizeof(long **) * I * J + sizeof(long) * I * J * K);
	long **ptrI = (long **)(arr + I);
	long *ptrJ = (long *)(arr + I + I * J);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
			for (long k = 0; k < K; ++k) {
				arr[i][j][k] = 0.0;
			}
		}
	}
	return arr;
}

/*
 * allocate a dynamic 3D array of doubles
 */
double ***dyn3DdblArray(long I, long J, long K)
{
	double ***arr = malloc(sizeof(double *) * I + sizeof(double **) * I * J + sizeof(double) * I * J * K);
	double **ptrI = (double **)(arr + I);
	double *ptrJ = (double *)(arr + I + I * J);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
			for (long k = 0; k < K; ++k) {
				arr[i][j][k] = 0.0;
			}
		}
	}
	return arr;
}

/*
 * allocate a dynamic 4D array of doubles
 */
double ****dyn4DdblArray(long I, long J, long K, long L)
{
	double ****arr = malloc(sizeof(double *) * I + sizeof(double **) * I * J + sizeof(double ***) * I * J * K + sizeof(double) * I * J * K * L);
	double ***ptrI = (double ***)(arr + I);
	double **ptrJ = (double **)(arr + I + I * J);
	double *ptrK = (double *)(arr + I + I * J + I * J * K);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
			for (long k = 0; k < K; ++k) {
				arr[i][j][k] = ptrK + J * K * L * i + K * L * j + L * k;
				for (long l = 0; l < L; ++l) {
					arr[i][j][k][l] = 0.0;
				}
			}
		}
	}
	return arr;
}

/*
 * allocate a dynamic array of strings
 */
char **dynStringArray(long I, long J)
{
	char **arr = malloc(sizeof(char *) * I + sizeof(char) * I * J);
	char *ptr = (char *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
		strcpy(arr[i], "");
	}
	return arr;
}
