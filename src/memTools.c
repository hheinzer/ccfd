/** \file
 *
 * \brief Memory management functions
 *
 * \author hhh
 * \date Fri 27 Mar 2020 02:38:03 PM CET
 *
 * \note
 * > Manual Memory Management:
 * >
 * > "The manual type involves malloc and free, and is where most of your
 * > segfaults happen. This memory model is why Jesus weeps when he has to
 * > code in C."
 * >
 * > \- Ben Klemens
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cgnslib.h"

/** \brief Allocate a dynamic 2D array of integers
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \return Pointer to a 2D integer array
 */
long **dyn2DintArray(long I, long J)
{
	long **arr = calloc(1, sizeof(long *) * I + sizeof(long) * I * J);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	long *ptr = (long *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
	}
	return arr;
}

/** \brief Allocate a dynamic 2D array of cgsize_t
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \return Pointer to a 2D cgsize_t array
 */
cgsize_t **dyn2DcgsizeArray(long I, long J)
{
	cgsize_t **arr = calloc(1, sizeof(cgsize_t *) * I + sizeof(cgsize_t) * I * J);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	cgsize_t *ptr = (cgsize_t *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
	}
	return arr;
}

/** \brief Allocate a dynamic 2D array of doubles
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \return Pointer to a 2D double array
 */
double **dyn2DdblArray(long I, long J)
{
	double **arr = calloc(1, sizeof(double *) * I + sizeof(double) * I * J);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	double *ptr = (double *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
	}
	return arr;
}

/** \brief Allocate a dynamic 3D array of integers
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \param[in] K Number of elements in the third dimension
 * \return Pointer to a 3D integer array
 */
long ***dyn3DintArray(long I, long J, long K)
{
	long ***arr = calloc(1, sizeof(long *) * I + sizeof(long **) * I * J + sizeof(long) * I * J * K);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	long **ptrI = (long **)(arr + I);
	long *ptrJ = (long *)(arr + I + I * J);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
		}
	}
	return arr;
}

/** \brief Allocate a dynamic 3D array of doubles
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \param[in] K Number of elements in the third dimension
 * \return Pointer to a 3D double array
 */
double ***dyn3DdblArray(long I, long J, long K)
{
	double ***arr = calloc(1, sizeof(double *) * I + sizeof(double **) * I * J + sizeof(double) * I * J * K);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	double **ptrI = (double **)(arr + I);
	double *ptrJ = (double *)(arr + I + I * J);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
		}
	}
	return arr;
}

/** \brief Allocate a dynamic 4D array of doubles
 * \param[in] I Number of elements in the first dimension
 * \param[in] J Number of elements in the second dimension
 * \param[in] K Number of elements in the third dimension
 * \param[in] L Number of elements in the fourth dimension
 * \return Pointer to a 4D double array
 */
double ****dyn4DdblArray(long I, long J, long K, long L)
{
	double ****arr = calloc(1, sizeof(double *) * I + sizeof(double **) * I * J + sizeof(double ***) * I * J * K + sizeof(double) * I * J * K * L);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	double ***ptrI = (double ***)(arr + I);
	double **ptrJ = (double **)(arr + I + I * J);
	double *ptrK = (double *)(arr + I + I * J + I * J * K);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptrI + J * i;
		for (long j = 0; j < J; ++j) {
			arr[i][j] = ptrJ + J * K * i + K * j;
			for (long k = 0; k < K; ++k) {
				arr[i][j][k] = ptrK + J * K * L * i + K * L * j + L * k;
			}
		}
	}
	return arr;
}

/** \brief Allocate a dynamic array of strings
 * \param[in] I Number of elements in the first dimension
 * \param[in] J String length of each element
 * \return Pointer to a string array
 */
char **dynStringArray(long I, long J)
{
	char **arr = calloc(1, sizeof(char *) * I + sizeof(char) * I * J);
	if (!arr) {
		printf("| ERROR: could not allocate arr\n");
		exit(1);
	}

	char *ptr = (char *)(arr + I);
	for (long i = 0; i < I; ++i) {
		arr[i] = ptr + J * i;
	}
	return arr;
}
