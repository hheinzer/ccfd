/** \file
 *
 * \brief Contains `outputTime_t` struct definition
 *
 * \author hhh
 * \date Mon 23 Mar 2020 10:37:31 PM CET
 */

#ifndef OUTPUT_H
#define OUTPUT_H

typedef struct outputTime_t outputTime_t;

#include <stdbool.h>
#include <stdlib.h>

/**
 * \brief Output times linked list
 */
struct outputTime_t {
	long iter;		/**< iteration number at output */
	double time;		/**< computational time at output */
	outputTime_t *next;	/**< pointer to next output time */
};

extern char strOutFile[STRLEN];
extern double IOtimeInterval;
extern int IOiterInterval;
extern int iVisuProg;
extern char parameterFile[STRLEN];
extern FILE *resFile;
extern bool doErrorOutput;
extern outputTime_t *outputTimes;

void initOutput(void);
void dataOutput(double time, long iter);
void finalizeDataOutput(void);
void cgnsWriteMesh(void);
void freeOutputTimes(void);

#endif
