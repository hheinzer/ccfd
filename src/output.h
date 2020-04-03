/*
 * output.h
 *
 * Created: Mon 23 Mar 2020 10:37:31 PM CET
 * Author : hhh
 */

#ifndef OUTPUT_H
#define OUTPUT_H

typedef struct outputTime_t outputTime_t;

#include <stdbool.h>
#include <stdlib.h>

extern char strOutFile[STRLEN];
extern double IOtimeInterval;
extern int IOiterInterval;
extern int iVisuProg;
extern char parameterFile[STRLEN];
extern FILE *resFile;
extern bool doErrorOutput;

struct outputTime_t {
	long iter;
	double time;
	outputTime_t *next;
};
extern outputTime_t *outputTimes;

void initOutput(void);
void dataOutput(double time, long iter);
void finalizeDataOutput(void);

#endif
