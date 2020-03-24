/*
 * output.h
 *
 * Created: Mon 23 Mar 2020 10:37:31 PM CET
 * Author : hhh
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdbool.h>

extern char strOutFile[STRLEN];
extern double IOtimeInterval;
extern int IOiterInterval;
extern int iVisuProg;
extern char parameterFile[STRLEN];
extern int resUnit;
extern bool doErrorOutput;

typedef struct outputTime_t outputTime_t;
struct outputTime_t {
	long iter;
	double time;
	outputTime_t *next;
};
extern outputTime_t outputTimes;

void initOutput(void);

#endif
