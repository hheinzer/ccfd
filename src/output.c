/*
 * output.c
 *
 * Created: Mon 23 Mar 2020 10:42:06 PM CET
 * Author : hhh
 */

#include <stdbool.h>

#include "main.h"
#include "output.h"

/* extern variables */
char strOutfile[STRLEN];
double IOtimeInterval;
int IOiterInterval;
int iVisuProg;
char parameterFile[STRLEN];
int resUnit;
bool doErrorOutput;
outputTime_t outputTimes;


