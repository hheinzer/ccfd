/*
 * output.c
 *
 * Created: Mon 23 Mar 2020 10:42:06 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <string.h>

#include "main.h"
#include "output.h"
#include "readInTools.h"

/* extern variables */
char strOutFile[STRLEN];
double IOtimeInterval;
int IOiterInterval;
int iVisuProg;
char parameterFile[STRLEN];
int resUnit;
bool doErrorOutput;
outputTime_t outputTimes;

/*
 * Initialize output
 */
void initOutput(void)
{
	printf("\nInitializing IO:\n");
	strcpy(strOutFile, getStr("fileName", NULL));
	IOtimeInterval = getDbl("IOtimeInterval", NULL);
	IOiterInterval = getDbl("IOiterInterval", NULL);
	iVisuProg = getInt("outputFormat", "1");
}
