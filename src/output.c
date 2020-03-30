/*
 * output.c
 *
 * Created: Mon 23 Mar 2020 10:42:06 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "main.h"
#include "output.h"
#include "readInTools.h"
#include "timeDiscretization.h"
#include "analyze.h"
#include "equation.h"
#include "exactFunction.h"

/* extern variables */
char strOutFile[STRLEN];
double IOtimeInterval;
int IOiterInterval;
int iVisuProg;
char parameterFile[STRLEN];
FILE *resFile;
bool doErrorOutput;
outputTime_t *outputTimes;

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

/*
 * tabular CSV output, only for 1D data
 */
void csvOutput(char fileName[STRLEN], double time, double iter, bool doExact)
{
	/* prepare data (only for equidistant grids) */
	double flowData[nElems][NVAR];
	if (doExact) {
		long iElem = 0;
		elem_t *aElem = firstElem;
		double pVar[NVAR];
		while (aElem) {
			exactFunc(intExactFunc, aElem->bary, time, pVar);
			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = pVar[RHO];
			flowData[iElem][2] = pVar[VX];
			flowData[iElem][3] = pVar[P];
			iElem++;
			aElem = aElem->next;
		}
	} else {
		long iElem = 0;
		elem_t *aElem = firstElem;
		while (aElem) {
			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = aElem->pVar[RHO];
			flowData[iElem][2] = aElem->pVar[VX];
			flowData[iElem][3] = aElem->pVar[P];
			iElem++;
			aElem = aElem->next;
		}
	}

	/* sort array */
	long n = nElems;
	bool isSwapped = true;
	double tmp[NVAR];
	while (isSwapped && (n > 0)) {
		isSwapped = false;
		for (long iElem = 0; iElem < n - 1; ++iElem) {
			if (flowData[iElem][0] > flowData[iElem + 1][0]) {
				memcpy(tmp, flowData[iElem + 1], NVAR * sizeof(double));
				memcpy(flowData[iElem + 1], flowData[iElem], NVAR * sizeof(double));
				memcpy(flowData[iElem], tmp, NVAR * sizeof(double));
				isSwapped = true;
			}
		}
		n--;
	}

	/* write data */
	FILE *csvFile = fopen(fileName, "w");
	fprintf(csvFile, "CoordinateX, Density, Velocity, Pressure\n");
	for (long iElem = 0; iElem < nElems; ++iElem) {
		fprintf(csvFile, "%15.9f,%15.9f,%15.9f,%15.9f\n", flowData[iElem][0],
			flowData[iElem][1], flowData[iElem][2], flowData[iElem][3]);
	}
	fclose(csvFile);
}

/*
 * data output in different formats
 */
void dataOutput(double time, long iter)
{
	/* output times */
	outputTime_t *outputTime = malloc(sizeof(outputTime_t));
	outputTime->time = time;
	outputTime->iter = iter;
	outputTime->next = outputTimes;
	outputTimes = outputTime;

	/* write flow solution */
	char fileName[2 * STRLEN];
	if (isStationary) {
		sprintf(fileName, "%s_%09ld", strOutFile, iter);
	} else {
		sprintf(fileName, "%s_%015.7f", strOutFile, time);
	}
	switch (iVisuProg) {
	case CGNS:
		strcat(fileName, ".cgns");
		//cgnsOutput(fileName, time, iter, false);
		break;
	case CURVE:
		strcat(fileName, ".curve");
		//curveOutput(fileName, time, iter, false);
		break;
	case DAT:
		strcat(fileName, ".csv");
		csvOutput(fileName, time, iter, false);
		break;
	default:
		printf("| ERROR: Output Format unknown\n");
		exit(1);
	}

	/* write exact solution, if applicable */
	if (hasExactSolution) {
		if (isStationary) {
			sprintf(fileName, "%s_ex_%09ld", strOutFile, iter);
		} else {
			sprintf(fileName, "%s_ex_%015.7f", strOutFile, time);
		}
		switch (iVisuProg) {
		case CGNS:
			strcat(fileName, ".cgns");
			//cgnsOutput(fileName, time, iter, true);
			break;
		case CURVE:
			strcat(fileName, ".curve");
			//curveOutput(fileName, time, iter, true);
			break;
		case DAT:
			strcat(fileName, ".csv");
			csvOutput(fileName, time, iter, true);
			break;
		default:
			printf("| ERROR: Output Format unknown\n");
			exit(1);
		}
	}
}
