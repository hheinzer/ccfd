/*
 * timeDiscretization.c
 *
 * Created: Sat 21 Mar 2020 07:52:42 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "timeDiscretization.h"
#include "readInTools.h"
#include "output.h"

/* extern variables */
double	cfl;
double	dfl;
double	t;

double	timeOverall;
double	ioTime;

int	timeOrder;
bool	isTimeStep1D;

bool	isStationary;
long	maxIter;
double	stopTime;
long	iniIterationNumber;
double	startTime;
double	abortResidual;
int	abortVariable;
char	abortVariableName[4];
double	clAbortResidual, cdAbortResidual;
bool	doAbortOnClResidual, doAbortOnCdResidual;
bool	isRestart;
double	restartTime;

int	printIter;
double	printTime;

bool	isRK;
int	nRKstages;
double	RKcoeff[5] = {0.0};
bool	isImplicit;

/*
 * Initialize the time discretization
 */
void initTimeDisc(void)
{
	printf("\nTime Discretization:\n");
	cfl = getDbl("CFL", "0.9");
	dfl = getDbl("DFL", "0.9");
	isTimeStep1D = getBool("timeStep1D", "F");
	timeOrder = getInt("timeOrder", "1");
	if (timeOrder > 3) {
		printf("| ERROR: Temporal Discretization Order must be either 1, 2 or 3\n");
		exit(1);
	}

	isImplicit = getBool("implicit", "F");

	nRKstages = getInt("nRKstages", "1");
	if (nRKstages > 5) {
		printf("| ERROR: No more than 5 RK Stages are Possible\n");
		exit(1);
	}
	if ((nRKstages == 1) && (timeOrder == 1)) {
		printf("| Time Stepping Scheme: Euler Time Stepping\n");
		RKcoeff[0] = 1.0;
	} else {
		printf("| Time Stepping Scheme: RK Time Stepping\n");
		switch (nRKstages) {
		case 3:
			switch (timeOrder) {
			case 1:
				RKcoeff[0] = 0.1418;
				RKcoeff[1] = 0.4;
				RKcoeff[2] = 1.0;
				break;
			case 2:
				RKcoeff[0] = 0.1918;
				RKcoeff[1] = 0.4929;
				RKcoeff[2] = 1.0;
				break;
			case 3:
				RKcoeff[0] = 0.2884;
				RKcoeff[1] = 0.5010;
				RKcoeff[2] = 1.0;
				break;
			}
			break;
		case 4:
			switch (timeOrder) {
			case 1:
				RKcoeff[0] = 0.0833;
				RKcoeff[1] = 0.2069;
				RKcoeff[2] = 0.4265;
				RKcoeff[3] = 1.0;
				break;
			case 2:
				RKcoeff[0] = 0.1084;
				RKcoeff[1] = 0.2602;
				RKcoeff[2] = 0.5052;
				RKcoeff[3] = 1.0;
				break;
			case 3:
				RKcoeff[0] = 0.1666;
				RKcoeff[1] = 0.3027;
				RKcoeff[2] = 0.5275;
				RKcoeff[3] = 1.0;
				break;
			}
			break;
		case 5:
			switch (timeOrder) {
			case 1:
				RKcoeff[0] = 0.0533;
				RKcoeff[1] = 0.1263;
				RKcoeff[2] = 0.2375;
				RKcoeff[3] = 0.4414;
				RKcoeff[4] = 1.0;
				break;
			case 2:
				RKcoeff[0] = 0.0695;
				RKcoeff[1] = 0.1602;
				RKcoeff[2] = 0.2898;
				RKcoeff[3] = 0.5060;
				RKcoeff[4] = 1.0;
				break;
			case 3:
				RKcoeff[0] = 0.1067;
				RKcoeff[1] = 0.1979;
				RKcoeff[2] = 0.3232;
				RKcoeff[3] = 0.5201;
				RKcoeff[4] = 1.0;
				break;
			}
			break;
		default:
			printf("| ERROR: Wrong Number of RK Coefficients\n");
			exit(1);
		}

		char line[STRLEN] = "";
		sprintf(line, "{%g", RKcoeff[0]);
		for (int iRK = 1; iRK < nRKstages; ++iRK) {
			sprintf(line + strlen(line), ", %g", RKcoeff[iRK]);
		}
		sprintf(line + strlen(line), "}");
		printf("| %16s[%d] = %27s\n", "RKcoeff", nRKstages, line);
	}

	/* stationary computation */
	if (isStationary) {
		doAbortOnClResidual = doAbortOnCdResidual = false;

		abortResidual = getDbl("abortResidual", "1e-6");

		abortVariable = getInt("abortVariable", "1");
		switch (abortVariable) {
		case RHO + 1:
			strcpy(abortVariableName, "RHO");
			break;
		case MX + 1:
			strcpy(abortVariableName, "MX");
			break;
		case MY + 1:
			strcpy(abortVariableName, "MY");
			break;
		case E + 1:
			strcpy(abortVariableName, "E");
			break;
		}

		clAbortResidual = getDbl("cl_abortResidual", "0.0");
		cdAbortResidual = getDbl("cd_abortResidual", "0.0");

		printf("| Stationary Problem\n");
		if ((clAbortResidual == 0.0) && (cdAbortResidual == 0)) {
			printf("| Abort Residual: %g\n", abortResidual);
		} else if (clAbortResidual > 0.0) {
			printf("| Cl Abort Residual: %g\n", clAbortResidual);
			doAbortOnClResidual = true;
		} else if (cdAbortResidual > 0.0) {
			printf("| Cd Abort Residual: %g\n", cdAbortResidual);
			doAbortOnCdResidual = true;
		} else {
			printf("| ERROR: Wrong Definition of Abort Residual\n");
			exit(1);
		}
	} else {
		printf("| Transient Problem\n");
	}

	maxIter = getInt("maxIter", "100000");
	stopTime = getDbl("tEnd", NULL);

	if (isRestart) {
		t = restartTime;
	} else {
		t = 0.0;
	}

	if (isStationary) {
		printf("| Initial Iteration Number: %ld\n", iniIterationNumber);
	} else {
		printf("| Start Time: %g\n", t);
	}

	printIter = (iniIterationNumber / IOiterInterval + 1) * IOiterInterval;
	printTime = (floor(t / IOtimeInterval) + 1) * IOtimeInterval;
}
