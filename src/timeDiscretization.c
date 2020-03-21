/*
 * timeDiscretization.c
 *
 * Created: Sat 21 Mar 2020 07:52:42 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "timeDiscretization.h"
#include "readInTools.h"

/* extern variables */
double	CFL;
double	DFL;
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
char	abortVariableName[3];
double	clAbortResidual, cdAbortResidual;
bool	doAbortOnClResidual, doAbortOnCdResidual;
bool	isRestart;
double	restartTime;

int	printInterval;
double	printTime;

bool	isRK;
int	nRKstages;
double	RKcoeff[6];
bool	isImplicit;

/*
 * Initialize the time discretization
 */
void initTimeDisc(void)
{
	printf("\nTime Discretization:\n");
}
