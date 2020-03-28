/*
 * timeDiscretization.h
 *
 * Created: Sat 21 Mar 2020 07:48:34 PM CET
 * Author : hhh
 */

#ifndef TIMEDISCRETIZATION_H
#define TIMEDISCRETIZATION_H

#include <stdbool.h>

extern double	cfl;
extern double	dfl;
extern double	t;

extern double	timeOverall;
extern double	ioTime;

extern int	timeOrder;
extern bool	isTimeStep1D;

extern bool	isStationary;
extern long	maxIter;
extern double	stopTime;
extern long	iniIterationNumber;
extern double	startTime;
extern double	abortResidual;
extern int	abortVariable;
extern char	abortVariableName[4];
extern double	clAbortResidual, cdAbortResidual;
extern bool	doAbortOnClResidual, doAbortOnCdResidual;
extern bool	isRestart;
extern double	restartTime;

extern int	printIter;
extern double	printTime;

extern bool	isRK;
extern int	nRKstages;
extern double	RKcoeff[5];
extern bool	isImplicit;

void initTimeDisc(void);

#endif
