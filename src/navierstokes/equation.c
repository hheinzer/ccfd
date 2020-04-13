/*
 * equation.c
 *
 * Created: Tue 24 Mar 2020 08:30:28 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "equation.h"
#include "readInTools.h"

/* extern variables */
double pi;

bool doCalcSource;
double R;
double gamma;
double gamma1;
double gamma2;
double gamma1q;
double cp;
double Pr;
double mu;

int iFlux;

int intExactFunc;
int sourceFunc;
double sqrt3q, sqrt3, sqrt2;

/*
 * Initialize equations
 */
void initEquation(void)
{
	printf("\nInitializing Constants:\n");
	printf("| Using NAVIER-STOKES equations\n");

	gamma = getDbl("gamma", "1.4");
	R = getDbl("R", "287");

	mu = getDbl("mu", "0.0");

	Pr = getDbl("Pr", "0.72");

	iFlux = getInt("fluxFunction", NULL);
	switch (iFlux) {
	case GOD:
		printf("| Flux Function: Godunov\n");
		break;
	case ROE:
		printf("| Flux Function: Roe\n");
		break;
	case HLL:
		printf("| Flux Function: HLL\n");
		break;
	case HLLE:
		printf("| Flux Function: HLLE\n");
		break;
	case HLLC:
		printf("| Flux Function: HLLC\n");
		break;
	case LXF:
		printf("| Flux Function: Lax-Friedrichs\n");
		break;
	case STW:
		printf("| Flux Function: Steger-Warming\n");
		break;
	case CEN:
		printf("| Flux Function: Central Discretization\n");
		break;
	case AUSMD:
		printf("| Flux Function: AUSMD\n");
		break;
	case AUSMDV:
		printf("| Flux Function: AUSMDV\n");
		break;
	case VANLEER:
		printf("| Flux Function: van Leer\n");
		break;
	default:
		printf("| ERROR: Unknown Flux Function:\n");
		printf("|  1: Godunov scheme\n");
		printf("|  2: Roe scheme\n");
		printf("|  3: HLL scheme\n");
		printf("|  4: HLLE scheme\n");
		printf("|  5: HLLC scheme\n");
		printf("|  6: Lax-Friedrichs scheme\n");
		printf("|  7: Steger-Warming scheme\n");
		printf("|  8: Central scheme\n");
		printf("|  9: AUSMD scheme\n");
		printf("| 10: AUSMDV scheme\n");
		printf("| 11: van Leer scheme\n");
		exit(1);
	}

	pi = acos(-1.0);
	gamma1 = gamma - 1.0;
	gamma2 = gamma - 2.0;
	gamma1q = 1.0 / gamma1;
	sqrt2 = sqrt(2.0);
	sqrt3 = sqrt(3.0);
	sqrt3q = 1 / sqrt3;

	doCalcSource = getBool("calcSource", "F");
	if (doCalcSource) {
		sourceFunc = getInt("sourceFunction", NULL);
	}
}
