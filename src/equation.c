/** \file
 *
 * \brief Contains the function for initializing the physical constants
 *
 * \date Tue 24 Mar 2020 08:30:28 AM CET
 * \author hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "equation.h"
#include "readInTools.h"

/* extern variables */
double pi;				/**< pi */

bool doCalcSource;			/**< calculate source flag */
double R;				/**< specific gas constant */
double gam;				/**< specific heat ratio */
double gam1;				/**< `gam` - 1 */
double gam2;				/**< `gam` - 2 */
double gam1q;				/**< 1.0 / (`gam` - 1) */
double cp;				/**< specific heat capacity */
double Pr;				/**< Prandtl number */
double mu;				/**< dynamic viscosity */

int iFlux;				/**< flux function control */

int intExactFunc;			/**< exact function control */
int sourceFunc;				/**< source function control */
double sqrt2;				/**< sqrt(2.0) */
double sqrt3;				/**< sqrt(3.0) */
double sqrt3q;				/**< 1.0 / sqrt(3.0) */

/**
 * \brief Initialize equations
 */
void initEquation(void)
{
	printf("\nInitializing Constants:\n");

	#ifdef NAVIERSTOKES
	printf("| Using NAVIER-STOKES equations\n");
	#endif

	#ifdef EULER
	printf("| Using EULER equations\n");
	#endif

	gam = getDbl("gamma", "1.4");
	R = getDbl("R", "287");

	#ifdef NAVIERSTOKES
	mu = getDbl("mu", "0.0");
	Pr = getDbl("Pr", "0.72");
	#endif

	#ifdef EULER
	mu = 0.0;
	Pr = 0.72;
	#endif

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
	gam1 = gam - 1.0;
	gam2 = gam - 2.0;
	gam1q = 1.0 / gam1;
	sqrt2 = sqrt(2.0);
	sqrt3 = sqrt(3.0);
	sqrt3q = 1 / sqrt3;

	doCalcSource = getBool("calcSource", "F");
	if (doCalcSource) {
		sourceFunc = getInt("sourceFunction", NULL);
	}
}
