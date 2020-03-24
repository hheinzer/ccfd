/*
 * equation.c
 *
 * Created: Tue 24 Mar 2020 08:30:28 AM CET
 * Author : hhh
 */

#include <stdio.h>
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
//double CFL;
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
	printf("| Using EULER equations\n");

	gamma = getDbl("gamma", "1.4");
	R = getDbl("R", "287");

	iFlux = getInt("fluxFunction", NULL);

	pi = acos(-1.0);
	gamma1 = gamma - 1.0;
	gamma2 = gamma - 2.0;
	gamma1q = 1.0 / gamma1;
	sqrt2 = sqrt(2.0);
	sqrt3 = sqrt(3.0);
	sqrt3q = 1 / sqrt3;

	mu = 0.0;

	Pr = 0.72;

	doCalcSource = getBool("calcSource", "F");
	if (doCalcSource) {
		sourceFunc = getInt("sourceFunction", NULL);
	}
}
