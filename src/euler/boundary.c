/*
 * boundary.c
 *
 * Created: Tue 24 Mar 2020 10:10:51 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "boundary.h"
#include "readInTools.h"
#include "equation.h"

/* extern variables */
boundary_t *firstBC;
int nBC;
bool isPeriodic;

void initBoundary(void)
{
	printf("\nInitializing Boundary Conditions:\n");

	nBC = getInt("nBC", NULL);

	firstBC = NULL;

	for (int iBC = 0; iBC < nBC; ++iBC) {
		boundary_t *aBC = malloc(sizeof(boundary_t));

		aBC->next = firstBC;
		firstBC = aBC;

		int intIn = getInt("BCtype", NULL);
		aBC->BCtype = intIn / 100;
		aBC->BCid = intIn % 100;

		double Ma, alpha, c, v;
		switch (aBC->BCtype) {
		case SLIPWALL:
			printf("| BC Type: Slipwall\n");
			break;
		case INFLOW:
			printf("| BC Type: Inflow\n");

			aBC->pVar[RHO] = getDbl("rho", NULL);

			Ma = getDbl("mach", NULL);

			alpha = getDbl("alpha", NULL);

			aBC->pVar[P] = getDbl("pressure", NULL);

			c = sqrt(gamma * aBC->pVar[P] / aBC->pVar[RHO]);
			v = Ma * c;
			aBC->pVar[VX] = v * cos(alpha * pi / 180.0);
			aBC->pVar[VY] = v * sin(alpha * pi / 180.0);
			break;
		case CHARACTERISTIC:
			printf("| BC Type: Characteristic\n");

			aBC->pVar[RHO] = getDbl("rho", NULL);

			Ma = getDbl("mach", NULL);

			alpha = getDbl("alpha", NULL);

			aBC->pVar[P] = getDbl("pressure", NULL);

			c = sqrt(gamma * aBC->pVar[P] / aBC->pVar[RHO]);
			v = Ma * c;
			aBC->pVar[VX] = v * cos(alpha * pi / 180.0);
			aBC->pVar[VY] = v * sin(alpha * pi / 180.0);
			break;
		case OUTFLOW:
			printf("| BC Type: Outflow\n");
			break;
		case EXACTSOL:
			printf("| BC Type: Exact Function\n");

			aBC->exactFunc = getInt("BCexactFunc", NULL);
			break;
		case PERIODIC:
			printf("| BC Type: Periodic\n");

			aBC->connection = getDblArray("connection", NDIM, NULL);
			break;
		default:
			printf("| ERROR: Illegal boundary condition!\n");
			exit(1);
		}
	}
}
