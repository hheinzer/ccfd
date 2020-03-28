/*
 * initialCondition.c
 *
 * Created: Fri 27 Mar 2020 02:31:24 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "main.h"
#include "readInTools.h"
#include "initialCondition.h"
#include "equation.h"
#include "memTools.h"
#include "timeDiscretization.h"
#include "mesh.h"
#include "exactFunction.h"

/* extern variables */
int icType;
int nDomains;
int *domainID;
double rp1Dinterface;
double alpha;
double **refState;

/*
 * get initial flow conditions from init file
 */
void initInitialCondition(void)
{
	printf("\nInitializing Initial Conditions:\n");
	icType = getInt("icType", NULL);
	switch (icType) {
	case 1:
		nDomains = getInt("nDomains", "1");
		domainID = malloc(nDomains * sizeof(int));
		refState = dyn2DdblArray(nDomains, NVAR);
		for (int i = 0; i < nDomains; ++i) {
			printf("| Domain No.%d:\n", i);

			domainID[i] = getInt("domainID", NULL);
			printf("|   Domain ID: %d\n", domainID[i]);

			refState[i][RHO] = getDbl("rho", NULL);

			double Ma = getDbl("mach", NULL);
			alpha = getDbl("alpha", NULL);

			refState[i][P] = getDbl("pressure", NULL);

			double c = sqrt(gamma * refState[i][P] / refState[i][RHO]);
			double v = Ma * c;
			refState[i][VX] = v * cos(alpha * pi / 180.0);
			refState[i][VY] = v * sin(alpha * pi / 180.0);
		}
		break;
	case 2:
		intExactFunc = getInt("exactFunc", NULL);
		if (intExactFunc == 5) {
			refState = dyn2DdblArray(2, NVAR);
			rp1Dinterface = getDbl("RP_1D_interface", "0.5");
			memcpy(refState[0], getDblArray("StateLeft", 4, NULL), 4 * sizeof(double));
			memcpy(refState[1], getDblArray("StateRight", 4, NULL), 4 * sizeof(double));
		}
		break;
	default:
		printf("| ERROR: Initial Condition '%d' not known\n", icType);
		exit(1);
	}
}

/*
 * set initial flow in all cells
 */
void setInitialCondition(void)
{
	printf("\nSetting Initial Conditions:\n");
	elem_t *aElem;
	if (isRestart) {
		// TODO: CGNS_readSolution();
	} else {
		switch (icType) {
		case 0:
			/* cell test */
			aElem = firstElem;
			while (aElem) {
				aElem->pVar[VX] = 0.0;
				aElem->pVar[VY] = 0.0;
				aElem->pVar[RHO] = 1.0 * aElem->id;
				aElem->pVar[P] = 1.0;

				aElem = aElem->next;
			}
			break;
		case 1:
			/* 1: mult-domain homogenous initial condition
			 * 2: homogenous initial condition over all domains */
			aElem = firstElem;
			while (aElem) {
				if (nDomains == 1) {
					memcpy(aElem->pVar, refState[0], 4 * sizeof(double));
				} else {
					memcpy(aElem->pVar, refState[aElem->domain], 4 * sizeof(double));
				}
				aElem = aElem->next;
			}
			break;
		case 2:
			/* exact function */
			aElem = firstElem;
			while (aElem) {
				exactFunc(intExactFunc, aElem->bary, 0.0, aElem->pVar);
				aElem = aElem->next;
			}
			break;
		default:
			printf("| ERROR: Illegal Initial Condition: %d\n", icType);
			exit(1);
		}
	}

	aElem = firstElem;
	while (aElem) {
		primCons(aElem->pVar, aElem->cVar);
		aElem = aElem->next;
	}

	printf("| Done.\n");
}
