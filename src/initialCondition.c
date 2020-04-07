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
#include "equationOfState.h"
#include "cgnslib.h"

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
 * read solution from CGNS file (for restart)
 */
void cgnsReadSolution(void)
{
	int indexFile;
	/* open CGNS file */
	if (cg_open(strIniCondFile, CG_MODE_READ, &indexFile))
		cg_error_exit();

	/* get number of elements */
	char zoneName[32];
	cgsize_t iSize[3];
	if (cg_zone_read(indexFile, 1, 1, zoneName, iSize))
		cg_error_exit();

	if (nElems != iSize[1]) {
		printf("| ERROR: Wrong Number of Elements in CGNS flow solution\n");
		exit(1);
	}

	/* allocate array for the flow solution */
	double rhoArr[nElems], vxArr[nElems], vyArr[nElems], pArr[nElems];
	cgsize_t rMin[1] = {1}, rMax[1] = {nElems};
	if (cg_field_read(indexFile, 1, 1, 1, "Density", RealDouble, rMin, rMax, rhoArr))
		cg_error_exit();
	if (cg_field_read(indexFile, 1, 1, 1, "VelocityX", RealDouble, rMin, rMax, vxArr))
		cg_error_exit();
	if (cg_field_read(indexFile, 1, 1, 1, "VelocityY", RealDouble, rMin, rMax, vyArr))
		cg_error_exit();
	if (cg_field_read(indexFile, 1, 1, 1, "Pressure", RealDouble, rMin, rMax, pArr))
		cg_error_exit();

	/* read iteration number, time, and wall clock time */
	if (cg_goto(indexFile, 1, "end"))
		cg_error_exit();

	char descriptorName[32], *text;
	if (cg_descriptor_read(1, descriptorName, &text))
		cg_error_exit();
	sscanf(text, "%lg %lg", &t, &timeOverall);

	/* close CGNS file */
	if (cg_close(indexFile))
		cg_error_exit();

	/* save CGNS solution into mesh */
	elem_t *aElem = firstElem;
	while (aElem) {
		aElem->pVar[RHO] = rhoArr[aElem->id];
		aElem->pVar[VX]  = vxArr[aElem->id];
		aElem->pVar[VY]  = vyArr[aElem->id];
		aElem->pVar[P]   = pArr[aElem->id];

		aElem = aElem->next;
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
		cgnsReadSolution();
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
