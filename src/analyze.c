/*
 * analyze.c
 *
 * Created: Sun 29 Mar 2020 06:29:36 PM CEST
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "analyze.h"
#include "readInTools.h"
#include "timeDiscretization.h"
#include "output.h"
#include "memTools.h"

/* extern variables */
bool doCalcWing;
wing_t wing;
recordPoint_t recordPoint;
bool hasExactSolution;

/*
 * initialize recording points
 */
void initRecordPoints(void)
{
	for (long iPt = 0; iPt < recordPoint.nPoints; ++iPt) {
		elem_t *aElem = firstElem;
		bool isInside;
		while (aElem) {
			isInside = true;
			for (int i = 0; i < aElem->elemType; ++i) {
				int iNode1 = i;
				int iNode2 = i + 1;
				if (iNode2 >= aElem->elemType) {
					iNode2 = 0;
				}
				double u[NDIM];
				u[X] = aElem->node[iNode2]->x[Y] - aElem->node[iNode1]->x[Y];
				u[Y] = aElem->node[iNode1]->x[X] - aElem->node[iNode2]->x[X];
				double tmp = sqrt(u[X] * u[X] + u[X] * u[X]);
				u[X] /= tmp;
				u[Y] /= tmp;

				double x[NDIM];
				x[X] = recordPoint.x[iPt][X] - aElem->node[iNode1]->x[X];
				x[Y] = recordPoint.x[iPt][Y] - aElem->node[iNode1]->x[Y];

				double projection = u[X] * x[X] + u[Y] * x[Y];
				if (projection > 0.0) {
					isInside = false;
				}
			}

			if (isInside) {
				recordPoint.elem[iPt] = aElem;
				break;
			}

			aElem = aElem->next;
		}

		if (!isInside) {
			printf("| ERROR: Record Point # %ld is not in Domain\n", iPt);
			exit(1);
		}

		/* open file for writing */
		char ioFileName[2 * STRLEN];
		sprintf(ioFileName, "%s_recordPoint_%ld.csv", strOutFile, iPt);
		recordPoint.ioFile[iPt] = fopen(ioFileName, "w");
		fprintf(recordPoint.ioFile[iPt], "Time, Density, VelocityX, VelocityY, Pressure\n");
	}
}

/*
 * initialize required data for calculation of cl and cd
 */
void initWing(void)
{
	/* get a pointer to the wing section's BC */
	boundary_t *aBC = firstBC, *wingBC;
	while (aBC) {
		if (aBC->BCtype * 100 + aBC->BCid == wing.wallId) {
			wingBC = aBC;
		}
		aBC = aBC->next;
	}

	/* build lists for pressure and suction side of the wing section */
	sidePtr_t *aBCsidePtr = firstBCside;
	while (aBCsidePtr) {
		if (wingBC == aBCsidePtr->side->BC) {
			bool isInserted = false;
			sidePtr_t *aInsSidePtr = malloc(sizeof(sidePtr_t));
			aInsSidePtr->side = aBCsidePtr->side->connection;
			double xIns[NDIM];
			xIns[X] = aInsSidePtr->side->GP[X] + aInsSidePtr->side->elem->bary[X];
			xIns[Y] = aInsSidePtr->side->GP[Y] + aInsSidePtr->side->elem->bary[Y];
			aInsSidePtr->next = NULL;

			/* check the direction of the side's normal vector */
			sidePtr_t *aSidePtr;
			if (aBCsidePtr->side->n[Y] <= 0.0) {
				if (!wing.firstSuctionSide) {
					wing.firstSuctionSide = aInsSidePtr;
					isInserted = true;
				}

				aSidePtr = wing.firstSuctionSide;
				double x[NDIM];
				x[X] = aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X];
				x[Y] = aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y];
				if (xIns[X] < x[X]) {
					aInsSidePtr->next = aSidePtr;
					wing.firstSuctionSide = aInsSidePtr;
					isInserted = true;
				}
				aSidePtr = wing.firstSuctionSide;
			} else {
				if (!wing.firstPressureSide) {
					wing.firstPressureSide = aInsSidePtr;
					isInserted = true;
				}

				aSidePtr = wing.firstPressureSide;
				double x[NDIM];
				x[X] = aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X];
				x[Y] = aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y];
				if (xIns[X] < x[X]) {
					aInsSidePtr->next = aSidePtr;
					wing.firstPressureSide = aInsSidePtr;
					isInserted = true;
				}
				aSidePtr = wing.firstPressureSide;
			}

			sidePtr_t *aLastSidePtr = aSidePtr;
			aSidePtr = aSidePtr->next;

			/* sort current side according to its coordinates into
			 * the list */
			while (aSidePtr && !isInserted) {
				double x[NDIM];
				x[X] = aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X];
				x[Y] = aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y];
				if (xIns[X] < x[X]) {
					aLastSidePtr->next = aInsSidePtr;
					aInsSidePtr->next = aSidePtr;
					isInserted = true;
				}
				aLastSidePtr = aSidePtr;
				aSidePtr = aSidePtr->next;
			}

			if (!isInserted) {
				aLastSidePtr->next = aInsSidePtr;
			}
		}

		aBCsidePtr = aBCsidePtr->next;
	}
}

/*
 * get wing parameters
 */
void readWing(void)
{
	printf("\nInitializing Wing:\n");
	wing.refLength = getDbl("referenceLength", NULL);
	wing.wallId = getInt("wall_id", NULL);
}

/*
 * initialize analyze
 */
void initAnalyze(void)
{
	printf("\nInitializing Analysis:\n");
	hasExactSolution = getBool("exactSolution", "F");
	doCalcWing = getBool("calcWing", "F");

	char resFileName[STRLEN];
	if (doCalcWing || isStationary) {
		strcat(strcpy(resFileName, strOutFile), "_analysis.csv");

		if (isRestart) {
			resFile = fopen(resFileName, "r");
			if (resFile) {
				/* file exists, close and append to it */
				fclose(resFile);
				resFile = fopen(resFileName, "a");
			} else {
				/* file does not exist, create it and write
				 * header */
				resFile = fopen(resFileName, "w");
				if (doCalcWing) {
					fprintf(resFile, "Iter, Time, Residual(%s), C_L, C_D\n",
							abortVariableName);
				} else {
					fprintf(resFile, "Iter, Time, Residual(RHO), Residual(VX), Residual(VY), Residual(E)\n");
				}
			}
		} else {
			resFile = fopen(resFileName, "w");
			if (doCalcWing) {
				fprintf(resFile, "Iter, Time, Residual(%s), C_L, C_D\n",
						abortVariableName);
			} else {
				fprintf(resFile, "Iter, Time, Residual(RHO), Residual(VX), Residual(VY), Residual(E)\n");
			}
		}
	}

	/* gnuplot file for residuals */
	char demFileName[STRLEN];
	if (doCalcWing) {
		readWing();

		/* residuals plot file */
		strcat(strcpy(demFileName, strOutFile), "_residuals.dem");
		FILE *demFile = fopen(demFileName, "w");
		fprintf(demFile, "set title 'Residual Plot'\n");
		fprintf(demFile, "set logscale y\n");
		fprintf(demFile, "plot '%s' using 1:3 title '%s' with lines\n",
				resFileName, abortVariableName);
		fprintf(demFile, "pause -1");
		fclose(demFile);

		/* clcd plot file */
		strcat(strcpy(demFileName, strOutFile), "_clcd.dem");
		demFile = fopen(demFileName, "w");
		fprintf(demFile, "set title 'cl/cd Plot'\n");
		fprintf(demFile, "plot '%s' using 1:4 title 'cl' with lines, \\\n",
				resFileName);
		fprintf(demFile, "     '%s' using 1:5 title 'cd' with lines\n",
				resFileName);
		fprintf(demFile, "pause -1");
		fclose(demFile);
	} else if (isStationary) {
		strcat(strcpy(demFileName, strOutFile), "_residuals.dem");
		FILE *demFile = fopen(demFileName, "w");
		fprintf(demFile, "set title 'Residual Plot'\n");
		fprintf(demFile, "set logscale y\n");
		fprintf(demFile, "plot '%s' using 1:3 title 'RHO' with lines, \\\n",
				resFileName);
		fprintf(demFile, "plot '%s' using 1:4 title 'MX' with lines, \\\n",
				resFileName);
		fprintf(demFile, "plot '%s' using 1:5 title 'MY' with lines, \\\n",
				resFileName);
		fprintf(demFile, "plot '%s' using 1:6 title 'E' with lines\n",
				resFileName);
		fprintf(demFile, "pause -1");
		fclose(demFile);
	}

	if (doCalcWing) {
		initWing();
	}

	if (recordPoint.nPoints > 0) {
		initRecordPoints();
	}
}

/*
 * compute aerodynamic coefficients and extract values at record points
 */
void analyze(double time, long iter, double resIter[NVAR + 2])
{

}

/*
 * calculate L1, L2, and Linf error norms
 */
void calcErrors(double time)
{

}
