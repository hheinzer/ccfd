/** \file
 *
 * \brief Contains functions for analyzing flow results
 *
 * \author hhh
 * \date Sun 29 Mar 2020 06:29:36 PM CEST
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
#include "reconstruction.h"
#include "exactFunction.h"
#include "equation.h"
#include "initialCondition.h"

/* extern variables */
bool doCalcWing;			/**< calculate CL CD flag */
wing_t wing;				/**< holds data for coefficient calculation */
recordPoint_t recordPoint;		/**< record flow field at a specific point */
bool hasExactSolution;			/**< exact solution existence flag */

/**
 * \brief Initialize recording points
 */
void initRecordPoints(void)
{
	recordPoint.x = dyn2DdblArray(recordPoint.nPoints, 2);

	for (long iPt = 0; iPt < recordPoint.nPoints; ++iPt) {
		double *coords = getDblArray("pointCoords", NDIM, NULL);
		recordPoint.x[iPt][X] = coords[X];
		recordPoint.x[iPt][Y] = coords[Y];
		free(coords);

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

/**
 * \brief Initialize required data for calculation of CL and CD
 */
void initWing(void)
{
	/* get a pointer to the wing section's BC */
	boundary_t *aBC = firstBC, *wingBC = NULL;
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
			if (!aInsSidePtr) {
				printf("| ERROR: could not allocate aInsSidePtr\n");
				exit(1);
			}

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

/**
 * \brief Get parameters for calculating aerodynamic coefficients
 */
void readWing(void)
{
	printf("\nInitializing Wing:\n");
	wing.refLength = getDbl("referenceLength", "1.0");
	wing.wallId = getInt("wall_id", NULL);
}

/**
 * \brief Initialize the analysis function
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
					fprintf(resFile, "Iter, Time, Residual(%s), CL, CD\n",
							abortVariableName);
				} else {
					fprintf(resFile, "Iter, Time, Residual(RHO), Residual(VX), Residual(VY), Residual(E)\n");
				}
			}
		} else {
			resFile = fopen(resFileName, "w");
			if (doCalcWing) {
				fprintf(resFile, "Iter, Time, Residual(%s), CL, CD\n",
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
		strcat(strcpy(demFileName, strOutFile), "_CLCD.dem");
		demFile = fopen(demFileName, "w");
		fprintf(demFile, "set title 'CL/CD Plot'\n");
		fprintf(demFile, "plot '%s' using 1:4 title 'CL' with lines, \\\n",
				resFileName);
		fprintf(demFile, "     '%s' using 1:5 title 'CD' with lines\n",
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
		fprintf(demFile, "     '%s' using 1:4 title 'MX' with lines, \\\n",
				resFileName);
		fprintf(demFile, "     '%s' using 1:5 title 'MY' with lines, \\\n",
				resFileName);
		fprintf(demFile, "     '%s' using 1:6 title 'E' with lines\n",
				resFileName);
		fprintf(demFile, "pause -1");
		fclose(demFile);
	}

	if (doCalcWing) {
		initWing();
	}

	recordPoint.nPoints = getInt("nRecordPoints", "0");
	if (recordPoint.nPoints > 0) {
		initRecordPoints();
	}
}

/**
 * \brief Calculate CL and CD around the specified wall
 */
void calcCoef(void)
{
	/* initialize values */
	double cl = 0.0, cd = 0.0;
	double v = refState[0][VX] / cos(alpha * pi / 180.0);
	double qInfQ = 1.0 / (refState[0][RHO] * 0.5 * v * v);
	double qInfLq = qInfQ / wing.refLength;
	double pInf = refState[0][P];

	/* open file for writing cp data */
	char pressureFileName[STRLEN];
	strcat(strcpy(pressureFileName, strOutFile), "_CP_pressureSide.csv");
	FILE *cp1File = fopen(pressureFileName, "w");
	if (!cp1File) {
		printf("| ERROR: Cannot open Output File for CP I/O\n");
		exit(1);
	}

	/* presure side: CL, CD, CP */
	fprintf(cp1File, "x, y, phi, CP_pressureSide\n");
	sidePtr_t *aSidePtr = wing.firstPressureSide;
	while (aSidePtr) {
		double p0 = aSidePtr->side->pVar[P];
		double n[NDIM];
		n[X]  = aSidePtr->side->n[X];
		n[Y]  = aSidePtr->side->n[Y];
		double len = aSidePtr->side->len;
		cl += n[Y] * p0 * len;
		cd += n[X] * p0 * len;
		double cp = (p0 - pInf) * qInfQ;
		fprintf(cp1File, "%15.7f,%15.7f,%15.7f,%15.7f\n",
			aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X],
			aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y],
			atan2(aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y],
			      aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X]),
			cp);
		aSidePtr = aSidePtr->next;
	}
	fclose(cp1File);

	/* open file for writing CP data */
	char suctionFileName[STRLEN];
	strcat(strcpy(suctionFileName, strOutFile), "_CP_suctionSide.csv");
	FILE *cp2File = fopen(suctionFileName, "w");
	if (!cp2File) {
		printf("| ERROR: Cannot open Output File for CP I/O\n");
		exit(1);
	}

	fprintf(cp2File, "x, y, phi, CP_suctionSide\n");
	aSidePtr = wing.firstSuctionSide;
	while (aSidePtr) {
		double p0 = aSidePtr->side->pVar[P];
		double n[NDIM];
		n[X]  = aSidePtr->side->n[X];
		n[Y]  = aSidePtr->side->n[Y];
		double len = aSidePtr->side->len;
		cl += n[Y] * p0 * len;
		cd += n[X] * p0 * len;
		double cp = (p0 - pInf) * qInfQ;
		fprintf(cp2File, "%15.7f,%15.7f,%15.7f,%15.7f\n",
			aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X],
			aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y],
			atan2(aSidePtr->side->GP[Y] + aSidePtr->side->elem->bary[Y],
			      aSidePtr->side->GP[X] + aSidePtr->side->elem->bary[X]),
			cp);
		aSidePtr = aSidePtr->next;
	}
	fclose(cp2File);

	cl *= qInfLq;
	cd *= qInfLq;

	/* saving as global variables */
	double alphaLoc = alpha * pi / 180.0;
	wing.cl = cl * cos(alphaLoc) - cd * sin(alphaLoc);
	wing.cd = cd * cos(alphaLoc) + cl * sin(alphaLoc);

	/* write gnuplot file for CP plot */
	char demFileName[STRLEN];
	strcat(strcpy(demFileName, strOutFile), "_CP.dem");
	FILE *demFile = fopen(demFileName, "r");
	if (!demFile) {
		demFile = fopen(demFileName, "w");
		fprintf(demFile, "unset key\n");
		fprintf(demFile, "set xlabel 'x'\n");
		fprintf(demFile, "set ylabel 'cp'\n");
		fprintf(demFile, "set yrange [*:*] reverse\n");
		fprintf(demFile, "plot '%s' using 1:4 w l lc rgb 'black', \\\n",
				pressureFileName);
		fprintf(demFile, "     '%s' using 1:4 w l lc rgb 'black'\n",
				suctionFileName);
		fprintf(demFile, "pause -1");
	}
	fclose(demFile);
}

/**
 * \brief Evaluation of recording points
 * \param[in] time Calculation time at output
 */
void evalRecordPoints(double time)
{
	for (long iPt = 0; iPt < recordPoint.nPoints; ++iPt) {
		elem_t *aElem = recordPoint.elem[iPt];
		fprintf(recordPoint.ioFile[iPt],
			"%20.12f,%20.12f,%20.12f,%20.12f,%20.12f\n",
			time + aElem->dt, aElem->pVar[RHO], aElem->pVar[VX],
			aElem->pVar[VY], aElem->pVar[P]);
	}
}

/**
 * \brief Compute aerodynamic coefficients and extract values at record points
 * \param[in] time Calculation time at output
 * \param[in] iter Iteration count at output
 * \param[in,out] resIter The residual vector containing the CL and CD
 *	residuals at 4th and 5th index position
 */
void analyze(double time, long iter, double resIter[NVAR + 2])
{
	/* record points */
	if (recordPoint.nPoints > 0) {
		evalRecordPoints(time);
	}

	/* aerodynamic coefficients */
	if (doCalcWing) {
		resIter[4] = wing.cl;
		resIter[5] = wing.cd;

		calcCoef();

		resIter[4] = fabs(resIter[4] - wing.cl) / firstElem->dt;
		resIter[5] = fabs(resIter[5] - wing.cd) / firstElem->dt;

		fprintf(resFile, "%7ld, %13.8f, %15.8e, %15.10f, %15.10f\n",
			iter, time + firstElem->dt, resIter[abortVariable],
			wing.cl, wing.cd);
	} else {
		if (isStationary) {
			fprintf(resFile, "%7ld, %13.8f, %15.8e, %15.8e, %15.8e, %15.8e\n",
				iter, time + firstElem->dt, resIter[RHO],
				resIter[VX], resIter[VY], resIter[E]);
		}
	}
}

/**
 * \brief Calculate L1, L2, and Linf error norms, between the flow solution and
 *	the exact solution
 * \param[in] time Calculation time at output
 */
void calcErrors(double time)
{
	double L1[NVAR] = {0.0}, L2[NVAR] = {0.0}, Linf[NVAR] = {0.0};

	spatialReconstruction(time);

	#pragma omp parallel for reduction(max:Linf), reduction(+:L1,L2)
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		for (int iGP = 0; iGP < aElem->nGP; ++iGP) {
			/* exact function */
			double pVarEx[NVAR];
			exactFunc(intExactFunc, aElem->xGP[iGP], time, pVarEx);

			/* space expansion to obtain pVar at current GP */
			double dx[NDIM] = {
				aElem->xGP[iGP][X] - aElem->bary[X],
				aElem->xGP[iGP][Y] - aElem->bary[Y]};

			double pVar[NVAR] = {
				aElem->pVar[RHO] + dx[X] * aElem->u_x[RHO] + dx[Y] * aElem->u_y[RHO],
				aElem->pVar[VX]  + dx[X] * aElem->u_x[VX]  + dx[Y] * aElem->u_y[VX],
				aElem->pVar[VY]  + dx[X] * aElem->u_x[VY]  + dx[Y] * aElem->u_y[VY],
				aElem->pVar[P]   + dx[X] * aElem->u_x[P]   + dx[Y] * aElem->u_y[P]};

			/* compute errors at GP */
			double err[NVAR] = {
				fabs(pVarEx[RHO] - pVar[RHO]),
				fabs(pVarEx[VX]  - pVar[VX]),
				fabs(pVarEx[VY]  - pVar[VY]),
				fabs(pVarEx[P]   - pVar[P])};

			/* update Linf error norm */
			Linf[RHO] = fmax(Linf[RHO], err[RHO]);
			Linf[VX]  = fmax(Linf[VX],  err[VX]);
			Linf[VY]  = fmax(Linf[VY],  err[VY]);
			Linf[P]   = fmax(Linf[P],   err[P]);

			/* update L1 error norm */
			L1[RHO] += err[RHO] * aElem->wGP[iGP];
			L1[VX]  += err[VX]  * aElem->wGP[iGP];
			L1[VY]  += err[VY]  * aElem->wGP[iGP];
			L1[P]   += err[P]   * aElem->wGP[iGP];

			/* update L2 error norm */
			L2[RHO] += err[RHO] * err[RHO] * aElem->wGP[iGP];
			L2[VX]  += err[VX]  * err[VX]  * aElem->wGP[iGP];
			L2[VY]  += err[VY]  * err[VY]  * aElem->wGP[iGP];
			L2[P]   += err[P]   * err[P]   * aElem->wGP[iGP];
		}
	}

	/* finalize L1 and L2 */
	L1[RHO] *= totalArea_q;
	L1[VX]  *= totalArea_q;
	L1[VY]  *= totalArea_q;
	L1[P]   *= totalArea_q;

	L2[RHO] = sqrt(L2[RHO] * totalArea_q);
	L2[VX]  = sqrt(L2[VX]  * totalArea_q);
	L2[VY]  = sqrt(L2[VY]  * totalArea_q);
	L2[P]   = sqrt(L2[P]   * totalArea_q);

	/* input and output */
	printf("\nError Analysis at t = %g\n", time);
	printf("|                RHO           VX           VY            P\n");
	printf("| L1  :  %7.5e  %7.5e  %7.5e  %7.5e\n", L1[RHO], L1[VX], L1[VY], L1[P]);
	printf("| L2  :  %7.5e  %7.5e  %7.5e  %7.5e\n", L2[RHO], L2[VX], L2[VY], L2[P]);
	printf("| Linf:  %7.5e  %7.5e  %7.5e  %7.5e\n", Linf[RHO], Linf[VX], Linf[VY], Linf[P]);
}

/**
 * \brief Calculate the global residual of the conservative variables
 * \param[in,out] resIter The residual vector containing the residuals
 *	of the conservative variables at the first four positions
 */
void globalResidual(double resIter[NVAR + 2])
{
	memset(resIter, 0, (NVAR + 2) * sizeof(double));

	#pragma omp parallel for reduction(+:resIter[:NVAR + 2])
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		resIter[RHO] += aElem->area * aElem->u_t[RHO] * aElem->u_t[RHO];
		resIter[MX]  += aElem->area * aElem->u_t[MX]  * aElem->u_t[MX];
		resIter[MY]  += aElem->area * aElem->u_t[MY]  * aElem->u_t[MY];
		resIter[E]   += aElem->area * aElem->u_t[E]   * aElem->u_t[E];
	}

	/* compute 2-Norm of the residual */
	resIter[RHO] = sqrt(resIter[RHO] * totalArea_q);
	resIter[MX]  = sqrt(resIter[MX]  * totalArea_q);
	resIter[MY]  = sqrt(resIter[MY]  * totalArea_q);
	resIter[E]   = sqrt(resIter[E]   * totalArea_q);
}

/**
 * \brief Free all memory that was allocated for the wing evaluation
 */
void freeAnalyze(void)
{
	if (doCalcWing) {
		sidePtr_t *aSidePtr = wing.firstSuctionSide;
		while (aSidePtr) {
			if (aSidePtr->next) {
				sidePtr_t *tmp = aSidePtr;
				aSidePtr = aSidePtr->next;
				free(tmp);
			} else {
				free(aSidePtr);
				break;
			}
		}

		aSidePtr = wing.firstPressureSide;
		while (aSidePtr) {
			if (aSidePtr->next) {
				sidePtr_t *tmp = aSidePtr;
				aSidePtr = aSidePtr->next;
				free(tmp);
			} else {
				free(aSidePtr);
				break;
			}
		}
	}
}
