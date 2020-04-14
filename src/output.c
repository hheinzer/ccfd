/*
 * output.c
 *
 * Created: Mon 23 Mar 2020 10:42:06 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "main.h"
#include "output.h"
#include "readInTools.h"
#include "timeDiscretization.h"
#include "analyze.h"
#include "equation.h"
#include "exactFunction.h"
#include "cgnslib.h"
#include "memTools.h"

/* extern variables */
char strOutFile[STRLEN];
double IOtimeInterval;
int IOiterInterval;
int iVisuProg;
char parameterFile[STRLEN];
FILE *resFile;
bool doErrorOutput;
outputTime_t *outputTimes;

/*
 * Initialize output
 */
void initOutput(void)
{
	printf("\nInitializing IO:\n");

	char *tmp = getStr("fileName", NULL);
	strcpy(strOutFile, tmp);
	free(tmp);

	IOtimeInterval = getDbl("IOtimeInterval", NULL);
	IOiterInterval = getDbl("IOiterInterval", NULL);
	iVisuProg = getInt("outputFormat", "1");
}

/*
 * tabular CSV output, only for 1D data
 */
void csvOutput(char fileName[STRLEN], double time, long iter, bool doExact)
{
	/* prepare data (only for equidistant grids) */
	double **flowData = dyn2DdblArray(nElems, NVAR);
	if (doExact) {
		long iElem = 0;
		elem_t *aElem = firstElem;
		double pVar[NVAR];
		while (aElem) {
			exactFunc(intExactFunc, aElem->bary, time, pVar);
			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = pVar[RHO];
			flowData[iElem][2] = pVar[VX];
			flowData[iElem][3] = pVar[P];
			iElem++;
			aElem = aElem->next;
		}
	} else {
		long iElem = 0;
		elem_t *aElem = firstElem;
		while (aElem) {
			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = aElem->pVar[RHO];
			flowData[iElem][2] = aElem->pVar[VX];
			flowData[iElem][3] = aElem->pVar[P];
			iElem++;
			aElem = aElem->next;
		}
	}

	/* sort array */
	long n = nElems;
	bool isSwapped = true;
	double tmp[NVAR];
	while (isSwapped && (n > 0)) {
		isSwapped = false;
		for (long iElem = 0; iElem < n - 1; ++iElem) {
			if (flowData[iElem][0] > flowData[iElem + 1][0]) {
				memcpy(tmp, flowData[iElem + 1], NVAR * sizeof(double));
				memcpy(flowData[iElem + 1], flowData[iElem], NVAR * sizeof(double));
				memcpy(flowData[iElem], tmp, NVAR * sizeof(double));
				isSwapped = true;
			}
		}
		n--;
	}

	/* write data */
	FILE *csvFile = fopen(fileName, "w");
	fprintf(csvFile, "CoordinateX, Density, Velocity, Pressure\n");
	for (long iElem = 0; iElem < nElems; ++iElem) {
		fprintf(csvFile, "%15.9f,%15.9f,%15.9f,%15.9f\n", flowData[iElem][0],
			flowData[iElem][1], flowData[iElem][2], flowData[iElem][3]);
	}
	fclose(csvFile);

	free(flowData);
}

/*
 * write solution to CGNS file
 */
void cgnsOutput(char fileName[STRLEN], double time, long iter, bool doExact)
{
	/* open solution file */
	int indexFile, indexBase, indexZone, indexSolution, indexField;
	if (cg_open(fileName, CG_MODE_WRITE, &indexFile))
		cg_error_exit();

	/* set up data for CGNS */
	cgsize_t iSize[3] = {nNodes, nElems, 0};

	/* create base */
	if (cg_base_write(indexFile, "Base", 2, 3, &indexBase))
		cg_error_exit();

	/* create zone */
	if (cg_zone_write(indexFile, indexBase, "Zone", iSize, Unstructured,
				&indexZone))
		cg_error_exit();

	/* link vertices and connectivity from the CGNS grid file */
	if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end"))
		cg_error_exit();

	if (cg_link_write("GridCoordinates", gridFile, "/Base/Zone/GridCoordinates"))
		cg_error_exit();

	if (nTrias > 0) {
		if (cg_link_write("Triangles", gridFile, "/Base/Zone/Triangles"))
			cg_error_exit();
	}

	if (nQuads > 0) {
		if (cg_link_write("Quadrilaterals", gridFile, "/Base/Zone/Quadrilaterals"))
			cg_error_exit();
	}

	/* create new solution zone */
	if (cg_sol_write(indexFile, indexBase, indexZone, "FlowSolution",
				CellCenter, &indexSolution))
		cg_error_exit();

	/* prepare density array, x-velocity array, y-velocity array, and
	 * pressure array */
	double *rhoArr = malloc(nElems * sizeof(double));
	double *vxArr  = malloc(nElems * sizeof(double));
	double *vyArr  = malloc(nElems * sizeof(double));
	double *vzArr  = malloc(nElems * sizeof(double));
	double *pArr   = malloc(nElems * sizeof(double));


	/* save solution in a CGNS compatible format */
	if (hasExactSolution) {
		elem_t *aElem = firstElem;

		while (aElem) {
			double pVar[NVAR];
			exactFunc(intExactFunc, aElem->bary, time, pVar);

			rhoArr[aElem->id] = pVar[RHO];
			vxArr[aElem->id]  = pVar[VX];
			vyArr[aElem->id]  = pVar[VY];
			vzArr[aElem->id]  = 0.0;
			pArr[aElem->id]   = pVar[P];

			aElem = aElem->next;
		}
	} else {
		elem_t *aElem = firstElem;

		while (aElem) {

			rhoArr[aElem->id] = aElem->pVar[RHO];
			vxArr[aElem->id]  = aElem->pVar[VX];
			vyArr[aElem->id]  = aElem->pVar[VY];
			vzArr[aElem->id]  = 0.0;
			pArr[aElem->id]   = aElem->pVar[P];

			aElem = aElem->next;
		}
	}

	/* write solution to CGNS file */
	if (cg_field_write(indexFile, indexBase, indexZone, indexSolution,
				RealDouble, "Density", rhoArr, &indexField))
		cg_error_exit();

	if (cg_field_write(indexFile, indexBase, indexZone, indexSolution,
				RealDouble, "VelocityX", vxArr, &indexField))
		cg_error_exit();

	if (cg_field_write(indexFile, indexBase, indexZone, indexSolution,
				RealDouble, "VelocityY", vyArr, &indexField))
		cg_error_exit();

	if (cg_field_write(indexFile, indexBase, indexZone, indexSolution,
				RealDouble, "VelocityZ", vzArr, &indexField))
		cg_error_exit();

	if (cg_field_write(indexFile, indexBase, indexZone, indexSolution,
				RealDouble, "Pressure", pArr, &indexField))
		cg_error_exit();

	/* write convergence information */
	if (cg_goto(indexFile, indexBase, "end"))
		cg_error_exit();

	char text[STRLEN];
	sprintf(text, "%20.12f %20.12f", time, timeOverall);

	if (cg_descriptor_write("ConvergenceInfo", text))
		cg_error_exit();

	/* close file */
	if (cg_close(indexFile))
		cg_error_exit();

	/* deallocate memory */
	free(rhoArr);
	free(vxArr);
	free(vyArr);
	free(vzArr);
	free(pArr);
}

/*
 * curved data output, only for 1D data
 */
void curveOutput(char fileName[STRLEN], double time, long iter, bool doExact)
{
	/* prepare data */
	double **flowData = dyn2DdblArray(nElems, NVAR);
	if (doExact) {
		long iElem = 0;
		elem_t *aElem = firstElem;
		while (aElem) {
			double pVar[NVAR];
			exactFunc(intExactFunc, aElem->bary, time, pVar);

			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = pVar[RHO];
			flowData[iElem][2] = pVar[VX];
			flowData[iElem][3] = pVar[P];

			iElem++;
			aElem = aElem->next;
		}
	} else {
		long iElem = 0;
		elem_t *aElem = firstElem;
		while (aElem) {
			flowData[iElem][0] = aElem->bary[X];
			flowData[iElem][1] = aElem->pVar[RHO];
			flowData[iElem][2] = aElem->pVar[VX];
			flowData[iElem][3] = aElem->pVar[P];

			iElem++;
			aElem = aElem->next;
		}
	}

	/* sort array */
	long n = nElems;
	bool isSwapped = true;
	double tmp[NVAR];
	while (isSwapped && (n > 0)) {
		isSwapped = false;
		for (long iElem = 0; iElem < n - 1; ++iElem) {
			if (flowData[iElem][0] > flowData[iElem + 1][0]) {
				memcpy(tmp, flowData[iElem + 1], NVAR * sizeof(double));
				memcpy(flowData[iElem + 1], flowData[iElem], NVAR * sizeof(double));
				memcpy(flowData[iElem], tmp, NVAR * sizeof(double));
				isSwapped = true;
			}
		}
		n--;
	}

	/* write data */
	FILE *curveFile = fopen(fileName, "w");

	fprintf(curveFile, "#Density\n");
	for (long iElem = 0; iElem < nElems; ++iElem) {
		fprintf(curveFile, "%15.9f,%15.9f\n", flowData[iElem][0],
				flowData[iElem][1]);
	}

	fprintf(curveFile, "\n#Velocity\n");
	for (long iElem = 0; iElem < nElems; ++iElem) {
		fprintf(curveFile, "%15.9f,%15.9f\n", flowData[iElem][0],
				flowData[iElem][2]);
	}

	fprintf(curveFile, "\n#Pressure\n");
	for (long iElem = 0; iElem < nElems; ++iElem) {
		fprintf(curveFile, "%15.9f,%15.9f\n", flowData[iElem][0],
				flowData[iElem][3]);
	}

	fclose(curveFile);

	free(flowData);
}

/*
 * data output in different formats
 */
void dataOutput(double time, long iter)
{
	/* output times */
	outputTime_t *outputTime = calloc(1, sizeof(outputTime_t));
	if (!outputTime) {
		printf("| ERROR: could not allocate outputTime\n");
		exit(1);
	}

	outputTime->time = time;
	outputTime->iter = iter;
	outputTime->next = outputTimes;
	outputTimes = outputTime;

	/* write flow solution */
	char fileName[2 * STRLEN];
	if (isStationary) {
		sprintf(fileName, "%s_%09ld", strOutFile, iter);
	} else {
		sprintf(fileName, "%s_%015.7f", strOutFile, time);
	}
	switch (iVisuProg) {
	case CGNS:
		strcat(fileName, ".cgns");
		cgnsOutput(fileName, time, iter, false);
		break;
	case CURVE:
		strcat(fileName, ".curve");
		curveOutput(fileName, time, iter, false);
		break;
	case DAT:
		strcat(fileName, ".csv");
		csvOutput(fileName, time, iter, false);
		break;
	default:
		printf("| ERROR: Output Format unknown\n");
		exit(1);
	}

	/* write exact solution, if applicable */
	if (hasExactSolution) {
		if (isStationary) {
			sprintf(fileName, "%s_ex_%09ld", strOutFile, iter);
		} else {
			sprintf(fileName, "%s_ex_%015.7f", strOutFile, time);
		}
		switch (iVisuProg) {
		case CGNS:
			strcat(fileName, ".cgns");
			cgnsOutput(fileName, time, iter, true);
			break;
		case CURVE:
			strcat(fileName, ".curve");
			curveOutput(fileName, time, iter, true);
			break;
		case DAT:
			strcat(fileName, ".csv");
			csvOutput(fileName, time, iter, true);
			break;
		default:
			printf("| ERROR: Output Format unknown\n");
			exit(1);
		}
	}
}

/*
 * finalize CGNS output
 */
void cgnsFinalizeOutput(void)
{
	/* count number of data outputs */
	cgsize_t nOutputs = 0;
	outputTime_t *outputTime = outputTimes;
	while (outputTime) {
		nOutputs++;
		outputTime = outputTime->next;
	}

	/* allocate and fill arrays */
	double times[nOutputs];
	long iters[nOutputs];
	char solutionNames[nOutputs][32];

	memset(times, 0, nOutputs * sizeof(double));
	memset(iters, 0, nOutputs * sizeof(long));
	memset(solutionNames, 0, nOutputs * 32 * sizeof(char));

	long iOutput = nOutputs;
	outputTime = outputTimes;
	while (outputTime) {
		iOutput--;
		times[iOutput] = outputTime->time;
		iters[iOutput] = outputTime->iter;
		sprintf(solutionNames[iOutput], "FlowSolution%09ld", iters[iOutput]);
		outputTime = outputTime->next;
	}

	int indexFile, indexBase, indexZone;
	/* open CGNS file */
	char masterFileName[STRLEN];
	strcat(strcpy(masterFileName, strOutFile), "_Master.cgns");
	if (cg_open(masterFileName, CG_MODE_WRITE, &indexFile))
		cg_error_exit();

	/* set up data for CGNS */
	cgsize_t iSize[3] = {nNodes, nElems, 0};

	/* create base */
	if (cg_base_write(indexFile, "Base", 2, 3, &indexBase))
		cg_error_exit();

	/* create zone */
	if (cg_zone_write(indexFile, indexBase, "Zone", iSize, Unstructured,
				&indexZone))
		cg_error_exit();

	/* link vertices and connectivity from the CGNS grid file */
	if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end"))
		cg_error_exit();

	if (cg_link_write("GridCoordinates", gridFile, "/Base/Zone/GridCoordinates"))
		cg_error_exit();

	if (nTrias > 0) {
		if (cg_link_write("Triangles", gridFile, "/Base/Zone/Triangles"))
			cg_error_exit();
	}

	if (nQuads > 0) {
		if (cg_link_write("Quadrilaterals", gridFile, "/Base/Zone/Quadrilaterals"))
			cg_error_exit();
	}

	/* link solutions */
	if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "end"))
		cg_error_exit();
	for (iOutput = 0; iOutput < nOutputs; iOutput++) {
		char solutionFileName[2 * STRLEN];
		if (isStationary) {
			sprintf(solutionFileName, "%s_%09ld.cgns", strOutFile,
					iters[iOutput]);
			if (cg_link_write(solutionNames[iOutput], solutionFileName,
					"/Base/Zone/FlowSolution"))
				cg_error_exit();
		} else {
			sprintf(solutionFileName, "%s_%015.7f.cgns", strOutFile,
					times[iOutput]);
			if (cg_link_write(solutionNames[iOutput], solutionFileName,
					"/Base/Zone/FlowSolution"))
				cg_error_exit();
		}
	}

	/* crete BaseIter node */
	if (cg_biter_write(indexFile, indexBase, "TimeIterValues", nOutputs))
		cg_error_exit();

	if (cg_goto(indexFile, indexBase, "BaseIterativeData_t", 1, "end"))
		cg_error_exit();

	cgsize_t tmp1[1] = {nOutputs};
	if (cg_array_write("TimeValues", RealDouble, 1, tmp1, times))
		cg_error_exit();

	if (cg_array_write("IterationValues", LongInteger, 1, tmp1, iters))
		cg_error_exit();

	/* create ZoneIter node */
	if (cg_ziter_write(indexFile, indexBase, indexZone, "ZoneIterativeData"))
		cg_error_exit();

	if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "ZoneIterativeData_t",
				1, "end"))
		cg_error_exit();

	cgsize_t tmp2[2] = {32, nOutputs};
	if (cg_array_write("FlowSolutionPointers", Character, 2, tmp2, solutionNames))
		cg_error_exit();

	/* close file */
	if (cg_close(indexFile))
		cg_error_exit();
}

/*
 * finalize data output
 */
void finalizeDataOutput(void)
{
	if (iVisuProg == CGNS) {
		cgnsFinalizeOutput();
	}
}

/*
 * write CGNS mesh
 */
void cgnsWriteMesh(void)
{
	/* set up data */
	cgsize_t iSize[3] = {nNodes, nElems, 0};

	/* allocate element arrays */
	cgsize_t **trias = dyn2DcgsizeArray(nTrias, 3);
	cgsize_t **quads = dyn2DcgsizeArray(nQuads, 4);
	cgsize_t **BCsides = dyn2DcgsizeArray(nBCsides, 2);
	double **nodes = dyn2DdblArray(NDIM + 1, nNodes);

	/* save verticies in a CGNS compatible format */
	long iNode = 0;
	node_t *aNode = firstNode;
	while (aNode) {
		nodes[X][iNode] = aNode->x[X];
		nodes[Y][iNode] = aNode->x[Y];
		nodes[2][iNode] = 0.0;

		aNode->id = iNode++;
		aNode = aNode->next;
	}

	/* save element connectivity in a CGNS compatible format */
	long iTria = 0, iQuad = 0;
	elem_t *aElem = firstElem;
	while (aElem) {
		switch (aElem->elemType) {
		case 3:
			trias[iTria  ][0] = aElem->node[0]->id + 1;
			trias[iTria  ][1] = aElem->node[1]->id + 1;
			trias[iTria++][2] = aElem->node[2]->id + 1;
			break;
		case 4:
			quads[iQuad  ][0] = aElem->node[0]->id + 1;
			quads[iQuad  ][1] = aElem->node[1]->id + 1;
			quads[iQuad  ][2] = aElem->node[2]->id + 1;
			quads[iQuad++][3] = aElem->node[3]->id + 1;
			break;
		}

		aElem = aElem->next;
	}

	/* save boundary elements */
	cgsize_t BCpartition[nBC + 1];
	if (!isPeriodic) {
		BCpartition[0] = (cgsize_t)nElems;
		cgsize_t nSide = nElems;
		long iBCside = 0;
		int iBC = 1;
		boundary_t *aBC = firstBC;
		while (aBC) {
			sidePtr_t *aBCside = firstBCside;
			while (aBCside) {
				if (aBCside->side->BC == aBC) {
					BCsides[iBCside  ][0] =
						aBCside->side->node[0]->id + 1;
					BCsides[iBCside++][1] =
						aBCside->side->node[1]->id + 1;
					nSide++;
				}

				aBCside = aBCside->next;
			}

			BCpartition[iBC++] = nSide;
			aBC = aBC->next;
		}
	}

	int indexFile, indexBase, indexZone, indexGrid, indexCoordinate,
	    indexSection, indexBounary;
	/* open CGNS grid file for writing */
	if (cg_open(gridFile, CG_MODE_WRITE, &indexFile))
		cg_error_exit();

	/* write coordinate base to CGNS file */
	if (cg_base_write(indexFile, "Base", 2, 3, &indexBase))
		cg_error_exit();

	/* write computational zone and grid to CGNS file */
	if (cg_zone_write(indexFile, indexBase, "Zone", iSize, Unstructured, &indexZone))
		cg_error_exit();
	if (cg_grid_write(indexFile, indexBase, indexZone, "GridCoordinates", &indexGrid))
		cg_error_exit();

	/* write cordinates to file */
	if (cg_coord_write(indexFile, indexBase, indexZone, RealDouble, "CoordinateX",
			nodes[X], &indexCoordinate))
		cg_error_exit();
	if (cg_coord_write(indexFile, indexBase, indexZone, RealDouble, "CoordinateY",
			nodes[Y], &indexCoordinate))
		cg_error_exit();
	if (cg_coord_write(indexFile, indexBase, indexZone, RealDouble, "CoordinateZ",
			nodes[2], &indexCoordinate))
		cg_error_exit();

	/* write the element connectivity to the CGNS file */
	if (nTrias > 0) {
		if (cg_section_write(indexFile, indexBase, indexZone, "Triangles", TRI_3,
				1, nTrias, 0, trias[0], &indexSection))
			cg_error_exit();
	}
	if (nQuads > 0) {
		if (cg_section_write(indexFile, indexBase, indexZone, "Quadrilaterals", QUAD_4,
				nTrias + 1, nTrias + nQuads, 0, quads[0], &indexSection))
			cg_error_exit();
	}

	/* write the boundary connectivity to the CGNS file */
	if (!isPeriodic) {
		if (cg_section_write(indexFile, indexBase, indexZone, "Boundaries",
					BAR_2, nElems + 1, nElems + nBCsides, 0,
					BCsides[0], &indexSection))
			cg_error_exit();

		/* write the boundary connectivity to the CGNS file */
		long iBC = 0;
		boundary_t *aBC = firstBC;
		while (aBC) {
			char cBC[STRLEN];
			sprintf(cBC, "%d", aBC->BCtype * 100 + aBC->BCid);
			cgsize_t iPoints[nBCsides];
			long nCount = 0;
			for (cgsize_t n = BCpartition[iBC] + 1; n <= BCpartition[iBC + 1]; ++n) {
				iPoints[nCount++] = n;
			}

			if (cg_boco_write(indexFile, indexBase, indexZone, cBC, BCGeneral,
					ElementList, nCount, iPoints, &indexBounary))
				cg_error_exit();

			aBC = aBC->next;
			iBC++;
		}
	}

	/* close CGNS file */
	if (cg_close(indexFile))
		cg_error_exit();

	free(nodes);
	free(trias);
	free(quads);
	free(BCsides);
}

/*
 * free all memory that was allocated for the output times
 */
void freeOutputTimes(void)
{
	outputTime_t *outputTime = outputTimes;
	while (outputTime) {
		if (outputTime->next) {
			outputTime_t *tmp = outputTime;
			outputTime = outputTime->next;
			free(tmp);
		} else {
			free(outputTime);
			break;
		}
	}
}
