/*
 * mesh.c
 *
 * Created: Tue 24 Mar 2020 12:45:57 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "main.h"
#include "mesh.h"
#include "cgnslib.h"
#include "readInTools.h"
#include "output.h"
#include "timeDiscretization.h"

/* extern variables */
char parameterFile[STRLEN],
     strMeshFormat[STRLEN],
     strMeshFile[STRLEN],
     strIniCondFile[STRLEN];

cartMesh_t cartMesh;
char gridFile[STRLEN];

int meshType;
int meshFormat;

long nNodes;

long nElems;
long nTrias;
long nQuads;

long nSides;
long nBCsides;
long nInnerSides;

double totalArea_q;
double xMin, xMax;
double yMin, yMax;
double dxRef;

elem_t *elem;
side_t *side;
side_t *BCside;

elem_t *firstElem;
node_t *firstNode;
side_t *firstSide;
side_t *firstBCside;

typedef struct sideList_t sideList_t;
struct sideList_t {
	long node[NDIM];
	int BC;
	side_t *side;
	bool isRotated;
};

/*
 * allocate a dynamic 2D array of integers
 */
long **dyn2DintArray(long I, int J)
{
	long **arr = malloc(sizeof(long *) * I + sizeof(long) * I * J);
	long *ptr = (long *)(arr + I);
	for (int i = 0; i < I; ++i) {
		arr[i] = (ptr + J * i);
	}
	return arr;
}

/*
 * allocate a dynamic 2D array of doubles
 */
double **dyn2DdblArray(unsigned long int I, int J)
{
	double **arr = malloc(sizeof(double *) * I + sizeof(double) * I * J);
	double *ptr = (double *)(arr + I);
	for (int i = 0; i < I; ++i) {
		arr[i] = (ptr + J * i);
	}
	return arr;
}

/*
 * create a cartesian mesh
 */
void createCartMesh(
	double ***vertex, long *nVertices, long ***BCedge, long *nBCedges, long ***quad)
{
	/* set up necessary variables */
	long iMax = cartMesh.iMax;
	long jMax = cartMesh.jMax;

	*nVertices = (iMax + 1) * (jMax + 1);
	nQuads = iMax * jMax;
	*nBCedges = 2 * (iMax + jMax);

	double dx = (xMax - xMin) / iMax;
	double dy = (yMax - yMin) / jMax;

	*quad = dyn2DintArray(nQuads, 5);
	*BCedge = dyn2DintArray(*nBCedges, 3);
	*vertex = dyn2DdblArray(*nVertices, 2);

	/* build vertices */
	double y = yMin, x;
	long k = 0;
	for (long j = 0; j < jMax + 1; ++j) {
		x = xMin;
		for (long i = 0; i < iMax + 1; ++i) {
			(*vertex)[k][X] = x;
			(*vertex)[k][Y] = y;
			//printf("%ld %g %g\n", k, (*vertex)[k][X], (*vertex)[k][Y]);
			k++;
			x += dx;
		}
		y += dy;
	}

	/* build elements */
	k = 0;
	for (long j = 0; j < jMax; ++j) {
		for (long i = 0; i < iMax; ++i) {
			(*quad)[k][0] = k + j;
			(*quad)[k][1] = (*quad)[k][0] + 1;
			(*quad)[k][2] = (*quad)[k][1] + iMax + 1;
			(*quad)[k][3] = (*quad)[k][2] - 1;
			(*quad)[k][4] = 1;
			//printf("%ld %ld %ld %ld %ld %ld\n", k, (*quad)[k][0], (*quad)[k][1], (*quad)[k][2], (*quad)[k][3], (*quad)[k][4]);
			k++;
		}
	}

	/* build BCedges */
	k = 0;
	/* bottom side */
	printf("| Bottom Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[BOTTOM]);
	for (int iBC = 0; iBC < cartMesh.nBC[BOTTOM]; ++iBC) {
		printf("|   BC Type: %d\n", cartMesh.BCtype[BOTTOM][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[BOTTOM][iBC][0], cartMesh.BCrange[BOTTOM][iBC][1]);
		for (int i = cartMesh.BCrange[BOTTOM][iBC][0] - 1; i < cartMesh.BCrange[BOTTOM][iBC][1]; ++i) {
			(*BCedge)[k][0] = i;
			(*BCedge)[k][1] = i + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[BOTTOM][iBC];
			//printf("%ld %ld %ld %ld\n", k, (*BCedge)[k][0], (*BCedge)[k][1], (*BCedge)[k][2]);
			k++;
		}
	}
	/* right side */
	printf("| Right Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[RIGHT]);
	for (int iBC = 0; iBC < cartMesh.nBC[RIGHT]; ++iBC) {
		printf("|   BC Type: %d\n", cartMesh.BCtype[RIGHT][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[RIGHT][iBC][0], cartMesh.BCrange[RIGHT][iBC][1]);
		for (int j = cartMesh.BCrange[RIGHT][iBC][0] - 1; j < cartMesh.BCrange[RIGHT][iBC][1]; ++j) {
			(*BCedge)[k][0] = (iMax + 1) * (j + 1) - 1;
			(*BCedge)[k][1] = (iMax + 1) * (j + 2) - 1;
			(*BCedge)[k][2] = cartMesh.BCtype[RIGHT][iBC];
			//printf("%ld %ld %ld %ld\n", k, (*BCedge)[k][0], (*BCedge)[k][1], (*BCedge)[k][2]);
			k++;
		}
	}
	/* top side */
	printf("| Top Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[TOP]);
	for (int iBC = 0; iBC < cartMesh.nBC[TOP]; ++iBC) {
		printf("|   BC Type: %d\n", cartMesh.BCtype[TOP][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[TOP][iBC][0], cartMesh.BCrange[TOP][iBC][1]);
		for (int i = cartMesh.BCrange[TOP][iBC][0] - 1; i < cartMesh.BCrange[TOP][iBC][1]; ++i) {
			(*BCedge)[k][0] = (iMax + 1) * jMax + i;
			(*BCedge)[k][1] = (*BCedge)[k][0] + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[TOP][iBC];
			//printf("%ld %ld %ld %ld\n", k, (*BCedge)[k][0], (*BCedge)[k][1], (*BCedge)[k][2]);
			k++;
		}
	}
	/* left side */
	printf("| Left Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[LEFT]);
	for (int iBC = 0; iBC < cartMesh.nBC[LEFT]; ++iBC) {
		printf("|   BC Type: %d\n", cartMesh.BCtype[LEFT][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[LEFT][iBC][0], cartMesh.BCrange[LEFT][iBC][1]);
		for (int j = cartMesh.BCrange[LEFT][iBC][0] - 1; j < cartMesh.BCrange[LEFT][iBC][1]; ++j) {
			(*BCedge)[k][0] = j * (iMax + 1);
			(*BCedge)[k][1] = (*BCedge)[k][0] + iMax + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[LEFT][iBC];
			//printf("%ld %ld %ld %ld\n", k, (*BCedge)[k][0], (*BCedge)[k][1], (*BCedge)[k][2]);
			k++;
		}
	}
}

/*
 * create a cartesian, structured mesh
 * read in of all supported mesh types:
 * *.msh, *.emc2, *.cgns
 */
void createMesh(void)
{
	nTrias = nQuads = 0;

	/* create cartesian mesh or read unstructured mesh from file */
	double **vertex;
	long **tria, **quad, **BCedge, *zoneConnect;
	long nVertices, nBCedges, nZones;
	switch (meshType) {
	case CARTESIAN:
		createCartMesh(&vertex, &nVertices, &BCedge, &nBCedges, &quad);
		zoneConnect = calloc(nVertices, sizeof(long));
		break;
	case UNSTRUCTURED:
		// TODO: read in unstructured mesh
		break;
	}

	/* generate mesh information */
	nElems = nTrias + nQuads;
	nInnerSides = (3 * nTrias + 4 * nQuads - nBCedges) / 2;
	nBCsides = nBCedges;
	nSides = nInnerSides + nBCedges;
	nNodes = 0;
}

/*
 * read in and store the mesh parameters from parameter file
 */
void readMesh(void)
{
	meshType = getInt("meshType", "1"); /* default is cartesian */
	switch (meshType) {
	case UNSTRUCTURED:
		printf("| Mesh Type is UNSTRUCTURED\n");

		strcpy(strMeshFormat, getStr("meshFormat", NULL));
		printf("| Mesh File Format is %s\n", strMeshFormat);

		strcpy(strMeshFile, getStr("meshFile", NULL));
		printf("| Mesh File: %s%s\n", strMeshFile, strMeshFormat);

		strcat(strMeshFile, strMeshFormat);
		break;
	case CARTESIAN:
		printf("| Mesh Type is CARTESIAN\n");

		cartMesh.iMax = getInt("nElemsX", NULL);
		cartMesh.jMax = getInt("nElemsY", NULL);

		double *X0 = getDblArray("X0", NDIM, NULL);
		xMin = X0[0];
		yMin = X0[1];
		double *DX = getDblArray("Xmax", NDIM, NULL);
		xMax = DX[0];
		yMax = DX[1];

		/* read boundary conditions */
		cartMesh.nBC = getIntArray("nBCsegments", 2 * NDIM, NULL);
		for (int iSide = 0; iSide < 2 * NDIM; ++iSide) {
			if (cartMesh.nBC[iSide] == 1) {
				cartMesh.BCtype[iSide][0] = getInt("meshBCtype", NULL);
				cartMesh.BCrange[iSide][0][0] = 1;
				if ((iSide % 2) == 0) {
					cartMesh.BCrange[iSide][0][1] = cartMesh.iMax;
				} else {
					cartMesh.BCrange[iSide][0][1] = cartMesh.jMax;
				}
			} else {
				for (int i = 0; i < cartMesh.nBC[i]; ++i) {
					cartMesh.BCtype[iSide][i] = getInt("meshBCtype", NULL);
					cartMesh.BCrange[iSide][i][0] = getInt("BCstart", NULL);
					cartMesh.BCrange[iSide][i][1] = getInt("BCend", NULL);
				}
			}
		}
		break;
	default:
		printf("ERROR: Mesh type can only be unstructured(=0) or cartesian(=1)\n");
		exit(1);
	}
}

/*
 * Initalize the mesh and call read mesh
 */
void initMesh(void)
{
	printf("\nInitializing Mesh:\n");
	readMesh();
	strcat(strcpy(gridFile, strOutFile), "_mesh.cgns");
	createMesh();
	if ((iVisuProg == CGNS) && (!isRestart)) {
		//cgnsWriteMesh();
	}
	dxRef = sqrt(1.0 / (totalArea_q * nElems));
}
