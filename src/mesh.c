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
	//createMesh();
	if ((iVisuProg == CGNS) && (!isRestart)) {
		//cgnsWriteMesh();
	}
	dxRef = sqrt(1.0 / (totalArea_q * nElems));
}
