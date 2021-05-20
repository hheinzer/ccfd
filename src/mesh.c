/** \file
 *
 * \brief Contains all the functions for reading and creating meshes
 *
 * \author hhh
 * \date Tue 24 Mar 2020 12:45:57 PM CET
 */

typedef struct sideList_t sideList_t;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

#include "main.h"
#include "mesh.h"
#include "readInTools.h"
#include "output.h"
#include "timeDiscretization.h"
#include "memTools.h"
#include "initialCondition.h"
#include "cgnslib.h"

/* extern variables */
char strMeshFormat[STRLEN];		/**< mesh format string */
char strMeshFile[STRLEN];		/**< mesh file base name */
char strIniCondFile[STRLEN];		/**< file name of the initial
						conditions file */

cartMesh_t cartMesh;			/**< cartesian mesh structure */
char gridFile[STRLEN];			/**< complete name of the output mesh file */

int meshType;				/**< code for the mesh type */
int meshFormat;				/**< code for the mesh format */

long nNodes;				/**< global number of nodes */

long nElems;				/**< global number of elements */
long nTrias;				/**< global number of triangles */
long nQuads;				/**< global number of quadrangles */

long nSides;				/**< global number of sides */
long nBCsides;				/**< global number of BC sides */
long nInnerSides;			/**< global number of non BC sides */

double totalArea_q;			/**< inverse of the global area of the mesh */
double xMin;				/**< maximum x-direction extension */
double xMax;				/**< minimum x-direction extension */
double yMin;				/**< maximum y-direction extension */
double yMax;				/**< minimum y-direction extension */
double dxRef;				/**< reference x length */

elem_t **elem;				/**< global element pointer array */
side_t **side;				/**< global side pointer array */
side_t **BCside;			/**< global BC side pointer array */

elem_t *firstElem;			/**< pointer to first element */
node_t *firstNode;			/**< pointer to first node */
side_t *firstSide;			/**< pointer to first side */
sidePtr_t *firstBCside;			/**< pointer to first BC side */

/**
 * \brief Helper structure for reading in the sides and deviding them into
 *	BC sides and non-BC sides
 */
struct sideList_t {
	long node[2];			/**< node ID array of the side */
	bool BC;			/**< flag for if the side is a BC side */
	side_t *side;			/**< pointer to the side */
	bool isRotated;			/**< flag for if the side is rotated */
};

/**
 * \brief Compute required vectors for reconstruction
 * \param[in] aElem A pointer to an element
 */
void createReconstructionInfo(elem_t *aElem)
{
	double r11, r12, r22;
	r11 = r12 = r22 = 0.0;

	side_t *aSide = aElem->firstSide;
	while (aSide) {
		aSide->baryBaryVec[X] = aSide->GP[X] - aSide->connection->GP[X];
		aSide->baryBaryVec[Y] = aSide->GP[Y] - aSide->connection->GP[Y];

		aSide->baryBaryDist = sqrt(aSide->baryBaryVec[X] * aSide->baryBaryVec[X] +
				           aSide->baryBaryVec[Y] * aSide->baryBaryVec[Y]);

		r11 += aSide->baryBaryVec[X] * aSide->baryBaryVec[X];
		r12 += aSide->baryBaryVec[X] * aSide->baryBaryVec[Y];
		r22 += aSide->baryBaryVec[Y] * aSide->baryBaryVec[Y];

		aSide = aSide->nextElemSide;
	}

	r11 = sqrt(r11);
	r12 = r12 / r11;
	r22 = sqrt(r22 - r12 * r12);

	aSide = aElem->firstSide;
	while (aSide) {
		double alpha1 = aSide->baryBaryVec[X] / (r11 * r11);
		double alpha2 = (aSide->baryBaryVec[Y] - r12 / r11 *
				 aSide->baryBaryVec[X]) / (r22 * r22);

		aSide->w[X] = alpha1 - r12 / r11 * alpha2;
		aSide->w[Y] = alpha2;

		aSide = aSide->nextElemSide;
	}

	/* gaussian integration points and weights for volume integrals */
	switch (aElem->elemType) {
	case 3:
		aElem->nGP = 3;
		aElem->wGP = malloc(aElem->nGP * sizeof(double));
		if (!aElem->wGP) {
			printf("| ERROR: could not allocate aElem->wGP\n");
			exit(1);
		}

		aElem->xGP = dyn2DdblArray(aElem->nGP, NDIM);

		aSide = aElem->firstSide;
		for (int iGP = 0; iGP < aElem->nGP; ++iGP) {
			aElem->xGP[iGP][X] = aElem->bary[X] + aSide->GP[X];
			aElem->xGP[iGP][Y] = aElem->bary[Y] + aSide->GP[Y];

			aElem->wGP[iGP] = aElem->area / 3.0;

			aSide = aSide->nextElemSide;
		}
		break;
	case 4:
		aElem->nGP = 5;
		aElem->wGP = malloc(aElem->nGP * sizeof(double));
		if (!aElem->wGP) {
			printf("| ERROR: could not allocate aElem->wGP\n");
			exit(1);
		}

		aElem->xGP = dyn2DdblArray(aElem->nGP, NDIM);

		aSide = aElem->firstSide;
		for (int iGP = 0; iGP < aElem->nGP - 1; ++iGP) {
			aElem->xGP[iGP][X] = aElem->bary[X] + aSide->GP[X];
			aElem->xGP[iGP][Y] = aElem->bary[Y] + aSide->GP[Y];

			aElem->wGP[iGP] = aElem->area / 6.0;

			aSide = aSide->nextElemSide;
		}
		aElem->xGP[4][X] = 0.5 * (aElem->node[0]->x[X] + aElem->node[2]->x[X]);
		aElem->xGP[4][Y] = 0.5 * (aElem->node[0]->x[Y] + aElem->node[2]->x[Y]);

		aElem->wGP[4] = aElem->area / 3.0;
		break;
	default:
		printf("| ERROR in createReconstructionInfo:\n");
		printf("| No quadrature rule for elements other than triangles or quadrilaterals\n");
		exit(1);
	}
}

/**
 * \brief Create connection for periodic BCs
 */
void connectPeriodicBC(void)
{
	isPeriodic = false;

	sidePtr_t *aBCside = firstBCside;
	double aGPpos[2], sGPpos[2];
	bool isConnected;
	while (aBCside) {
		side_t *aSide = aBCside->side;
		if ((aSide->BC->BCtype == PERIODIC) &&
		    ((aSide->BC->BCid % 10) == 1)) {
			isPeriodic = true;
			isConnected = false;

			aGPpos[X] = aSide->GP[X] + aSide->elem->bary[X];
			aGPpos[Y] = aSide->GP[Y] + aSide->elem->bary[Y];

			int nPeriodicBC = aSide->BC->BCid / 10;

			sidePtr_t *sBCside = firstBCside;
			while (sBCside) {
				side_t *sSide = sBCside->side;

				/* check for potential connection */
				if ((sSide->BC->BCtype == PERIODIC) &&
				    ((sSide->BC->BCid % 10) == 2) &&
				    ((sSide->BC->BCid / 10) == nPeriodicBC)) {

					/* check if the connection works out */
					sGPpos[X] = sSide->GP[X] + sSide->elem->bary[X];
					sGPpos[Y] = sSide->GP[Y] + sSide->elem->bary[Y];

					if ((fabs(aGPpos[X] + aSide->BC->connection[X] - sGPpos[X]) +
					     fabs(aGPpos[Y] + aSide->BC->connection[Y] - sGPpos[Y])) <= 1e-7) {

						isConnected = true;

						side_t *urSide = aSide->connection;
						side_t *targetSide = sSide->connection;
						urSide->connection = targetSide;
						targetSide->connection = urSide;

						/* reorganise lists: target side has to be removed */
						nSides--;
						if (firstSide == targetSide) {
							firstSide = firstSide->next;
						} else {
							side_t *sSide = firstSide;
							while (sSide) {
								if (sSide->next == targetSide) {
									sSide->next = targetSide->next;
									break;
								}
								sSide = sSide->next;
							}
						}
						break;
					}
				}

				sBCside = sBCside->next;
			}

			if (!isConnected) {
				printf("| ERROR in connectPeriodicBC: No connection found\n");
				printf("| Side GP was: %g %g\n", aGPpos[X], aGPpos[Y]);
				exit(1);
			}
		}

		aBCside = aBCside->next;
	}
}

/**
 * \brief Create side info: normal vector, side length, and ghost cells
 * \param[in] aSide A pointer to a side
 */
void createSideInfo(side_t *aSide)
{
	aSide->n[X] = aSide->node[1]->x[Y] - aSide->node[0]->x[Y];
	aSide->n[Y] = aSide->node[0]->x[X] - aSide->node[1]->x[X];

	aSide->len = sqrt(aSide->n[X] * aSide->n[X] + aSide->n[Y] * aSide->n[Y]);

	aSide->n[X] /= aSide->len;
	aSide->n[Y] /= aSide->len;

	double GP[NDIM];
	GP[X] = 0.5 * (aSide->node[0]->x[X] + aSide->node[1]->x[X]);
	GP[Y] = 0.5 * (aSide->node[0]->x[Y] + aSide->node[1]->x[Y]);
	aSide->GP[X] = GP[X] - aSide->elem->bary[X];
	aSide->GP[Y] = GP[Y] - aSide->elem->bary[Y];

	/* compute dot-product of vector from barycenter to side midpoint with
	 * normal vector, if it is less than zero, the direction must be
	 * switched */
	if ((aSide->n[X] * aSide->GP[X] + aSide->n[Y] * aSide->GP[Y]) < 0.0) {
		aSide->n[X] *= -1.0;
		aSide->n[Y] *= -1.0;
	}

	/* compute barycenter of ghost cell in case it is a boundary side */
	if (aSide->connection->BC) {
		double tmp = 2.0 * fabs(aSide->GP[X] * aSide->n[X]
				      + aSide->GP[Y] * aSide->n[Y]);
		aSide->connection->elem->bary[X] = aSide->elem->bary[X]
			+ tmp * aSide->n[X];
		aSide->connection->elem->bary[Y] = aSide->elem->bary[Y]
			+ tmp * aSide->n[Y];
	}

	aSide->connection->GP[X] = GP[X] - aSide->connection->elem->bary[X];
	aSide->connection->GP[Y] = GP[Y] - aSide->connection->elem->bary[Y];
	aSide->connection->n[X]  = - aSide->n[X];
	aSide->connection->n[Y]  = - aSide->n[Y];
	aSide->connection->len   = aSide->len;
}

/**
 * \brief Compute the cell specific values: barycenter, area, projection length
 *	of an element
 * \param[in] aElem A pointer to an element
 */
void createElemInfo(elem_t *aElem)
{
	/* compute barycenter */
	aElem->bary[X] = aElem->bary[Y] = 0.0;
	for (int i = 0; i < aElem->elemType; ++i) {
		aElem->bary[X] += aElem->node[i]->x[X] / aElem->elemType;
		aElem->bary[Y] += aElem->node[i]->x[Y] / aElem->elemType;
	}

	/* element area and its inverse */
	double vec[2] = {0}, n[4][2] = {0}, len[4] = {0};
	for (int i = 0; i < aElem->elemType; ++i) {
		int j = (i + 1) % aElem->elemType;
		vec[X] = aElem->node[j]->x[X] - aElem->node[i]->x[X];
		vec[Y] = aElem->node[j]->x[Y] - aElem->node[i]->x[Y];

		len[i] = sqrt(vec[X] * vec[X] + vec[Y] * vec[Y]);

		n[i][X] = aElem->node[j]->x[Y] - aElem->node[i]->x[Y];
		n[i][Y] = aElem->node[i]->x[X] - aElem->node[j]->x[X];

		double tmp = sqrt(n[i][X] * n[i][X] + n[i][Y] * n[i][Y]);
		n[i][X] /= tmp;
		n[i][Y] /= tmp;
	}

	aElem->area = 0.5 * len[0] * len[1] * fabs(n[0][X]*n[1][Y] - n[1][X]*n[0][Y]);
	if (aElem->elemType == 4) {
		aElem->area +=
			0.5 * len[2] * len[3] * fabs(n[2][X]*n[3][Y] - n[3][X]*n[2][Y]);
	}
	aElem->areaq = 1.0 / aElem->area;

	if (totalArea_q != 0.0) {
		totalArea_q = 1.0 / (1.0 / totalArea_q + aElem->area);
	} else {
		totalArea_q = aElem->areaq;
	}

	/* projection of cell onto x- and y-axis */
	aElem->sx = aElem->sy = 0;
	for (int i = 0; i < aElem->elemType; ++i) {
		aElem->sx += fabs(n[i][X]) * len[i] / 2.0;
		aElem->sy += fabs(n[i][Y]) * len[i] / 2.0;
	}
}

/** \brief Compare two elements of the `sideList` list, used for sorting the list
 *
 * Sort sideList array in order to retrieve the connectivity info more
 * efficiently. Basically we want to sort the sides on three different levels:
 * first, by their first node ID, secondly, by their second node ID, and thirdly
 * by if they are a boundary condition side. The result might look something
 * like this:
 * 3 3 1       1 2 0
 * 2 1 1       1 3 0
 * 1 2 0  -->  2 1 1
 * 3 1 0       2 2 1
 * 2 2 1       3 3 0
 * 3 3 0       3 3 1
 *
 * \param[in] a Pointer an element of the side list
 * \param[in] b Pointer an element of the side list
 * \return The difference of a and b, if the values are the same, the move to
 *	next sorting metric
 */
int compare(const void *a, const void *b)
{
	sideList_t *A = (sideList_t *)a;
	sideList_t *B = (sideList_t *)b;
	if (A->node[0] == B->node[0]) {
		if (A->node[1] == B->node[1]) {
			return (A->BC ? 1 : 0) - (B->BC ? 1 : 0);
		} else {
			return A->node[1] - B->node[1];
		}
	} else {
		return A->node[0] - B->node[0];
	}
}

/**
 * \brief Create a cartesian mesh
 * \param[in,out] vertex Pointer to 2D array, used for the vertices
 * \param[in,out] nVertices Pointer to the number of total vertices in `vertex`
 * \param[in,out] BCedge Pointer to 2D array, used for the BC edges
 * \param[in,out] nBCedges Pointer to the number of total BC edges
 * \param[in,out] quad Pointer to a 2D array, used for the quadrangles
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
			k++;
		}
	}

	/* build BCedges */
	k = 0;
	/* bottom side */
	printf("| Bottom Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[BOTTOM]);
	for (int iBC = 0; iBC < cartMesh.nBC[BOTTOM]; ++iBC) {
		printf("|   BC Type : %d\n", cartMesh.BCtype[BOTTOM][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[BOTTOM][iBC][0], cartMesh.BCrange[BOTTOM][iBC][1]);
		for (int i = cartMesh.BCrange[BOTTOM][iBC][0] - 1; i < cartMesh.BCrange[BOTTOM][iBC][1]; ++i) {
			(*BCedge)[k][0] = i;
			(*BCedge)[k][1] = i + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[BOTTOM][iBC];
			k++;
		}
	}
	/* right side */
	printf("| Right Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[RIGHT]);
	for (int iBC = 0; iBC < cartMesh.nBC[RIGHT]; ++iBC) {
		printf("|   BC Type : %d\n", cartMesh.BCtype[RIGHT][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[RIGHT][iBC][0], cartMesh.BCrange[RIGHT][iBC][1]);
		for (int j = cartMesh.BCrange[RIGHT][iBC][0] - 1; j < cartMesh.BCrange[RIGHT][iBC][1]; ++j) {
			(*BCedge)[k][0] = (iMax + 1) * (j + 1) - 1;
			(*BCedge)[k][1] = (iMax + 1) * (j + 2) - 1;
			(*BCedge)[k][2] = cartMesh.BCtype[RIGHT][iBC];
			k++;
		}
	}
	/* top side */
	printf("| Top Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[TOP]);
	for (int iBC = 0; iBC < cartMesh.nBC[TOP]; ++iBC) {
		printf("|   BC Type : %d\n", cartMesh.BCtype[TOP][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[TOP][iBC][0], cartMesh.BCrange[TOP][iBC][1]);
		for (int i = cartMesh.BCrange[TOP][iBC][0] - 1; i < cartMesh.BCrange[TOP][iBC][1]; ++i) {
			(*BCedge)[k][0] = (iMax + 1) * jMax + i;
			(*BCedge)[k][1] = (*BCedge)[k][0] + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[TOP][iBC];
			k++;
		}
	}
	/* left side */
	printf("| Left Side:\n");
	printf("|   # of BCs: %d\n", cartMesh.nBC[LEFT]);
	for (int iBC = 0; iBC < cartMesh.nBC[LEFT]; ++iBC) {
		printf("|   BC Type : %d\n", cartMesh.BCtype[LEFT][iBC]);
		printf("|   BC Range: %d -- %d\n", cartMesh.BCrange[LEFT][iBC][0], cartMesh.BCrange[LEFT][iBC][1]);
		for (int j = cartMesh.BCrange[LEFT][iBC][0] - 1; j < cartMesh.BCrange[LEFT][iBC][1]; ++j) {
			(*BCedge)[k][0] = j * (iMax + 1);
			(*BCedge)[k][1] = (*BCedge)[k][0] + iMax + 1;
			(*BCedge)[k][2] = cartMesh.BCtype[LEFT][iBC];
			k++;
		}
	}

	free(cartMesh.nBC);
}

/**
 * \brief Read in a gmsh mesh file
 * \param[in] fileName Name of the mesh file
 * \param[in,out] vertex Pointer to 2D array, used for the vertices
 * \param[in,out] nVertices Pointer to the number of total vertices in `vertex`
 * \param[in,out] BCedge Pointer to 2D array, used for the BC edges
 * \param[in,out] nBCedges Pointer to the number of total BC edges
 * \param[in,out] tria Pointer to a 2D array, used for the triangles
 * \param[in,out] quad Pointer to a 2D array, used for the quadrangles
 */
void readGmsh(char fileName[STRLEN], double ***vertex, long *nVertices, long ***BCedge,
		long *nBCedges, long ***tria, long ***quad)
{
	/* open the mesh file */
	FILE *meshFile = fopen(fileName, "r");
	if (!meshFile) {
		printf("| ERROR: Could not find Mesh File '%s'\n", fileName);
		exit(1);
	}

	char line[STRLEN];

	/* read in mesh format */
	int mshFmt;
	rewind(meshFile);
	while (fgets(line, sizeof(line), meshFile)) {
		if (!strcmp(line, "$MeshFormat\n")) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%d", &mshFmt);
			break;
		}
	}

	switch (mshFmt) {
	case 2: {
		/* read in number of nodes */
		*nVertices = 0;
		rewind(meshFile);
		while (fgets(line, sizeof(line), meshFile)) {
			if (!strcmp(line, "$Nodes\n")) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%ld", nVertices);
				break;
			}
		}
		if (*nVertices == 0) {
			printf("| ERROR: No Nodes in Mesh file\n");
			exit(1);
		}

		/* read in nodes */
		long id;
		double x, y;
		*vertex = dyn2DdblArray(*nVertices, 2);
		for (long iVert = 0; iVert < *nVertices; ++iVert) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%ld %lg %lg", &id, &x, &y);
			if (iVert == id - 1) {
				(*vertex)[iVert][X] = x;
				(*vertex)[iVert][Y] = y;
			} else {
				printf("| ERROR: NodeID %ld does not match Node Position\n", id);
				exit(1);
			}
		}

		*nBCedges = nTrias = nQuads = 0;

		/* read in number of elements */
		long nElem = 0;
		rewind(meshFile);
		while (fgets(line, sizeof(line), meshFile)) {
			if (!strcmp(line, "$Elements\n")) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%ld", &nElem);
				break;
			}
		}
		if (nElem == 0) {
			printf("| ERROR: No Elements in Mesh file\n");
			exit(1);
		}

		/* allocate space for the elements, more than actually necessary,
		 * because the number of individual elements is not known, will be
		 * freeed later on though */
		long **BCedgeTmp = dyn2DintArray(nElem, 3);
		long **triaTmp = dyn2DintArray(nElem, 4);
		long **quadTmp = dyn2DintArray(nElem, 5);
		if (!BCedgeTmp || !triaTmp || !quadTmp) {
			printf("| ERROR: Not enough memory available\n");
			exit(1);
		}

		/* read in elements */
		int type, nTags, pGroup;
		char *token, sep[] = " ";
		for (long iElem = 0; iElem < nElem; ++iElem) {
			fgets(line, sizeof(line), meshFile);

			/* read element id */
			token = strtok(line, sep);
			if (token) {
				id = strtol(token, NULL, 10);
			} else {
				printf("| ERROR: Reading error in gmsh file\n");
				exit(1);
			}

			/* read element type */
			token = strtok(NULL, sep);
			if (token) {
				type = strtol(token, NULL, 10);
			} else {
				printf("| ERROR: Reading error in gmsh file\n");
				exit(1);
			}

			/* read number of integer tags */
			token = strtok(NULL, sep);
			if (token) {
				nTags = strtol(token, NULL, 10);
			} else {
				printf("| ERROR: Reading error in gmsh file\n");
				exit(1);
			}

			/* read physical group, first integer tag */
			token = strtok(NULL, sep);
			if (token) {
				pGroup = strtol(token, NULL, 10);
			} else {
				printf("| ERROR: Reading error in gmsh file\n");
				exit(1);
			}

			/* discard rest of tags */
			for (int iTag = 1; iTag < nTags; ++iTag) {
				token = strtok(NULL, sep);
				if (!token) {
					printf("| ERROR: Reading error in gmsh file\n");
					exit(1);
				}
			}

			/* handle different element types */
			switch (type) {
			case 1:
				/* line */
				if (pGroup > 100) {
					/* boundary condition */
					BCedgeTmp[*nBCedges][2] = pGroup;

					/* read two nodes */
					for (int i = 0; i < 2; ++i) {
						token = strtok(NULL, sep);
						if (!token) {
							printf("| ERROR: Reading error in gmsh file\n");
							exit(1);
						}
						BCedgeTmp[*nBCedges][i] = strtol(token, NULL, 10) - 1;
					}

					(*nBCedges)++;
				}
				break;
			case 2:
				/* triangle */
				triaTmp[nTrias][3] = pGroup;

				/* read three nodes */
				for (int i = 0; i < 3; ++i) {
					token = strtok(NULL, sep);
					if (!token) {
						printf("| ERROR: Reading error in gmsh file\n");
						exit(1);
					}
					triaTmp[nTrias][i] = strtol(token, NULL, 10) - 1;
				}

				nTrias++;

				break;
			case 3:
				/* quadrilateral */
				quadTmp[nQuads][4] = pGroup;

				/* read three nodes */
				for (int i = 0; i < 4; ++i) {
					token = strtok(NULL, sep);
					if (!token) {
						printf("| ERROR: Reading error in gmsh file\n");
						exit(1);
					}
					quadTmp[nQuads][i] = strtol(token, NULL, 10) - 1;
				}

				nQuads++;

				break;
			}
		}

		/* allocate right amount of space and free temporary space */
		*BCedge = dyn2DintArray(*nBCedges, 3);
		if (!*BCedge) {
			printf("| ERROR: Not enough Memory\n");
			exit(1);
		} else {
			for (long iBCedge = 0; iBCedge < *nBCedges; ++iBCedge) {
				(*BCedge)[iBCedge][0] = BCedgeTmp[iBCedge][0];
				(*BCedge)[iBCedge][1] = BCedgeTmp[iBCedge][1];
				(*BCedge)[iBCedge][2] = BCedgeTmp[iBCedge][2];
			}
		}

		free(BCedgeTmp);

		if (nTrias > 0) {
			*tria = dyn2DintArray(nTrias, 4);
			if (!*tria) {
				printf("| ERROR: Not enough Memory\n");
				exit(1);
			} else {
				for (long iTria = 0; iTria < nTrias; ++iTria) {
					(*tria)[iTria][0] = triaTmp[iTria][0];
					(*tria)[iTria][1] = triaTmp[iTria][1];
					(*tria)[iTria][2] = triaTmp[iTria][2];
					(*tria)[iTria][3] = triaTmp[iTria][3];
				}
			}
		}

		free(triaTmp);

		if (nQuads > 0) {
			*quad = dyn2DintArray(nQuads, 5);
			if (!*quad) {
				printf("| ERROR: Not enough Memory\n");
				exit(1);
			} else {
				for (long iQuad = 0; iQuad < nQuads; ++iQuad) {
					(*quad)[iQuad][0] = quadTmp[iQuad][0];
					(*quad)[iQuad][1] = quadTmp[iQuad][1];
					(*quad)[iQuad][2] = quadTmp[iQuad][2];
					(*quad)[iQuad][3] = quadTmp[iQuad][3];
					(*quad)[iQuad][4] = quadTmp[iQuad][4];
				}
			}
		}

		free(quadTmp);

		break;
	}
	case 4: {
		/* read curve number of entity tags */
		long nPEntities, nCEntities;
		rewind(meshFile);
		while (fgets(line, sizeof(line), meshFile)) {
			if (!strcmp(line, "$Entities\n")) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%ld %ld", &nPEntities, &nCEntities);
				break;
			}
		}
		if (nCEntities == 0) {
			printf("| ERROR: No curves in mesh file\n");
			exit(1);
		}

		/* read curve entity tags */
		int pTag[nCEntities];
		int iCentitiy = 0;
		char *token, sep[] = " ";
		for (int l = 0; l < nPEntities; l++) {
			/* skip point tags */
			fgets(line, sizeof(line), meshFile);
		}
		for (int l = 0; l < nCEntities; l++) {
			fgets(line, sizeof(line), meshFile);

			token = strtok(line, sep);
			for (int k = 0; k < 7; ++k) {
				token = strtok(NULL, sep);
			}
			int numPtags = strtol(token, NULL, 10);

			if (numPtags) {
				token = strtok(NULL, sep);
				pTag[iCentitiy++] = strtol(token, NULL, 10);
			}
		}

		/* read in number of entity block and number of nordes */
		*nVertices = 0;
		int nEblocks;
		rewind(meshFile);
		while (fgets(line, sizeof(line), meshFile)) {
			if (!strcmp(line, "$Nodes\n")) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%d %ld", &nEblocks, nVertices);
				break;
			}
		}
		if (*nVertices == 0) {
			printf("| ERROR: No Nodes in Mesh file\n");
			exit(1);
		}

		/* read in nodes */
		long id[*nVertices];
		long nNinB, iVert = 0, iId = 0;
		double x, y;
		*vertex = dyn2DdblArray(*nVertices, 2);
		for (int e = 0; e < nEblocks; ++e) {
			fgets(line, sizeof(line), meshFile);

			/* read number of nodes in block */
			sscanf(line, "%*d %*d %*d %ld", &nNinB);

			/* read ids */
			for (long i = 0; i < nNinB; ++i) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%ld", &id[iId++]);
			}

			/* read coordinates */
			for (long i = 0; i < nNinB; ++i) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%lg %lg", &x, &y);

				if (iVert == id[iVert] - 1) {
					(*vertex)[iVert  ][X] = x;
					(*vertex)[iVert++][Y] = y;
				} else {
					printf("| ERROR: NodeID %ld does not match Node Position\n", id[iVert]);
					exit(1);
				}
			}
		}

		*nBCedges = nTrias = nQuads = 0;

		/* read in number of element blocks and number of elements */
		long nElem = 0;
		rewind(meshFile);
		while (fgets(line, sizeof(line), meshFile)) {
			if (!strcmp(line, "$Elements\n")) {
				fgets(line, sizeof(line), meshFile);
				sscanf(line, "%d %ld", &nEblocks, &nElem);
				break;
			}
		}
		if (nElem == 0) {
			printf("| ERROR: No Elements in Mesh file\n");
			exit(1);
		}

		/* allocate space for the elements, more than actually necessary,
		 * because the number of individual elements is not known, will be
		 * freeed later on though */
		long **BCedgeTmp = dyn2DintArray(nElem, 3);
		long **triaTmp = dyn2DintArray(nElem, 4);
		long **quadTmp = dyn2DintArray(nElem, 5);
		if (!BCedgeTmp || !triaTmp || !quadTmp) {
			printf("| ERROR: Not enough memory available\n");
			exit(1);
		}

		/* read in elements */
		int type;
		for (int e = 0; e < nEblocks; ++e) {
			fgets(line, sizeof(line), meshFile);

			/* read element type and number of elements in block */
			sscanf(line, "%*d %*d %d %ld", &type, &nNinB);

			/* handle different element types */
			switch (type) {
			case 1:
				/* line */
				for (long i = 0; i < nNinB; ++i) {
					/* boundary condition */
					BCedgeTmp[*nBCedges][2] = pTag[e];

					/* read two nodes */
					fgets(line, sizeof(line), meshFile);
					sscanf(line, "%*d %ld %ld",
							&BCedgeTmp[*nBCedges][0],
							&BCedgeTmp[*nBCedges][1]);

					BCedgeTmp[*nBCedges][0]--;
					BCedgeTmp[*nBCedges][1]--;

					(*nBCedges)++;
				}
				break;
			case 2:
				/* triangle */
				for (long i = 0; i < nNinB; ++i) {
					triaTmp[nTrias][3] = 1;

					/* read three nodes */
					fgets(line, sizeof(line), meshFile);
					sscanf(line, "%*d %ld %ld %ld",
							&triaTmp[nTrias][0],
							&triaTmp[nTrias][1],
							&triaTmp[nTrias][2]);

					triaTmp[nTrias][0]--;
					triaTmp[nTrias][1]--;
					triaTmp[nTrias][2]--;

					nTrias++;
				}

				break;
			case 3:
				/* quadrilateral */
				for (long i = 0; i < nNinB; ++i) {
					quadTmp[nQuads][4] = 1;

					/* read three nodes */
					fgets(line, sizeof(line), meshFile);
					sscanf(line, "%*d %ld %ld %ld %ld",
							&quadTmp[nQuads][0],
							&quadTmp[nQuads][1],
							&quadTmp[nQuads][2],
							&quadTmp[nQuads][3]);

					quadTmp[nQuads][0]--;
					quadTmp[nQuads][1]--;
					quadTmp[nQuads][2]--;
					quadTmp[nQuads][3]--;

					nQuads++;
				}

				break;
			}
		}

		/* allocate right amount of space and free temporary space */
		*BCedge = dyn2DintArray(*nBCedges, 3);
		if (!*BCedge) {
			printf("| ERROR: Not enough Memory\n");
			exit(1);
		} else {
			for (long iBCedge = 0; iBCedge < *nBCedges; ++iBCedge) {
				(*BCedge)[iBCedge][0] = BCedgeTmp[iBCedge][0];
				(*BCedge)[iBCedge][1] = BCedgeTmp[iBCedge][1];
				(*BCedge)[iBCedge][2] = BCedgeTmp[iBCedge][2];
			}
		}

		free(BCedgeTmp);

		if (nTrias > 0) {
			*tria = dyn2DintArray(nTrias, 4);
			if (!*tria) {
				printf("| ERROR: Not enough Memory\n");
				exit(1);
			} else {
				for (long iTria = 0; iTria < nTrias; ++iTria) {
					(*tria)[iTria][0] = triaTmp[iTria][0];
					(*tria)[iTria][1] = triaTmp[iTria][1];
					(*tria)[iTria][2] = triaTmp[iTria][2];
					(*tria)[iTria][3] = triaTmp[iTria][3];
				}
			}
		}
		free(triaTmp);

		if (nQuads > 0) {
			*quad = dyn2DintArray(nQuads, 5);
			if (!*quad) {
				printf("| ERROR: Not enough Memory\n");
				exit(1);
			} else {
				for (long iQuad = 0; iQuad < nQuads; ++iQuad) {
					(*quad)[iQuad][0] = quadTmp[iQuad][0];
					(*quad)[iQuad][1] = quadTmp[iQuad][1];
					(*quad)[iQuad][2] = quadTmp[iQuad][2];
					(*quad)[iQuad][3] = quadTmp[iQuad][3];
					(*quad)[iQuad][4] = quadTmp[iQuad][4];
				}
			}
		}
		free(quadTmp);

		break;
	}
	default:
		printf("| ERROR: Wrong gmsh Mesh Format '%d'\n", mshFmt);
		exit(1);
	}

	fclose(meshFile);

	printf("| %7ld Nodes read\n", *nVertices);
	printf("| %7ld Triangles read\n", nTrias);
	printf("| %7ld Quadrangles read\n", nQuads);
	printf("| %7ld Boundary Edges read\n", *nBCedges);
}

/**
 * \brief Read in EMC2 mesh file
 * \param[in] fileName Name of the mesh file
 * \param[in,out] vertex Pointer to 2D array, used for the vertices
 * \param[in,out] nVertices Pointer to the number of total vertices in `vertex`
 * \param[in,out] BCedge Pointer to 2D array, used for the BC edges
 * \param[in,out] nBCedges Pointer to the number of total BC edges
 * \param[in,out] tria Pointer to a 2D array, used for the triangles
 * \param[in,out] quad Pointer to a 2D array, used for the quadrangles
 */
void readEMC2(char fileName[STRLEN], double ***vertex, long *nVertices, long ***BCedge,
		long *nBCedges, long ***tria, long ***quad)
{
	/* open the mesh file */
	FILE *meshFile = fopen(fileName, "r");
	if (!meshFile) {
		printf("| ERROR: Could not find Mesh File '%s'\n", fileName);
		exit(1);
	}

	char line[STRLEN];

	/* read number of vertices */
	rewind(meshFile);
	while (fgets(line, sizeof(line), meshFile)) {
		if (strstr(line, "Vertices")) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%ld", nVertices);
			break;
		} else if (strstr(line, "End")) {
			break;
		}
	}

	/* allocate memory for vertices */
	*vertex = dyn2DdblArray(*nVertices, 2);
	for (long iVert = 0; iVert < *nVertices; ++iVert) {
		fgets(line, sizeof(line), meshFile);
		double x, y;
		sscanf(line, "%lg %lg", &x, &y);
		(*vertex)[iVert][X] = x;
		(*vertex)[iVert][Y] = y;
	}

	/* read number of edges */
	long nEdges;
	rewind(meshFile);
	while (fgets(line, sizeof(line), meshFile)) {
		if (strstr(line, "Edges")) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%ld", &nEdges);
			break;
		} else if (strstr(line, "End")) {
			break;
		}
	}

	/* read edge information, sort out edges with ID 0 */
	(*nBCedges) = 0;
	long **tmpEdges = dyn2DintArray(nEdges, 3);
	for (long iEdge = 0; iEdge < nEdges; ++iEdge) {
		fgets(line, sizeof(line), meshFile);

		long tmp;
		sscanf(line, "%*d %*d %ld", &tmp);

		if (tmp != 0) {
			long tmp0, tmp1, tmp2;
			sscanf(line, "%ld %ld %ld", &tmp0, &tmp1, &tmp2);
			tmpEdges[*nBCedges][0] = tmp0;
			tmpEdges[*nBCedges][1] = tmp1;
			tmpEdges[*nBCedges][2] = tmp2;
			(*nBCedges)++;
		}
	}

	*BCedge = dyn2DintArray(*nBCedges, 3);
	for (long iEdge = 0; iEdge < *nBCedges; ++iEdge) {
		(*BCedge)[iEdge][0] = tmpEdges[iEdge][0] - 1;
		(*BCedge)[iEdge][1] = tmpEdges[iEdge][1] - 1;
		(*BCedge)[iEdge][2] = tmpEdges[iEdge][2];
	}

	free(tmpEdges);

	/* read number of triangles */
	nTrias = 0;
	rewind(meshFile);
	while (fgets(line, sizeof(line), meshFile)) {
		if (strstr(line, "Triangles")) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%ld", &nTrias);
			break;
		} else if (strstr(line, "End")) {
			break;
		}
	}

	/* read in triangles */
	if (nTrias > 0) {
		*tria = dyn2DintArray(nTrias, 4);
		for (long iTria = 0; iTria < nTrias; ++iTria) {
			fgets(line, sizeof(line), meshFile);
			long tmp0, tmp1, tmp2, tmp3;
			sscanf(line, "%ld %ld %ld %ld", &tmp0, &tmp1, &tmp2, &tmp3);
			(*tria)[iTria][0] = tmp0 - 1;
			(*tria)[iTria][1] = tmp1 - 1;
			(*tria)[iTria][2] = tmp2 - 1;
			(*tria)[iTria][3] = tmp3;
		}
	}

	/* read number of quadrangles */
	nQuads = 0;
	rewind(meshFile);
	while (fgets(line, sizeof(line), meshFile)) {
		if (strstr(line, "Quadrangles")) {
			fgets(line, sizeof(line), meshFile);
			sscanf(line, "%ld", &nQuads);
			break;
		} else if (strstr(line, "End")) {
			break;
		}
	}

	/* read in quadrangles */
	if (nQuads > 0) {
		*quad = dyn2DintArray(nQuads, 5);
		for (long iQuad = 0; iQuad < nQuads; ++iQuad) {
			fgets(line, sizeof(line), meshFile);
			long tmp0, tmp1, tmp2, tmp3, tmp4;
			sscanf(line, "%ld %ld %ld %ld %ld",
					&tmp0, &tmp1, &tmp2, &tmp3, &tmp4);
			(*quad)[iQuad][0] = tmp0 - 1;
			(*quad)[iQuad][1] = tmp1 - 1;
			(*quad)[iQuad][2] = tmp2 - 1;
			(*quad)[iQuad][3] = tmp3 - 1;
			(*quad)[iQuad][4] = tmp4;
		}
	}

	fclose(meshFile);

	printf("| %7ld Nodes read\n", *nVertices);
	printf("| %7ld Triangles read\n", nTrias);
	printf("| %7ld Quadrangles read\n", nQuads);
	printf("| %7ld Boundary Edges read\n", *nBCedges);
}

/**
 * \brief Read in a CGNS mesh
 * \param[in] fileName Name of the mesh file
 * \param[in,out] vertex Pointer to 2D array, used for the vertices
 * \param[in,out] nVertices Pointer to the number of total vertices in `vertex`
 * \param[in,out] BCedge Pointer to 2D array, used for the BC edges
 * \param[in,out] nBCedges Pointer to the number of total BC edges
 * \param[in,out] tria Pointer to a 2D array, used for the triangles
 * \param[in,out] quad Pointer to a 2D array, used for the quadrangles
 */
void readCGNS(char fileName[STRLEN], double ***vertex, long *nVertices, long ***BCedge,
		long *nBCedges, long ***tria, long ***quad)
{
	int indexFile, indexBase = 1, indexZone = 1;
	/* open CGNS file, only one base and one zone */
	if (cg_open(fileName, CG_MODE_READ, &indexFile))
		cg_error_exit();

	/* get nodes */
	char zoneName[32];
	cgsize_t iSize[3];
	if (cg_zone_read(indexFile, indexBase, indexZone, zoneName, iSize))
		cg_error_exit();

	*nVertices = iSize[0];
	*vertex = dyn2DdblArray(*nVertices, 2);

	cgsize_t rMin[1] = {1}, rMax[2] = {*nVertices};
	double x[*nVertices], y[*nVertices];
	if (cg_coord_read(indexFile, indexBase, indexZone, "CoordinateX",
				RealDouble, rMin, rMax, x))
		cg_error_exit();
	if (cg_coord_read(indexFile, indexBase, indexZone, "CoordinateY",
				RealDouble, rMin, rMax, y))
		cg_error_exit();

	for (long iVert = 0; iVert < *nVertices; ++iVert) {
		(*vertex)[iVert][X] = x[iVert];
		(*vertex)[iVert][Y] = y[iVert];
	}

	/* get number of sections */
	int nSections;
	if (cg_nsections(indexFile, indexBase, indexZone, &nSections))
		cg_error_exit();

	cgsize_t *elems;
	nTrias = nQuads = 0;
	for (int indexSection = 1; indexSection <= nSections; ++indexSection) {
		/* get number of elements in section */
		char sectionName[32];
		ElementType_t sectionType;
		cgsize_t elStart, elEnd;
		int nbndry, parentFlag;
		if (cg_section_read(indexFile, indexBase, indexZone, indexSection,
				sectionName, &sectionType, &elStart, &elEnd,
				&nbndry, &parentFlag))
			cg_error_exit();

		/* get number of entries in element connectivity */
		cgsize_t dataSize;
		if (cg_ElementDataSize(indexFile, indexBase, indexZone,
				indexSection, &dataSize))
			cg_error_exit();

		elems = malloc(dataSize * sizeof(cgsize_t));
		if (!elems) {
			printf("| ERROR: could not allocate elems\n");
			exit(1);
		}

		/* read element connectivity */
		if (cg_elements_read(indexFile, indexBase, indexZone, indexSection,
				elems, NULL))
			cg_error_exit();

		/* handle different section types */
		switch (sectionType) {
		case BAR_2: {
			*nBCedges = (long)dataSize / 2;

			*BCedge = dyn2DintArray(*nBCedges, 3);

			long iBCedge = 0;
			for (long iElem = 0; iElem < (long)dataSize;) {
				(*BCedge)[iBCedge  ][0] = elems[iElem + 0] - 1;
				(*BCedge)[iBCedge++][1] = elems[iElem + 1] - 1;
				iElem += 2;
			}

			if (cg_nbocos(indexFile, indexBase, indexZone, &nbndry))
				cg_error_exit();

			iBCedge = 0;
			for (cgsize_t iBC = 1; iBC <= nbndry; ++iBC) {
				char bocoName[32];
				BCType_t BCtype;
				PointSetType_t ptsetType;
				cgsize_t nBCpoints, normalListSize;
				int indexNormal, nDataSets;
				DataType_t normalDataType;
				if (cg_boco_info(indexFile, indexBase, indexZone,
						iBC, bocoName, &BCtype, &ptsetType,
						&nBCpoints, &indexNormal, &normalListSize,
						&normalDataType, &nDataSets))
					cg_error_exit();

				for (cgsize_t iBCp = 0; iBCp < nBCpoints; ++iBCp) {
					int BCcode;
					sscanf(bocoName, "%d", &BCcode);
					(*BCedge)[iBCedge++][2] = BCcode;
				}
			}

			break;
		}
		case TRI_3: {
			nTrias = (long)dataSize / 3;

			*tria = dyn2DintArray(nTrias, 4);

			long iTria = 0;
			for (long iElem = 0; iElem < (long)dataSize;) {
				(*tria)[iTria  ][0] = elems[iElem + 0] - 1;
				(*tria)[iTria  ][1] = elems[iElem + 1] - 1;
				(*tria)[iTria  ][2] = elems[iElem + 2] - 1;
				(*tria)[iTria++][3] = indexZone;
				iElem += 3;
			}

			break;
		}
		case QUAD_4: {
			nQuads = (long)dataSize / 4;

			*quad = dyn2DintArray(nQuads, 5);

			long iQuad = 0;
			for (long iElem = 0; iElem < (long)dataSize;) {
				(*quad)[iQuad  ][0] = elems[iElem + 0] - 1;
				(*quad)[iQuad  ][1] = elems[iElem + 1] - 1;
				(*quad)[iQuad  ][2] = elems[iElem + 2] - 1;
				(*quad)[iQuad  ][3] = elems[iElem + 3] - 1;
				(*quad)[iQuad++][4] = indexZone;
				iElem += 4;
			}

			break;
		}
		case MIXED: {
			for (long iElem = 0; iElem < (long)dataSize;) {
				switch (elems[iElem]) {
				case BAR_2:
					(*nBCedges)++;
					iElem += 3;
					break;
				case TRI_3:
					nTrias++;
					iElem += 4;
					break;
				case QUAD_4:
					nQuads++;
					iElem += 5;
					break;
				default:
					printf("| ERROR: Only TRI_3 and QUAD_4 Elements allowed\n");
					exit(1);
				}
			}

			/* allocate connectivity arrays */
			if (!(*BCedge) && (*nBCedges > 0)) {
				*BCedge = dyn2DintArray(*nBCedges, 3);
			}
			if (!(*tria) && (nTrias > 0)) {
				*tria = dyn2DintArray(nTrias, 4);
			}
			if (!(*quad) && (nQuads > 0)) {
				*quad = dyn2DintArray(nQuads, 5);
			}

			long iTria, iQuad, iBCedge;
			iTria = iQuad = iBCedge = 0;
			for (long iElem = 0; iElem < (long)dataSize;) {
				switch (elems[iElem]) {
				case BAR_2:
					(*BCedge)[iBCedge  ][0] = elems[iElem + 1] - 1;
					(*BCedge)[iBCedge++][1] = elems[iElem + 2] - 1;
					iElem += 3;
					break;
				case TRI_3:
					(*tria)[iTria  ][0] = elems[iElem + 1] - 1;
					(*tria)[iTria  ][1] = elems[iElem + 2] - 1;
					(*tria)[iTria  ][2] = elems[iElem + 3] - 1;
					(*tria)[iTria++][3] = indexZone;
					iElem += 4;
					break;
				case QUAD_4:
					(*quad)[iQuad  ][0] = elems[iElem + 1] - 1;
					(*quad)[iQuad  ][1] = elems[iElem + 2] - 1;
					(*quad)[iQuad  ][2] = elems[iElem + 3] - 1;
					(*quad)[iQuad  ][3] = elems[iElem + 4] - 1;
					(*quad)[iQuad++][4] = indexZone;
					iElem += 5;
					break;
				}
			}

			/* add the boundary codes */
			if (*nBCedges > 0) {
				if (cg_nbocos(indexFile, indexBase, indexZone, &nbndry))
					cg_error_exit();

				iBCedge = 0;
				for (cgsize_t iBC = 1; iBC <= nbndry; ++iBC) {
					char bocoName[32];
					BCType_t BCtype;
					PointSetType_t ptsetType;
					cgsize_t nBCpoints, normalListSize;
					int indexNormal, nDataSets;
					DataType_t normalDataType;
					if (cg_boco_info(indexFile, indexBase, indexZone,
							iBC, bocoName, &BCtype, &ptsetType,
							&nBCpoints, &indexNormal, &normalListSize,
							&normalDataType, &nDataSets))
						cg_error_exit();

					for (cgsize_t iBCp = 0; iBCp < nBCpoints; ++iBCp) {
						int BCcode;
						sscanf(bocoName, "%d", &BCcode);
						(*BCedge)[iBCedge++][2] = BCcode;
					}
				}
			}

			break;
		}
		default:
			printf("| ERROR: Only TRI_3 and QUAD_4 Elements allowed\n");
			exit(1);
		}

		nElems = nTrias + nQuads;

		free(elems);
	}

	printf("| %7ld Nodes read\n", *nVertices);
	printf("| %7ld Triangles read\n", nTrias);
	printf("| %7ld Quadrangles read\n", nQuads);
	printf("| %7ld Boundary Edges read\n", *nBCedges);
}

/** \brief Create a cartesian or structured mesh
 *
 * Read in of all supported mesh types:
 *	- *.msh
 *	- *.msh2
 *	- *.msh4
 *	- *.emc2
 *	- *.cgns
 */
void createMesh(void)
{
	nTrias = nQuads = 0;

	/* create cartesian mesh or read unstructured mesh from file */
	double **vertex = NULL;
	long **tria, **quad, **BCedge;
	tria = quad = BCedge = NULL;
	long nVertices, nBCedges;
	switch (meshType) {
	case CARTESIAN:
		createCartMesh(&vertex, &nVertices, &BCedge, &nBCedges, &quad);
		break;
	case UNSTRUCTURED:
		if (!strcmp(strMeshFormat, ".msh") ||
				!strcmp(strMeshFormat, ".msh2") ||
				!strcmp(strMeshFormat, ".msh4")) {
			printf("| Reading gmsh File:\n");

			readGmsh(strMeshFile, &vertex, &nVertices, &BCedge,
					&nBCedges, &tria, &quad);

		} else if (!strcmp(strMeshFormat, ".mesh")) {
			printf("| Reading EMC2 File:\n");

			readEMC2(strMeshFile, &vertex, &nVertices, &BCedge,
					&nBCedges, &tria, &quad);

		} else if (!strcmp(strMeshFormat, ".cgns")) {
			printf("| Reading CGNS File:\n");

			readCGNS(strMeshFile, &vertex, &nVertices, &BCedge,
					&nBCedges, &tria, &quad);

		} else {
			printf("| ERROR: Unknown Mesh Format\n");
			exit(1);
		}
		break;
	}

	/* generate mesh information */
	nElems = nTrias + nQuads;
	nInnerSides = (3 * nTrias + 4 * nQuads - nBCedges) / 2;
	nBCsides = nBCedges;
	nSides = nInnerSides + nBCedges;
	nNodes = 0;

	/* create nodes */
	node_t **vertexPtr = calloc(nVertices, sizeof(node_t *));

	xMin = vertex[0][X];
	xMax = vertex[0][X];
	yMin = vertex[0][Y];
	yMax = vertex[0][Y];

	firstNode = NULL;
	for (long iNode = nVertices - 1; iNode >= 0; --iNode) {
		node_t *aNode = calloc(1, sizeof(node_t));
		if (!aNode) {
			printf("| ERROR: could not allocate aNode\n");
			exit(1);
		}

		aNode->id = iNode;
		aNode->x[X] = vertex[iNode][X];
		aNode->x[Y] = vertex[iNode][Y];

		xMin = fmin(aNode->x[X], xMin);
		xMax = fmax(aNode->x[X], xMax);
		yMin = fmin(aNode->x[Y], yMin);
		yMax = fmax(aNode->x[Y], yMax);

		aNode->next = firstNode;
		firstNode = aNode;

		vertexPtr[iNode] = aNode;

		nNodes++;
	}
	free(vertex);

	/* elements and side pointers */
	sideList_t *sideList = calloc(2 * nSides, sizeof(sideList_t));

	firstElem = NULL;
	long iElem = 0, iSidePtr = 0;

	/* loop over all triangles */
	elem_t *prevElem = NULL;
	for (long iTria = 0; iTria < nTrias; ++iTria) {
		elem_t *aElem = calloc(1, sizeof(elem_t));
		if (!aElem) {
			printf("| ERROR: could not allocate aElem\n");
			exit(1);
		}

		aElem->id = iElem++;
		aElem->elemType = 3;
		aElem->domain = tria[iTria][aElem->elemType];

		if (!firstElem) {
			firstElem = aElem;
			prevElem = aElem;
		} else {
			prevElem->next = aElem;
			prevElem = aElem;
		}
		aElem->next = NULL;

		aElem->node = calloc(aElem->elemType, sizeof(node_t *));
		if (!aElem->node) {
			printf("| ERROR: could not allocate aElem->node\n");
			exit(1);
		}

		for (int iNode = 0; iNode < aElem->elemType; ++iNode) {
			aElem->node[iNode] = vertexPtr[tria[iTria][iNode]];
		}

		aElem->bary[X] = aElem->bary[Y] = 0;
		for (int iNode = 0; iNode < aElem->elemType; ++iNode) {
			aElem->bary[X] += aElem->node[iNode]->x[X] / aElem->elemType;
			aElem->bary[Y] += aElem->node[iNode]->x[Y] / aElem->elemType;
		}

		aElem->firstSide = NULL;
		for (int iSide = 0; iSide < aElem->elemType; ++iSide) {
			side_t *aSide = calloc(1, sizeof(side_t));
			if (!aSide) {
				printf("| ERROR: could not allocate aSide\n");
				exit(1);
			}

			aSide->id = iSidePtr;
			aSide->connection = NULL;
			aSide->nextElemSide = NULL;
			aSide->next = NULL;
			aSide->BC = NULL;
			aSide->node[0] = aElem->node[iSide];
			aSide->node[1] = aElem->node[(iSide + 1) % aElem->elemType];
			aSide->elem = aElem;

			aSide->nextElemSide = aElem->firstSide;
			aElem->firstSide = aSide;

			long iNode1 = aElem->node[iSide]->id;
			long iNode2 = aElem->node[(iSide + 1) % aElem->elemType]->id;

			sideList[iSidePtr].node[0] = fmin(iNode1, iNode2);
			sideList[iSidePtr].node[1] = fmax(iNode1, iNode2);
			sideList[iSidePtr].BC = false;
			if (fmin(iNode1, iNode2) == iNode2) {
				sideList[iSidePtr].isRotated = true;
			} else {
				sideList[iSidePtr].isRotated = false;
			}
			sideList[iSidePtr].side = aSide;

			iSidePtr++;
		}
	}
	if (nTrias > 0) {
		free(tria);
	}

	/* loop over all quadrilaterals */
	for (long iQuad = 0; iQuad < nQuads; ++iQuad) {
		elem_t *aElem = calloc(1, sizeof(elem_t));
		if (!aElem) {
			printf("| ERROR: could not allocate aElem\n");
			exit(0);
		}

		aElem->id = iElem++;
		aElem->elemType = 4;
		aElem->domain = quad[iQuad][aElem->elemType];

		if (!firstElem) {
			firstElem = aElem;
			prevElem = aElem;
		} else {
			prevElem->next = aElem;
			prevElem = aElem;
		}
		aElem->next = NULL;

		aElem->node = calloc(aElem->elemType, sizeof(node_t *));
		if (!aElem->node) {
			printf("| ERROR: could not allocate aElem->node\n");
			exit(0);
		}

		for (int iNode = 0; iNode < aElem->elemType; ++iNode) {
			aElem->node[iNode] = vertexPtr[quad[iQuad][iNode]];
		}

		aElem->bary[X] = aElem->bary[Y] = 0;
		for (int iNode = 0; iNode < aElem->elemType; ++iNode) {
			aElem->bary[X] += aElem->node[iNode]->x[X] / aElem->elemType;
			aElem->bary[Y] += aElem->node[iNode]->x[Y] / aElem->elemType;
		}

		aElem->firstSide = NULL;
		for (int iSide = 0; iSide < aElem->elemType; ++iSide) {
			side_t *aSide = calloc(1, sizeof(side_t));
			if (!aSide) {
				printf("| ERROR: could not allocate aSide\n");
				exit(0);
			}

			aSide->id = iSidePtr;
			aSide->connection = NULL;
			aSide->nextElemSide = NULL;
			aSide->next = NULL;
			aSide->BC = NULL;
			aSide->node[0] = aElem->node[iSide];
			aSide->node[1] = aElem->node[(iSide + 1) % aElem->elemType];
			aSide->elem = aElem;

			aSide->nextElemSide = aElem->firstSide;
			aElem->firstSide = aSide;

			long iNode1 = aElem->node[iSide]->id;
			long iNode2 = aElem->node[(iSide + 1) % aElem->elemType]->id;

			sideList[iSidePtr].node[0] = fmin(iNode1, iNode2);
			sideList[iSidePtr].node[1] = fmax(iNode1, iNode2);
			sideList[iSidePtr].BC = false;
			if (fmin(iNode1, iNode2) == iNode2) {
				sideList[iSidePtr].isRotated = true;
			} else {
				sideList[iSidePtr].isRotated = false;
			}
			sideList[iSidePtr].side = aSide;

			iSidePtr++;
		}

	}
	if (nQuads > 0) {
		free(quad);
	}

	/* sides and connectivity */
	/* save all BCedges into the sideList array */
	for (long iSide = 0; iSide < nBCsides; ++iSide) {
		side_t *aSide = calloc(1, sizeof(side_t));
		if (!aSide) {
			printf("| ERROR: could not allocate aSide\n");
			exit(0);
		}

		aSide->id = iSidePtr;
		aSide->connection = NULL;
		aSide->nextElemSide = NULL;
		aSide->next = NULL;
		aSide->node[0] = vertexPtr[BCedge[iSide][0]];
		aSide->node[1] = vertexPtr[BCedge[iSide][1]];

		elem_t *aElem = calloc(1, sizeof(elem_t));
		if (!aElem) {
			printf("| ERROR: could not allocate aElem\n");
			exit(0);
		}

		aElem->id = -1;
		aSide->elem = aElem;

		long iNode1 = vertexPtr[BCedge[iSide][0]]->id;
		long iNode2 = vertexPtr[BCedge[iSide][1]]->id;

		sideList[iSidePtr].node[0] = fmin(iNode1, iNode2);
		sideList[iSidePtr].node[1] = fmax(iNode1, iNode2);
		sideList[iSidePtr].BC = true;
		if (fmin(iNode1, iNode2) == iNode2) {
			sideList[iSidePtr].isRotated = true;
		} else {
			sideList[iSidePtr].isRotated = false;
		}

		int BCtype = BCedge[iSide][2] / 100;
		int BCid = BCedge[iSide][2] % 100;

		boundary_t *aBC = firstBC;
		while (aBC) {
			if ((aBC->BCtype == BCtype) && (aBC->BCid == BCid)) {
				aSide->BC = aBC;
				break;
			}
			aBC = aBC->next;
		}
		if ((BCtype != 0) && (!aBC)) {
			printf("| ERROR: BC %d from Gridgen Meshfile not defined in Parameter File\n", BCtype * 100 + BCid);
			exit(1);
		}

		sideList[iSidePtr].side = aSide;

		iSidePtr++;
	}
	free(BCedge);
	free(vertexPtr);

	/* sort sideList */
	qsort(sideList, 2 * nSides, sizeof(sideList[0]), compare);

	/* initialize side lists in mesh */
	firstSide = NULL;
	firstBCside = NULL;
	for (long iSide = 0; iSide < 2 * nSides; iSide += 2) {
		side_t *aSide = sideList[iSide].side;
		side_t *bSide = sideList[iSide + 1].side;

		/* grid file read check */
		if ((sideList[iSide].node[0] - sideList[iSide + 1].node[0] +
		     sideList[iSide].node[1] - sideList[iSide + 1].node[1])
				!= 0) {
			printf("| ERROR: Problem in Connect\n");
			exit(1);
		}

		/* save sides into global element array */
		if ((aSide->BC) || (bSide->BC)) {
			/* boundary side: make sure aSide is the physical side
			 * and bSide the ghost side and save only aSide into
			 * list */
			if (aSide->BC) {
				side_t *cSide = aSide;
				aSide = bSide;
				bSide = cSide;
			}

			/* save boundary side (bSide) into boundary side list */
			sidePtr_t *aBCside = calloc(1, sizeof(sidePtr_t));
			if (!aBCside) {
				printf("| ERROR: could not allocate aBCside\n");
				exit(1);
			}

			aBCside->side = bSide;
			aBCside->next = firstBCside;
			firstBCside = aBCside;
		}

		aSide->next = firstSide;
		firstSide = aSide;

		aSide->connection = bSide;
		bSide->connection = aSide;
	}
	free(sideList);

	/* extend element info: area and projection of cell onto axes */
	totalArea_q = 0.0;
	elem_t *aElem = firstElem;
	while (aElem) {
		createElemInfo(aElem);
		aElem = aElem->next;
	}

	/* generate virtual barycenters and ghostcells */
	sidePtr_t *aBCside = firstBCside;
	while (aBCside) {
		side_t *aSide = aBCside->side;

		double tmp = 2.0 *
			fabs(aSide->connection->GP[X] * aSide->connection->n[X] +
			     aSide->connection->GP[Y] * aSide->connection->n[Y]);

		aSide->elem->bary[X] = aSide->connection->elem->bary[X] + tmp +
			aSide->connection->n[X];
		aSide->elem->bary[Y] = aSide->connection->elem->bary[Y] + tmp +
			aSide->connection->n[Y];

		aBCside = aBCside->next;
	}

	/* extend side info: vector between barycenters, gaussian integration
	 * points and normal vectors */
	side_t *aSide = firstSide;
	while (aSide) {
		createSideInfo(aSide);
		aSide = aSide->next;
	}

	/* periodic BCs */
	connectPeriodicBC();

	/* variables for reconstruction */
	aElem = firstElem;
	while (aElem) {
		createReconstructionInfo(aElem);
		aElem = aElem->next;
	}

	/* element and side lists */
	side = calloc(nSides, sizeof(side_t *));
	if (!side) {
		printf("| ERROR: could not allocate side\n");
		exit(1);
	}

	elem = calloc(nElems, sizeof(elem_t *));
	if (!elem) {
		printf("| ERROR: could not allocate elem\n");
		exit(1);
	}

	BCside = calloc(nBCsides, sizeof(side_t *));
	if (!BCside) {
		printf("| ERROR: could not allocate BCside\n");
		exit(1);
	}

	aElem = firstElem;
	iElem = 0;
	while (aElem) {
		elem[iElem++] = aElem;
		aElem = aElem->next;
	}

	aSide = firstSide;
	long iSide = 0;
	while (aSide) {
		side[iSide++] = aSide;
		aSide = aSide->next;
	}

	aBCside = firstBCside;
	iSide = 0;
	while (aBCside) {
		if (aBCside->side->BC->BCtype != PERIODIC) {
			BCside[iSide++] = aBCside->side;
		}

		aBCside = aBCside->next;
	}
	nBCsides = iSide;
}

/**
 * \brief Read in and store the mesh parameters from the parameter file
 */
void readMesh(void)
{
	meshType = getInt("meshType", "1"); /* default is cartesian */
	switch (meshType) {
	case UNSTRUCTURED:
		printf("| Mesh Type is UNSTRUCTURED\n");

		char *tmp = getStr("meshFormat", NULL);
		strcpy(strMeshFormat, tmp);
		printf("| Mesh File Format is %s\n", strMeshFormat);
		free(tmp);

		tmp = getStr("meshFile", NULL);
		strcpy(strMeshFile, tmp);
		printf("| Mesh File: %s%s\n", strMeshFile, strMeshFormat);
		free(tmp);

		strcat(strMeshFile, strMeshFormat);
		break;
	case CARTESIAN:
		printf("| Mesh Type is CARTESIAN\n");

		cartMesh.iMax = getInt("nElemsX", NULL);
		cartMesh.jMax = getInt("nElemsY", NULL);

		double *X0 = getDblArray("X0", NDIM, NULL);
		xMin = X0[0];
		yMin = X0[1];
		free(X0);

		double *DX = getDblArray("Xmax", NDIM, NULL);
		xMax = DX[0];
		yMax = DX[1];

		free(DX);

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
				for (int i = 0; i < cartMesh.nBC[iSide]; ++i) {
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

/**
 * \brief Initialize the mesh and call `readMesh`
 */
void initMesh(void)
{
	printf("\nInitializing Mesh:\n");
	readMesh();
	strcat(strcpy(gridFile, strOutFile), "_mesh.cgns");
	createMesh();
	if ((iVisuProg == CGNS) && (!isRestart)) {
		cgnsWriteMesh();
	}
	dxRef = sqrt(1.0 / (totalArea_q * nElems));
}

/**
 * \brief Free all allocated memory of the mesh
 */
void freeMesh(void)
{
	/* free all nodes */
	node_t *aNode = firstNode;
	while (aNode) {
		if (aNode->next) {
			node_t *tmp = aNode;
			aNode = aNode->next;
			free(tmp);
		} else {
			free(aNode);
			break;
		}
	}

	/* free all elements and their corresponding sides */
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		free(aElem->xGP);
		free(aElem->wGP);
		free(aElem->node);

		side_t *aSide = aElem->firstSide;
		while (aSide) {
			if (aSide->nextElemSide) {
				side_t *tmp = aSide;
				aSide = aSide->nextElemSide;
				free(tmp);
			} else {
				free(aSide);
				break;
			}
		}

		free(aElem);
	}
	free(elem);
	free(side);

	/* free all BC sides */
	sidePtr_t *aBCside = firstBCside;
	while (aBCside) {
		if (aBCside->next) {
			sidePtr_t *tmp = aBCside;
			aBCside = aBCside->next;
			free(tmp->side->elem);
			free(tmp->side);
			free(tmp);
		} else {
			free(aBCside->side->elem);
			free(aBCside->side);
			free(aBCside);
			break;
		}
	}
}
