/*
 * mesh.h
 *
 * Created: Tue 24 Mar 2020 12:45:49 PM CET
 * Author : hhh
 */

#ifndef MESH_H
#define MESH_H

#include "main.h"
#include "boundary.h"

typedef struct node_t node_t;
struct node_t {
	long id;
	node_t *next;
	double x[NDIM];
};

typedef struct elem_t elem_t;
typedef struct side_t side_t;
struct side_t {
	long id;
	int BCtype;
	int BCid;
	boundary_t *BC;
	double pVar[NVAR];
	double n[NDIM];
	double len;
	double baryBaryVec[NDIM];
	double baryBaryDist;
	double GP[NDIM];
	double w[NDIM];
	double flux[NVAR];
	side_t *connection;
	side_t *nextElemSide;
	side_t *next;
	node_t *node[2];
	elem_t *elem;
};

struct elem_t {
	int elemType;
	long id;
	int domain;
	double bary[NDIM];
	double sx;
	double sy;
	double area;
	double areaq;
	double pVar[NVAR];
	double cVar[NVAR];
	double cVarStage[NVAR];
	double u_x[NVAR];
	double u_y[NVAR];
	double u_t[NVAR];
	double source[NVAR];
	double dt;
	double dtLoc;
	double venkEps_sq;
	int innerSides;
	int nGP;
	double **xGP;
	double *wGP;
	side_t *firstSide;
	elem_t *next;
	node_t **node;
};

typedef struct cartMesh_t cartMesh_t;
struct cartMesh_t {
	int iMax, jMax;
	int *nBC;
	int BCtype[2 * NDIM][NBC];
	int BCrange[2 * NDIM][NBC][2];
};

extern char parameterFile[STRLEN],
	    strMeshFormat[STRLEN],
	    strMeshFile[STRLEN],
	    strIniCondFile[STRLEN];

extern cartMesh_t cartMesh;
extern char gridFile[STRLEN];

extern int meshType;
extern int meshFormat;

extern long nNodes;

extern long nElems;
extern long nTrias;
extern long nQuads;

extern long nSides;
extern long nBCsides;
extern long nInnerSides;

extern double totalArea_q;
extern double xMin, xMax;
extern double yMin, yMax;
extern double dxRef;

extern elem_t **elem;
extern side_t **side;
extern side_t **BCside;

extern elem_t *firstElem;
extern node_t *firstNode;
extern side_t *firstSide;
extern side_t *firstBCside;

void initMesh(void);

#endif
