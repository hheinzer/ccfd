/** \file
 *
 * \brief Contains the definitions of all structs for the mesh handling
 *
 * \author hhh
 * \date Tue 24 Mar 2020 12:45:49 PM CET
 */

#ifndef MESH_H
#define MESH_H

typedef struct node_t node_t;
typedef struct elem_t elem_t;
typedef struct side_t side_t;
typedef struct sidePtr_t sidePtr_t;
typedef struct cartMesh_t cartMesh_t;

#include "main.h"
#include "boundary.h"

/**
 * \brief Structure for a single node in a linked list of nodes
 */
struct node_t {
	long id;			/**< unique node ID */
	double x[NDIM];			/**< coordinates of the node */
	node_t *next;			/**< next node in the list */
};

/**
 * \brief Structure for a single side in the global side list
 */
struct side_t {
	long id;			/**< unique side ID */
	int BCtype;			/**< boundary condition type */
	int BCid;			/**< boundary condition Sub-ID */
	boundary_t *BC;			/**< pointer to the boundary condition */
	double pVar[NVAR];		/**< primitive variables state at side */
	double n[NDIM];			/**< normal vector of side */
	double len;			/**< length of the side */
	double baryBaryVec[NDIM];	/**< vector from element barycenter to
						barycenter of neighbor element */
	double baryBaryDist;		/**< length of `baryBaryVec` */
	double GP[NDIM];		/**< vector from element barycenter to
						the Gaussian point of the side */
	double w[NDIM];			/**< omegaX and omegaY entries for 2nd
						order gradient reconstruction */
	double flux[NVAR];		/**< numerical flux of the side */
	side_t *connection;		/**< neighbor side */
	side_t *nextElemSide;		/**< pointer to the next side of the
						element */
	side_t *next;			/**< point to the next side in the global
						side list */
	node_t *node[2];		/**< pointer array to the nodes of the side */
	elem_t *elem;			/**< pointer to the element of the side */
};

struct sidePtr_t {
	side_t *side;
	sidePtr_t *next;
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
extern sidePtr_t *firstBCside;

void initMesh(void);
void freeMesh(void);

#endif
