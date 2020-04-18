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

/**
 * \brief Secondary side lists used for various things
 */
struct sidePtr_t {
	side_t *side;			/**< pointer to a side */
	sidePtr_t *next;		/**< pointer to the next side in the
						secondary list */
};

/**
 * \brief Structure for a single element in the global element list
 */
struct elem_t {
	int elemType;			/**< element type: triangle (3) or
						quadrangle (4) */
	long id;			/**< unique element Id */
	int domain;			/**< flow domain number */
	double bary[NDIM];		/**< coordinates ob element barycenter */
	double sx;			/**< cell extension in x-direction */
	double sy;			/**< cell extension in y-direction */
	double area;			/**< area of the element */
	double areaq;			/**< inverse of element area */
	double pVar[NVAR];		/**< primitive variables of element */
	double cVar[NVAR];		/**< conservative variables of element */
	double cVarStage[NVAR];		/**< conservative variables at initial
						Runge-Kutta stage */
	double u_x[NVAR];		/**< x-gradient of primitive variables */
	double u_y[NVAR];		/**< y-gradient of primitive variables */
	double u_t[NVAR];		/**< t-gradient of primitive variables */
	double source[NVAR];		/**< source term */
	double dt;			/**< element time step */
	double dtLoc;			/**< local element time step */
	double venkEps_sq;		/**< Venkatakrishnan limiter constant
						for element */
	int innerSides;			/**< number of non-BC sides of element */
	int nGP;			/**< number of Gaussian integration points */
	double **xGP;			/**< Gaussian points for volume integral */
	double *wGP;			/**< Gaussian weights */
	side_t *firstSide;		/**< pointer to the first side of the element */
	elem_t *next;			/**< pointer to the next element in
						global element list */
	node_t **node;			/**< pointer array of the element's nodes */
};

/**
 * \brief Structure holding the information for a cartesian mesh
 */
struct cartMesh_t {
	int iMax;			/**< number of cells in x-direction */
	int jMax;			/**< number of cells in y-direction */
	int *nBC;			/**< number of different BC per side */
	int BCtype[2 * NDIM][NBC];	/**< list of BC types per side */
	int BCrange[2 * NDIM][NBC][2];	/**< list of BC ranges per side */
};

extern char parameterFile[STRLEN];
extern char strMeshFormat[STRLEN];
extern char strMeshFile[STRLEN];
extern char strIniCondFile[STRLEN];

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
extern double xMin;
extern double xMax;
extern double yMin;
extern double yMax;
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
