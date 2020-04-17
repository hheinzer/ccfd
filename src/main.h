/** \file
 *
 * Contains the global constants and definitions
 *
 * \author hhh
 * \date Sat 21 Mar 2020 10:44:50 AM CET
 */

#ifndef MAIN_H
#define MAIN_H

#define STRLEN 256		/**< string length */

/**
 * Index aliases for the conservative variables vector
 */
enum conservativeVariables {
	RHO,			/**< density */
	MX,			/**< momentum in x-direction */
	MY,			/**< momentum in y-direction */
	E			/**< energy */
};

/**
 * Index aliases for the primitive variables vector
 */
enum primitiveVariables {
	VX = 1,			/**< velocity in x-direction */
	VY,			/**< velocity in y-direction */
	P			/**< pressure */
};

/**
 * Index aliases for the x- and y-components of a vector
 */
enum directions {
	X,			/**< x-direction */
	Y			/**< y-direction */
};

/**
 * Aliases for the different boundary condition types
 */
enum boundaryConditionType {
	SLIPWALL = 1,		/**< slip wall, or Euler wall */
	WALL,			/**< no slip wall, or Navier-Stokes wall */
	INFLOW,			/**< supersonic inflow */
	OUTFLOW,		/**< supersonic outflow */
	CHARACTERISTIC,		/**< characteristic BC, or subsonic inflow */
	EXACTSOL,		/**< exact solution */
	PERIODIC,		/**< periodic BC */
	PRESSURE_OUT		/**< subsonic pressure outflow */
};

/**
 * Aliases for the sides of a cartesian mesh
 */
enum cartesianMeshSides {
	BOTTOM,			/**< bottom side */
	RIGHT,			/**< right side */
	TOP,			/**< top side */
	LEFT			/**< left side */
};

/**
 * Flag for the mesh type
 */
enum meshType {
	UNSTRUCTURED,		/**< unstructured */
	CARTESIAN		/**< cartesian */
};

/**
 * Aliases for the different flux functions
 */
enum fluxFunction {
	GOD = 1,		/**< Godunov flux function */
	ROE,			/**< Roe flux function */
	HLL,			/**< HLL flux function */
	HLLE,			/**< HLLE flux function */
	HLLC,			/**< HLLC flux function */
	LXF,			/**< Lax-Friedrichs flux function */
	STW,			/**< Steger-Warming flux function */
	CEN,			/**< central flux function */
	AUSMD,			/**< AUSMD flux function */
	AUSMDV,			/**< AUSMDV flux function */
	VANLEER			/**< Van Leer flux function */
};

/**
 * Flag for the limiter function
 */
enum limiterFunction {
	BARTHJESPERSEN = 1,	/**< Barth & Jespersen limiter */
	VENKATAKRISHNAN		/**< Venkatakrishnan limiter */
};

/**
 * General parameters for the Program
 */
enum generalParameters {
	NDIM = 2,		/**< number of dimensions, cannot be changed */
	NVAR = 4,		/**< number of variables, primitive or conservative */
	NBC = 20		/**< maximum number of boundary conditions */
};

/**
 * Output format for the results
 */
enum ioFormat {
	CGNS = 1,		/**< .CGNS file format */
	CURVE,			/**< .curve file format */
	CSV			/**< .csv file format */
};

/**
 * Index aliases for the residual vector
 */
enum clcdResiduals {
	CL = 4,			/**< lift coefficient */
	CD			/**< drag coefficient */
};

#endif
