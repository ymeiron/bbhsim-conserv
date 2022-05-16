#ifndef DEFINITIONS_H
#define DEFINITIONS_H
/******************************************************************************
 * definitions.h *
 * Contains standard constants and control tags *
 **************************************************************************** */

/*** Simulation Name ***/
#define Name				"Initial Conditions Test"
#define ShortName			"km101226b"
#define Documentation		"The initial position of the two BHs is R0=10; the orbit loaded is sim101117g."


/*** Parameters for Kinematical Maps ***/
#define XMAX				22.0
#define RES					65				// Must be odd!
#define VBINS				100
#define INTEGRATION_TIME	10000.0
#define SNAPSHOTS {9925.2,  9933. ,  9940.7,  9948.5,  9956.2,  9963.9,  9971.7,  9979.4,  9987.2,  9994.9}
//#define SNAPSHOTS {9933.5,  9940.7,  9947.9,  9955.1,  9962.3,  9969.6,  9976.8,  9984. ,  9991.1,  9998.3}

// #define SOFT 0.004

/*** Main Simulation Parameters ***/
#define M					1.0
//#define q					1
//#define a					2 // NOTE: 'a' is the maximal separation and not semi-major axis!

//#define ORBIT_INFINITY		210.0
//#define RMAX				60.0
#define TIDAL_RADIUS		0.0001
#define CLOSENESS           0.01

/*** Flags ***/
// Possible flags are: fNOBULGE, fINNER, fOUTER, fUNIFORM_VELOCITY, fADDITIONAL_DATA, fEccentricity, fCIRCULAR
//#define fEccentricity
//#define fCIRCULAR
//#define fFAST_SEARCH
//#define fNOBULGE

/*** Accuracy Levels ***/
// be careful here!
#define H					1
#define EPS					1e-10


/*** STANDARD CONSTANTS ***/
// DO NOT TOUCH!!!

#define ORBIT_STAT_STABLE         0
#define ORBIT_STAT_CLOSETOBH1     1
#define ORBIT_STAT_CLOSETOBH2     2
#define ORBIT_STAT_VERY_FAR       3

#define ORBIT_STAT_TERMINATE    128

#define ORBIT_STAT_EXCEED_MAXIT 251
#define ORBIT_STAT_DIVERGENT    252
#define ORBIT_STAT_FALLTOBH2    253
#define ORBIT_STAT_FALLTOBH1    254
#define ORBIT_INTEGRATOR_FAIL   255


#endif
