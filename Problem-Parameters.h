/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky

This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/

//=========================== Problem Parameters =========================
#define PP_N		3			// Number of dimensions
#define PP_MAX_M	50					// Maximal  number of inequalities
#define PP_MAX_N	50					// Maximal  space dimension
#define PP_ETA		14					// Rank of receptive field 		
#define PP_DELTA	1					// Density of receptive field
#define PP_MAX_K	100000				// Maximal number of receptive field points
#define PP_PATH		"C:/Users/Akella/source/repos/nikolay-olkhovskiy/LP-visualization-MPI/"	// Working directory of the application
#define PP_LPP_FILE	"input_example.txt"	// File with initial data
#define PP_OUT_FILE "image.txt"			// File with output results
#define PP_INPUT_FILE	"point.txt"		// File with coordinates of target point
#define PP_PICTURE_FILE	"picture.txt"	//File with coordinates suitable to display
#define PP_PICTURE_OUT					// Flag to create additional file with image points coordinates

//-------------------------- Macroses ---------------------------
#define PF_MIN(x,y) (x<y?x:y)