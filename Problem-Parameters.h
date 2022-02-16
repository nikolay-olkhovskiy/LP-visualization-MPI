/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky

This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/

//=========================== Problem Parameters =========================
#define PP_N		3									// Number of dimensions
#define PP_MAX_M	250									// Maximal  number of inequalities
#define PP_MAX_N	120									// Maximal  space dimension
#define PP_ETA		2									// Rank of receptive field 		
#define PP_DELTA	1									// Density of receptive field
#define PP_MAX_K	18446744073709550592ull				// Maximal number of receptive field points
//#define PP_PATH		"C:/Users/Akella/source/repos/nikolay-olkhovskiy/LP-visualization/"	// Working directory of the application
#define PP_PATH		"/home/olkhovskiina/"	// Working directory of the application
#define PP_LPP_FILE	"input_example.txt"	// File with initial data
//#define PP_LPP_FILE	"python/lpp7_100.txt"
#define PP_IMAGE_FILE "image.txt"		// File with output results
#define PP_INPUT_FILE	"point.txt"		// File with coordinates of target point
#define PP_PICTURE_FILE	"picture.txt"	// File with coordinates suitable to display
#define PP_RECEPTIVE_FIELD_FILE "receptive_field.txt" // File with coordinates of receptive field points
//#define PP_IMAGE_OUT					// Flag to save list of target distances to the file
//#define PP_RECEPTIVE_FIELD_OUT			// Flag to save receptive field points to the file

//-------------------------- Macroses ---------------------------
#define PF_MIN(x,y) (x<y?x:y)