/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Nikolay A. Olkhovsky
This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Types.h"			// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int		PD_n;		// Space dimension
static int		PD_m;		// Number of inequalities
static int		PD_K;		// Number of receptive fied points // (must be < 2 147 483 647)	
static int		PD_retina_k = 0;	// Currend number of retina point
static int		PD_recept_k = 0;	// Currend number of point of receptive field
//========================== Problem data structures ==============================
static PT_matrix_T PD_A;	// Matrix of the system Ax <= b
static PT_column_T PD_b;	// Column of the constant terms of the system Ax <= b
static PT_vector_T PD_c;	// Coefficients of the objective function <c,x>
static PT_matrix_T PD_E;	// Matrix of vectors e(i) forming basis othogonal to objective function
static PT_point_T PD_g;		// Point of retina
static PT_point_T PD_z = {0., 0., 200.};	// Center of retina

static PT_image_T PD_I; // Retina
//========================== Files ================================================
static string PD_lppFile; /* LPP file in the following format:
------------ begin of file -------------
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
c_1 c_2 ... c_n
------------ end of file----------------*/
static string PD_outFile; /* OUT file in the following format:
------------ begin of file -------------
m
I_1
I_2
...
I_m
------------ end of file----------------*/