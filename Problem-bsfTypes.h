/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-bsfTypes.h (Predefined Problem-depended BSF Types)
Prefix: PT_bsf
Author: Nikolay A. Olkhovsky
This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Order parameters
	int k;						// Number of receptive point
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	int inequalityNo;
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_float_T objectiveDistance;
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Optional filling. Do not delete!
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Optional filling. Do not delete!
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Optional filling. Do not delete!
};