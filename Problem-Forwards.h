/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-bsf-Forwards.h (Problem Function Forwards)
Author: Nikolay A. Olkhovsky 
This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-bsfTypes.h"	// Predefined BSF types
#include "Problem-Types.h"		// Problem Types
//====================== Problem Functions ===========================
inline PT_point_T floatsToValarray(PT_float_T arr[]);
inline void basis_Init();
inline void basis_Print();
inline void G(PT_bsf_parameter_T* parameter);
inline bool parameterOutOfRetina(PT_bsf_parameter_T* parameter);

// Helper functions for MapF implementation
inline bool isInnerPoint(PT_point_T point);
inline PT_point_T targetProjection(int i, PT_point_T x);
inline PT_float_T targetDistance(PT_point_T x);