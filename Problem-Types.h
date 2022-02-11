/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-Types.h (Problem Types)
Prefix: PT
Author: Nikolay A. Olkhovsky
This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef float						PT_float_T;
typedef unsigned long long int		PT_integer_T;
typedef PT_float_T					PT_vector_T[PP_MAX_N];
typedef PT_vector_T					PT_matrix_T[PP_MAX_M];
//typedef PT_float_T				PT_point_T[PP_MAX_N];
typedef PT_float_T					PT_column_T[PP_MAX_M];
#ifdef PP_FILES_OUT
//typedef PT_float_T					PT_image_T[PP_MAX_K][PP_MAX_N + 1];
typedef PT_float_T**				PT_image_T;
#endif