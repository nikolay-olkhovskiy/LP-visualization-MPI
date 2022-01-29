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
typedef float								PT_float_T;
typedef int									PT_integer_T;
typedef std::valarray<PT_float_T>			PT_vector_T;
typedef std::valarray<PT_vector_T>			PT_matrix_T;
typedef std::valarray<PT_float_T>			PT_point_T;
typedef std::valarray<PT_float_T>			PT_column_T;
typedef std::pair<PT_point_T, PT_float_T>	PT_pair_T;
typedef std::vector<PT_pair_T>				PT_image_T;