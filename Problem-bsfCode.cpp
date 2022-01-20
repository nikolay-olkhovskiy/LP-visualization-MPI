/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: Problem-bsfCode.cpp (Problem-dependent Code)
Prefix: PC
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

//----------------------- Predefined problem-dependent functions -----------------
void PC_bsf_Init(bool* success) {

}

void PC_bsf_SetListSize(int* listSize) {

}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {

}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {	// For Job 0

}

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {// For Job 1
	// Optional filling. Do not delete!
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {// For Job 2
	// Optional filling. Do not delete!
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {// For Job 3
	// Optional filling. Do not delete!
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) {			// For Job 0
	// z = x (+) y
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {	// For Job 1
	// Optional filling. Do not delete!
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {	// For Job 2
	// Optional filling. Do not delete!
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {	// For Job 3
	// optional filling
}

void PC_bsf_ProcessResults(		// For Job 0
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
) {

}

void PC_bsf_ProcessResults_1(	// For Job 1	
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
) {
	// Optional filling. Do not delete!
}

void PC_bsf_ProcessResults_2(	// For Job 2
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
	) {
	// Optional filling. Do not delete!
}

void PC_bsf_ProcessResults_3(	// For Job 3
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
	) {
	// Optional filling. Do not delete!
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit
) {
	// Optional filling. Do not delete!
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Problem parameters ====================================================" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP

}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 0
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;


}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 1
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Optional filling. Do not delete!

}

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 2
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Optional filling. Do not delete!

}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 3
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Optional filling. Do not delete!

}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 0

}

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 1
	// Optional filling. Do not delete!
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 2
	// Optional filling. Do not delete!
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 3
	// Optional filling. Do not delete!
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {

}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {

}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//----------------------------- User functions -----------------------------
