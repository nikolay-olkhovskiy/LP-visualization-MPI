/*==============================================================================
Project: Linear programming by deep neural network
Theme: Visualization of Linear Programming Problem (ViLiPP)
Module: Problem-bsfCode.cpp (Problem-dependent Code)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code is developed based on the BSF skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

//----------------------- Predefined problem-dependent functions -----------------
void PC_bsf_Init(bool* success) {
	FILE* stream;
	PT_float_T buf;
	const char* lppFile;
	// ------------- Load LPP data -------------------

	PD_lppFile = PP_PATH;
	PD_lppFile += PP_LPP_FILE;
	lppFile = PD_lppFile.c_str();
	stream = fopen(lppFile, "r");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << lppFile << "'.\n";
		*success = false; 
//		system("pause");
		return;
	}

	if (fscanf(stream, "%d%d", &PD_m, &PD_n) == 0) { 
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster) 
			cout << "Unexpected end of file" << endl; 
		*success = false; 
//		system("pause");
		return; 
	}

	PD_A.resize(PD_m);
	PD_b.resize(PD_m);
	for (int i = 0; i < PD_m; i++) {
		PD_A[i].resize(PD_n);
		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%f", &buf) == 0) { 
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster) 
					cout << "Unexpected end of file" << endl; 
				*success = false; 
//				system("pause");
				return; 
			}
			PD_A[i][j] = buf;
		}
		if (fscanf(stream, "%f", &buf) == 0) { 
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) 
				cout << "Unexpected end of file" << endl; 
			*success = false; 
//			system("pause");
			return; 
		}
		PD_b[i] = buf;
	}

	PD_c.resize(PD_n);
	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%f", &buf) == 0) { 
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) 
				cout << "Unexpected end of file" << endl; 
			*success = false; 
//			system("pause");
			return; 
		}
		PD_c[j] = buf;
	}
	fclose(stream);

	// ------------- Load target point coordinates -------------------

	PD_lppFile = PP_PATH;
	PD_lppFile += PP_INPUT_FILE;
	lppFile = PD_lppFile.c_str();
	stream = fopen(lppFile, "r");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << lppFile << "'.\n";
		*success = false;
		//		system("pause");
		return;
	}

	PD_z.resize(PD_n);
	for (int i = 0; i < PD_n; i++) {
		if (fscanf(stream, "%f", &buf) == 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file" << endl;
			*success = false;
			//				system("pause");
			return;
		}
		PD_z[i] = buf;
	}

	fclose(stream);

	basis_Init();

	PD_K = (int)powf(2 * PP_ETA + 1, (PT_float_T)PD_n - 1);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = (int)PD_m;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	parameterOutP->pointNo = parameterIn.pointNo;
	for (int i = 0; i < PD_n; i++)
		parameterOutP->receptivePoint[i] = parameterIn.receptivePoint[i];
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {	// For Job 0
	PT_point_T g = floatsToValarray(BSF_sv_parameter.receptivePoint);
	int i = mapElem->inequalityNo;
	if((PD_A[i] * PD_c).sum() > 0 && isInnerPoint(g))
		reduceElem->objectiveDistance = targetDistance(targetProjection(i, g));
	else
		reduceElem->objectiveDistance = FLT_MAX;
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
	if (isfinite(x->objectiveDistance) && isfinite(y->objectiveDistance))
		z->objectiveDistance = PF_MIN(x->objectiveDistance, y->objectiveDistance);
	else if (isfinite(x->objectiveDistance))
		z->objectiveDistance = x->objectiveDistance;
	else if (isfinite(y->objectiveDistance))
		z->objectiveDistance = y->objectiveDistance;
	else
		z->objectiveDistance = FLT_MAX;
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
	PT_vector_T g = floatsToValarray(parameter->receptivePoint);
	PD_I.push_back(make_pair(g, reduceResult->objectiveDistance));
	do {
		parameter->pointNo += 1;
		G(parameter);
	} while (parameter->pointNo < PD_K && parameterOutOfRetina(parameter));

	*exit = (parameter->pointNo >= PD_K);
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
	cout << "Dimensions: " << PD_n << endl;
	cout << "Number of inequalities: " << PD_m << endl;
	cout << "Receptive field rank: " << PP_ETA << endl;
	cout << "Receptive field density: " << PP_DELTA << endl;
	cout << "Maximum number of points: " << PD_K << endl;
	cout << "Receptive field coordinates: ";
	for (int i = 0; i < PD_n; i++) {
		cout << PD_z[i] << " ";
	}
	cout << endl;
	basis_Print();
//	system("pause");
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 0
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Point number:\t" << parameter.pointNo << endl;
	cout << "Point coordinates:\t";
	copy(parameter.receptivePoint, parameter.receptivePoint + PD_n, ostream_iterator<PT_float_T>(cout, " "));
	cout << endl;
	cout << "Z coordinates:\t";
	print_Point(PD_z);
	cout << endl;
	cout << "Field distance:\t" << sqrt(pow(floatsToValarray(parameter.receptivePoint) - PD_z, 2.0f).sum()) << endl;
	cout << "Receptive field rank:\t" << PP_ETA * PP_DELTA << endl;
//	system("pause");
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
	FILE* stream;
	const char* fileName;
	int m = (int)PD_I.size();
	int n = PD_n;

	//--------------- Output results -----------------//
	PD_outFile = PP_PATH;
	PD_outFile += PP_OUT_FILE;
	fileName = PD_outFile.c_str();
	cout << "-----------------------------------" << endl;
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d\n", m);

	for (int i = 0; i < m; i++) {
		fprintf(stream, "%.4f\n", PD_I[i].second);
	}
	fclose(stream);
	cout << "Image is saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;
#ifdef PP_PICTURE_OUT
	//-------------- Output Coordinates -------------//
	PD_outFile = PP_PATH;
	PD_outFile += PP_PICTURE_FILE;
	fileName = PD_outFile.c_str();
	
	cout << "-----------------------------------" << endl;
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d\t%d\n", m, n);

	for (int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++)
			fprintf(stream, "%.4f\t", PD_I[i].first[j]);
		fprintf(stream, "%.4f\n", PD_I[i].second);
	}
	fclose(stream);
	cout << "Coordinates are saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;
#endif
//	system("pause");
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
	parameter->pointNo = 0;
	G(parameter);
	while (parameterOutOfRetina(parameter)) {
		parameter->pointNo += 1;
		G(parameter);
	}
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->inequalityNo = i;
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
inline PT_point_T floatsToValarray(PT_float_T arr[]) {
	PT_vector_T result(PD_n);
	for (int i = 0; i < PD_n; i++)
		result[i] = arr[i];
	return result;
}
inline void basis_Init() {
	//PD_c
	int j;
	PT_float_T length;
	PT_float_T tailSum;
	PT_vector_T PD_c2 = pow(PD_c, 2.0f);
	PD_E.resize(PD_n);
	PD_E[0] = PD_c;
	for (int i = 1; i < PD_n; i++) {
		PD_E[i].resize(PD_n);
		for(j = 0; j < i; j++)	PD_E[i][j] = 0;
		tailSum = vector_Sum(PD_c2, i);
		if (tailSum == 0) {
			PD_E[i][i - 1] = 0;
			PD_E[i][i] = 1;
			j++;
			for (; j < PD_n; j++) { PD_E[i][j] = 0; }
		}
		else if (PD_c[i - 1] == 0.) {
			PD_E[i][i - 1] = 1.;
			for (; j < PD_n; j++) { PD_E[i][j] = 0; }
		}
		else {
			PD_E[i][i - 1] = (PT_float_T)((-1. * tailSum) / PD_c[i - 1]); //Possible division by zero!
			for (; j < PD_n; j++) { PD_E[i][j] = PD_c[j]; }
		}
		length = sqrt(pow(PD_E[i], 2.0f).sum());
		PD_E[i] /= length;
	}
}
inline void print_Point(PT_point_T x) {
	int N = (int)x.size();
	for (int i = 0; i < N; i++)
		cout << x[i] << " ";
}
inline void print_Vector(PT_vector_T x) {
	int N = (int)x.size();
	for (int i = 0; i < N; i++)
		cout << x[i] << " ";
}
inline void basis_Print() {
	for (int i = 0; i < (int)PD_E.size(); i++) {
		print_Vector(PD_E[i]);
		cout << endl;
	}
}
inline PT_float_T vector_Sum(PT_vector_T v, int start) {
	PT_float_T result = 0.0f;
	int N = (int)v.size();
	for (int i = start; i < N; i++) {
		result += v[i];
	}
	return result;
}
inline void G(PT_bsf_parameter_T *parameter) {
	PT_point_T tempPoint;
	PT_integer_T dimensionPointsNumber;
	vector<PT_integer_T> i;
	int pointNo = parameter->pointNo;

	i.resize(PD_n - 1);
	for (int j = PD_n - 1; j > 0; j--) {
		dimensionPointsNumber = (PT_integer_T)powf(2 * PP_ETA + 1, (PT_float_T)j - 1); //Possible overfilling!
		i[j - 1] = pointNo / dimensionPointsNumber;
		pointNo = pointNo % dimensionPointsNumber;
	}
	tempPoint = PD_z;
	for (int j = 1; j < PD_n; j++) {
		tempPoint += PD_E[j] * (float)(i[j - 1] * PP_DELTA - PP_ETA * PP_DELTA);
	}
	for(int i = 0; i < PD_n; i++)
		parameter->receptivePoint[i] = tempPoint[i];
};
inline bool parameterOutOfRetina(PT_bsf_parameter_T* parameter) {
	PT_float_T distanceToZ = sqrt(pow(floatsToValarray(parameter->receptivePoint) - PD_z, 2.0f).sum());
	return distanceToZ > PP_ETA * PP_DELTA;
}

inline bool isInnerPoint(PT_point_T point) {
	bool result = true;
	for (int i = 0; i < PD_m; i++)
		if ((PD_A[i] * point).sum() > PD_b[i])
			result = false;
	return result;
}

// Projection of receptive field point to recessive subspace gamma_i(x)
inline PT_point_T targetProjection(int i, PT_point_T x) {
	return x - (((PD_A[i] * x).sum() - PD_b[i]) / (PD_A[i] * PD_c).sum()) * PD_c;
}

// Distance from projection to retina rho_c(x)
inline PT_float_T targetDistance(PT_point_T x) {
	return (PD_c * (PD_z - x)).sum() / sqrt(pow(PD_c, 2.0f).sum());
}