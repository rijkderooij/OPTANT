#ifndef LINEARBUCKLINGSOLVER_H
#define LINEARBUCKLINGSOLVER_H

#include "fei_iostream.hpp" 
#include "LinearStaticSolver.h"

// Forward Declarations
class Domain;
class OutputRequest;
class ArpackOperations;

template <class T1,class T2>
class ARSymStdEig;

template <class T1,class T2, class T3>
class ARSymGenEig;

class LinearBucklingSolver : public LinearStaticSolver{
	
	public:
		// Constructor
		LinearBucklingSolver(Domain* domPtr, int globalLoadCaseID, int gPSLoadCaseID, const char* outputFile,  bool verbose);
		~LinearBucklingSolver();
		
		// Functions
		int Initialize(MPI_Comm &comm);
		int Prepare();
		int Solve();
		int WriteOutput(bool writeLoadCaseInfo=true);
		
		// Getters
		ARSymStdEig<double,ArpackOperations>* ArpackProblem(){return ProbSymStdEigB_;};
		int NumConvergedModes(){return nconvStdB_;}; 
	
	protected:
	
	private:
		// Create random initial vector for buckling analysis
		void InitialVector(int n, double* initVec);
		
		int nModes_;									// Number of desired modes
	
		// Geometric matrices/vectors
		feiMatrixPtr GeomStiffnesMatrix_;
		feiVectorPtr GeomDisplacementVector_;
		feiVectorPtr GeomForceVector_;
		feiLinearSystemPtr GeomLinSys_;
		
		// Prestress variables
		Solver* psSolver_;
		
		
		
		// Arpack problems
		ARSymStdEig<double,ArpackOperations> *ProbSymStdEigB_;
		int nconvStdB_;									// Number of converged values of ProbSymStdEigB_
		
		
		ARSymGenEig<double , ArpackOperations,ArpackOperations> *ProbSymGenEigAB_;
		int nconvGenAB_;								// Number of converged values of ProbSymGenEigAB_
		
		
};
#endif
