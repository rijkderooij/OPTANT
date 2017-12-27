#ifndef LINEARSTATICSOLVER_H
#define LINEARSTATICSOLVER_H

#include "fei_iostream.hpp" 
#include "Solver.h"

// Forward Declarations
class Domain;
class OutputRequest;

// Type definitions
typedef fei::SharedPtr<fei::Vector> feiVectorPtr;
typedef fei::SharedPtr<fei::Matrix> feiMatrixPtr;
typedef fei::SharedPtr<fei::LinearSystem> feiLinearSystemPtr;
typedef fei::SharedPtr<fei::Solver> feiSolverPtr;

class LinearStaticSolver : public Solver{
	
	public:
		// Constructor and Deconstructor
		LinearStaticSolver(Domain* domPtr, int globalLoadCaseID, const char* outputFile, bool verbose);
		virtual ~LinearStaticSolver();
		
		// Functions
		virtual int Initialize(MPI_Comm &comm);
		virtual int Prepare();
		virtual int Solve();
		virtual int WriteOutput(bool writeLoadCaseInfo=true);
		
		// Getters
		feiVectorPtr DisplacementVector(){return DisplacementVector_;}
	
	protected:
		
		feiMatrixPtr StiffnesMatrix_;
		feiVectorPtr DisplacementVector_;
		feiVectorPtr ForceVector_;
		feiLinearSystemPtr LinSys_;
		feiSolverPtr Solver_;
		
	
	private:
	
};
#endif
