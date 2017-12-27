#ifndef SOLVER_H
#define SOLVER_H

#include "fei_iostream.hpp" 

// Forward Declarations
class Domain;
class OutputRequest;

// Type definitions
typedef fei::SharedPtr<fei::Factory> feiFactoryPtr;
typedef fei::SharedPtr<fei::MatrixGraph> feiMatrixGraphPtr;
typedef fei::SharedPtr<fei::VectorSpace> feiVectorSpacePtr;
typedef fei::SharedPtr<fei::Vector> feiVectorPtr;

class Solver{
	
	public:
		// Constructor and destructor
		Solver(Domain* domPtr, int globalLoadCaseID, const char* outputFile, bool verbose);
		virtual ~Solver();
		
		// Functions
		virtual int Initialize(MPI_Comm &comm);
		virtual int Prepare()=0;
		virtual int Solve() = 0;
		virtual int WriteOutput(bool writeLoadCaseInfo=true);
		
		// Getters
		feiVectorSpacePtr NodeSpace(){return NodeSpace_;}
		feiMatrixGraphPtr MatrixGraph(){return MatrixGraph_;}
		virtual feiVectorPtr DisplacementVector()=0;
	
	protected:
		Domain* DomainPtr_;
		int GlobalLoadCaseID_;
		const char* OutputFile_;
		bool verbose_;		  // True if this is processor 1, false otherwise
		
		
		feiFactoryPtr Factory_;
		feiVectorSpacePtr NodeSpace_;
		feiMatrixGraphPtr MatrixGraph_;
		
	
	private:
		
	
};
#endif
