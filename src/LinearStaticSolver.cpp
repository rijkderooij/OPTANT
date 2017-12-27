#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 


// Constructor
LinearStaticSolver::LinearStaticSolver(Domain* domPtr, int globalLoadCaseID,const char* outputFile,  bool verbose)
							: Solver(domPtr,globalLoadCaseID,outputFile,verbose) {
}

// Destructor
LinearStaticSolver::~LinearStaticSolver(){
	Solver_.reset();
}

// Initialize LinearStatic system
int LinearStaticSolver::Initialize(MPI_Comm &comm){
	Solver::Initialize(comm);
}

// Prepare LinearStatic system
int LinearStaticSolver::Prepare(){
	
	// ======================================================== //
	// 		C O M P U T E  E L E M E N T  M A T R I C E S		//
	// -------------------------------------------------------- //
	
	// Compute element matrices only if this has not been done yet
	if(!DomainPtr_->ElementMatricesBuilt){
	
		if(verbose_) FEI_COUT<< FEI_ENDL << "BUILDING ELEMENT MATRICES..." << FEI_ENDL;
		
		// Initialize timer for matrix creation
		double fei_creatematrix_start_time = fei::utils::cpu_time();
		
		// Build all element stiffness matrices
		if(verbose_) CHK_ERR(DomainPtr_->BuildElementMatrices());
		
		// Compute time for creating the matrix and print
		double fei_create_element_matrices_time = fei::utils::cpu_time() - fei_creatematrix_start_time;
		if (verbose_) {
			FEI_COUT << "DONE! " << FEI_ENDL;
			FEI_COUT << "CPU time to build element matrices:   " << fei_create_element_matrices_time << FEI_ENDL;
		}
	}
	
	// ======================================================== //
	// 			S E T  U P  L I N E A R  S Y S T E M			//
	// -------------------------------------------------------- //
	
	if(verbose_) FEI_COUT<< FEI_ENDL << "SETTING UP LINEAR SYSTEM..." << FEI_ENDL;
			
	// Initialize timer to load the matrix
	double start_load_time = fei::utils::cpu_time();
	
	// Create StiffnesMatrix, DisplacementVector, and ForceVector based on matrixGraph
	StiffnesMatrix_ 	= Factory_->createMatrix(MatrixGraph_);
	DisplacementVector_ = Factory_->createVector(MatrixGraph_, true);
	ForceVector_  		= Factory_->createVector(MatrixGraph_);
	LinSys_ 			= Factory_->createLinearSystem(MatrixGraph_);
	
	// Apply solverParams to linear system	
	CHK_ERR( LinSys_->parameters(DomainPtr_->SolverParams));
	
	// Set StiffnesMatrix, DisplacementVector, and ForceVector to LinSys
	LinSys_->setMatrix(StiffnesMatrix_);
	LinSys_->setSolutionVector(DisplacementVector_);
	LinSys_->setRHS(ForceVector_);
	
	if(verbose_) DomainPtr_->AssembleStiffnessMatrix(MatrixGraph_.get(),StiffnesMatrix_.get());
	
	// Apply constraints and loads for this load case: id=GlobalLoadCaseID_
	if(verbose_){
		CHK_ERR(DomainPtr_->ApplyConstraints(LinSys_.get(),GlobalLoadCaseID_));
		CHK_ERR(DomainPtr_->ApplyLoads(MatrixGraph_.get(),ForceVector_.get(),GlobalLoadCaseID_));
	}
	// Loading of LinSys is completed
	CHK_ERR( LinSys_->loadComplete() );
	
	// Compute loading time
	double fei_load_time = fei::utils::cpu_time() - start_load_time;
	
	// Print output
	if (verbose_) {
		FEI_COUT << "DONE! " << FEI_ENDL;
		FEI_COUT << "CPU time to load element data into linear system:    " << fei_load_time << FEI_ENDL;
	}
	
	return 0;
}

// Solve LinearStatic system
int LinearStaticSolver::Solve(){
	
	// ======================================================== //
	// 			S O L V E  L I N E A R  S Y S T E M  			//
	// -------------------------------------------------------- //
	
	if (verbose_) FEI_COUT << FEI_ENDL<< "SOLVING LINEAR SYSTEM..." << FEI_ENDL;
	
	// First, check whether the solverParams specify a "Trilinos_Solver" 
	// string (e.g. MUMPS). If not, it won't matter but if so, the 
	// factory may switch based on this value, when creating the solver 
	// object instance.
	
	std::string solver_name_value;
	DomainPtr_->SolverParams.getStringParamValue("Trilinos_Solver", solver_name_value);

	// Create pointer to solver_name_value, depending on whether it is empty
	const char* charptr_solvername =
	solver_name_value.empty() ? NULL : solver_name_value.c_str();

	// Create solver
	Solver_ = Factory_->createSolver(charptr_solvername);
	
	// Start solution procedure
	int status;
	int itersTaken = 0;

	double start_solve_time = fei::utils::cpu_time();

	int err = Solver_->solve(LinSys_.get(),
			  NULL, 			//preconditioningMatrix
			  DomainPtr_->SolverParams,
			  itersTaken,		//output: # of iterations
			  status); 			//output: 0 -> successful

	// Compute solver time
	double solve_time = fei::utils::cpu_time()-start_solve_time;
	
	// Check for error and print
	if (err!=0 && verbose_) {
		FEI_COUT << "solve returned err: " << err <<", status: "
			   << status << FEI_ENDL;
	}

	// Print solver time
	if (verbose_) {
		FEI_COUT << "DONE! " << FEI_ENDL;
		FEI_COUT << "CPU time for solving: " << solve_time << FEI_ENDL;
	}

	// Scatter data from the underlying non-overlapping data decomposition 
	// to the overlapping data decomposition. In other words, update 
	// values for shared indices from underlying uniquely owned data.
	CHK_ERR(DisplacementVector_->scatterToOverlap() );
	
	
	return 0;
}

// WriteOutput for LinearStatic system
int LinearStaticSolver::WriteOutput(bool writeLoadCaseInfo){
	
	// Call writeoutput of solver
	if(writeLoadCaseInfo){
		Solver::WriteOutput();
	}
	
	// Write output for this loadcase
	if(verbose_){
		// Write outputRequests  
		CHK_ERR(DomainPtr_->WriteOutput(this,OutputFile_,GlobalLoadCaseID_));
	}
	
	
	
	return 0;
}
