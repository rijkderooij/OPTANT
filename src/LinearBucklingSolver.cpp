#include "../include/PCH_OPTANT.h"
#include <random>
#include "fei_iostream.hpp" 

#include "arssym.h"
#include "argsym.h"


// Constructors
LinearBucklingSolver::LinearBucklingSolver(Domain* domPtr,int globalLoadCaseID, int gPSLoadCaseID, const char* outputFile,  bool verbose)
						: LinearStaticSolver(domPtr,globalLoadCaseID,outputFile,verbose){
	
	// Create linear solver for prestresses
	if(gPSLoadCaseID==0)
		psSolver_ = NULL;			// No prestress
	else
		psSolver_ = new LinearStaticSolver(domPtr,gPSLoadCaseID,outputFile,verbose);
		
	if(verbose_){						
		// Extract nModes from SolverParams in Domain
		DomainPtr_->SolverParams.getIntParamValue("nModes", nModes_);
	}	
}

// Destructor
LinearBucklingSolver::~LinearBucklingSolver(){
	if(psSolver_!=NULL) delete psSolver_;
}

// Initialize LinearBuckling system
int LinearBucklingSolver::Initialize(MPI_Comm &comm){
	
	// Initialize, prepare and solve PreStress Linear System
	if(psSolver_!=NULL){
		if(verbose_) FEI_COUT << "Initialize Prestress Linear System!" << FEI_ENDL;
		psSolver_->Initialize(comm);
	}
	
	Solver::Initialize(comm);
}

// Prepare LinearBuckling system
int LinearBucklingSolver::Prepare(){
	// Write general output for the solve
	Solver::WriteOutput();
	
	// Initialize, prepare and solve PreStress Linear System
	if(psSolver_!=NULL){
		if(verbose_)
			FEI_COUT << "- LINEAR BUCKLING SOLVE - PREPARE PRESTRESS -" << FEI_ENDL;
		psSolver_->Prepare();
		psSolver_->Solve();
		psSolver_->WriteOutput(false);  // Write output, but without load case info
		
		
		// ================================================================ //
		// 			  COMPUTE PRESTRESS GEOMETRIC ELEMENT MATRICES			//
		// ---------------------------------------------------------------- //
		
		// Compute geometric element matrices
		if(verbose_) FEI_COUT<< FEI_ENDL << "BUILDING PRESTRESS GEOMETRIC ELEMENT MATRICES..." << FEI_ENDL;
		
		// Initialize timer for matrix creation
		double fei_creatematrix_start_time = fei::utils::cpu_time();
		
		// Build all element stiffness matrices. Use makeNegative=true, 
		// since this differential stiffness will be subtracted from the 
		// standard stiffness matrix
		if(verbose_) CHK_ERR(DomainPtr_->BuildGeometricElementMatrices(psSolver_,true));
		
		// Compute time for creating the matrix and print
		double fei_create_element_matrices_time = fei::utils::cpu_time() - fei_creatematrix_start_time;
		if (verbose_) {
			FEI_COUT << "DONE! " << FEI_ENDL;
			FEI_COUT << "CPU time to build prestress geometric element matrices:   " << fei_create_element_matrices_time << FEI_ENDL;
		}
		
		// Reset all element temperatures
		DomainPtr_->ResetElementTemperatures();
			
	}
	
	// Prepare and solve Linear system
	LinearStaticSolver::Prepare();
	LinearStaticSolver::Solve();
	
	
	
	if(psSolver_!=NULL){
		
		// ============================================================ //
		// 	  MODIFY LINEAR STATIC STIFFNESS MATRIX WITH PRESTRESSES	//
		// ------------------------------------------------------------ //
		
		// Assemble (negative of) geometric element stiffness matrices
		// into the global stiffness matrix. This results in the effec-
		// tive stiffness matrix including prestresses
		if(verbose_) DomainPtr_->AssembleStiffnessMatrix(MatrixGraph_.get(),StiffnesMatrix_.get(),true);
		
		// Apply constraints for this load case: id=GlobalLoadCaseID_
		if(verbose_){
			CHK_ERR(DomainPtr_->ApplyConstraints(LinSys_.get(),GlobalLoadCaseID_));
		}
		
		// Complete the loading  of LinSys_ again to apply the constraints
		CHK_ERR( LinSys_->loadComplete() );
		
	
	}
	
	if(verbose_)
		FEI_COUT << "------ LINEAR BUCKLING SOLVE - PREPARE ------" << FEI_ENDL;
	
	// ================================================================ //
	// C O M P U T E  G E O M E T R I C  E L E M E N T  M A T R I C E S	//
	// ---------------------------------------------------------------- //
	
	// Compute geometric element matrices
	if(verbose_) FEI_COUT<< FEI_ENDL << "BUILDING GEOMETRIC ELEMENT MATRICES..." << FEI_ENDL;
	
	// Initialize timer for matrix creation
	double fei_creatematrix_start_time = fei::utils::cpu_time();
	
	// Build all element stiffness matrices
	if(verbose_) CHK_ERR(DomainPtr_->BuildGeometricElementMatrices(this));
	
	// Compute time for creating the matrix and print
	double fei_create_element_matrices_time = fei::utils::cpu_time() - fei_creatematrix_start_time;
	if (verbose_) {
		FEI_COUT << "DONE! " << FEI_ENDL;
		FEI_COUT << "CPU time to build geometric element matrices:   " << fei_create_element_matrices_time << FEI_ENDL;
	}
	
	// ============================================================ //
	// S E T  U P  L I N E A R  S Y S T E M  F O R  B U C K L I N G	//
	// ------------------------------------------------------------ //
	
	if(verbose_) FEI_COUT<< FEI_ENDL << "SETTING UP LINEAR SYSTEM FOR BUCKLING..." << FEI_ENDL;
			
	// Initialize timer to load the matrix
	double start_load_time = fei::utils::cpu_time();
	
	// Create StiffnesMatrix, DisplacementVector, and ForceVector based on matrixGraph
	GeomStiffnesMatrix_ 	= Factory_->createMatrix(MatrixGraph_);
	GeomDisplacementVector_ = Factory_->createVector(MatrixGraph_, true);
	GeomForceVector_  		= Factory_->createVector(MatrixGraph_);
	GeomLinSys_ 			= Factory_->createLinearSystem(MatrixGraph_);
	
	// Set GeomStiffnesMatrix, GeomDisplacementVector, and GeomForceVector Geomto LinSys
	GeomLinSys_->setMatrix(GeomStiffnesMatrix_);
	GeomLinSys_->setSolutionVector(GeomDisplacementVector_);
	GeomLinSys_->setRHS(GeomForceVector_);
	
	// Assemble GeomStiffnessMatrix_ (Use useGeomStiffness=true to 
	// indicate the geometric element matrices are to be used)
	if(verbose_) DomainPtr_->AssembleStiffnessMatrix(MatrixGraph_.get(),GeomStiffnesMatrix_.get(),true);
		
	// Sum random vector into GeomForceVector_ (used for buckling initialization)
	if(verbose_){
		int nDOF = DomainPtr_->NumTotalNodes * DomainConstants::DISP_FIELD_SIZE;
		double* initVec = new double[nDOF];
		InitialVector(nDOF,initVec);
		for(int id = 0; id<nDOF; id++){
			CHK_ERR( GeomForceVector_->copyIn(1, &id, &initVec[id], 0));
		}
	}
	
	// Apply constraints for this load case: id=GlobalLoadCaseID_
	if(verbose_){
		CHK_ERR(DomainPtr_->ApplyConstraints(GeomLinSys_.get(),GlobalLoadCaseID_));
	}
	
	// Loading of GeomLinSys is completed
	CHK_ERR( GeomLinSys_->loadComplete() );
	
	// Compute loading time
	double fei_load_time = fei::utils::cpu_time() - start_load_time;
	
	// Print output
	if (verbose_) {
		FEI_COUT << "DONE! " << FEI_ENDL;
		FEI_COUT << "CPU time to load element data into linear system for buckling:    " << fei_load_time << FEI_ENDL;
	}
	
	return 0;
}

// Solve LinearBuckling system
int LinearBucklingSolver::Solve(){
	
	if(verbose_)
	FEI_COUT << "------- LINEAR BUCKLING SOLVE - SOLVE -------" << FEI_ENDL;
	if (verbose_) FEI_COUT << FEI_ENDL<< "SOLVING BUCKLING SYSTEM..." << FEI_ENDL;
	
	// ======================================================== //
	// 	      S O L V E  B U C K L I N G  P R O B L E M 		//
	// -------------------------------------------------------- //
	
	double start_solve_time = fei::utils::cpu_time();
	
	// Compute number of DOF and share among processes
	int nDOF;
	if(verbose_) nDOF = DomainPtr_->NumTotalNodes * DomainConstants::DISP_FIELD_SIZE;
	CHK_ERR(MPI_Bcast(&nDOF,1,MPI_INT,0,MPI_COMM_WORLD));
	
	// Create ArpackOperations object
	ArpackOperations ao(Factory_,StiffnesMatrix_,GeomStiffnesMatrix_,nDOF,verbose_);
	
	// Initialize positiveness of stiffness matrix to be true
	bool isPosDef = true;
	
	// Perform the buckling analysis
	
	if(verbose_){
		if (verbose_) FEI_COUT <<  "Solving Eigenvalue Problem..." << FEI_ENDL;
		// First compute initialization vector
		double* initVec = new double[nDOF];
		for(int id=0;id<nDOF;id++) GeomForceVector_->copyOut(1, &id, &initVec[id], 0);		
		
		//--------------------------------------------------------------------//
		//------------ S O L V E  (B i A - l a m b d a * I ) x = 0 -----------//
		//--------------------------------------------------------------------//
		
		// Create the ARPACK++ solver for (A-lambda*i)x=0
		ProbSymStdEigB_ = new ARSymStdEig<double,ArpackOperations>(nDOF,nModes_,&ao,&ArpackOperations::OpBiA,
			"LM",0,0.0,0,initVec);
				
		//ProbSymStdEigB_->Trace();
		
		// Solve eigenvalue problem
		nconvStdB_ = ProbSymStdEigB_->FindEigenvectors();
		
		// Loop through eigenvalue to find a negative. In that case the
		// system was not positive
		for(std::size_t i=0;i<nconvStdB_;i++){
			if(ProbSymStdEigB_->Eigenvalue(i)<0) isPosDef = false;
		}
		
		/*
		//----------------------------------------------------------------//
		//------------ S O L V E  (A - l a m b d a * B ) x = 0 -----------//
		//----------------------------------------------------------------//
		std::cout << std::endl << "--------------- SOLVE (A-lambda*B)x=0 ---------------"<< std::endl<< std::endl;
		
		// Set Initial vector again
		for(int id=0;id<nDOF;id++) GeomForceVector_->copyOut(1, &id, &initVec[id], 0);		
		FEI_COUT << "initVec: " << FEI_ENDL;
		for(std::size_t i=0;i<nDOF;i++){ 
			//FEI_COUT << initVec[i]<< "  ";
		}
		FEI_COUT << FEI_ENDL;
		
		// Create the ARPACK++ solver for (A-lambda*B)x=0
		ProbSymGenEigAB_ = new ARSymGenEig<double , ArpackOperations,ArpackOperations>
		          (nDOF,nModes_,&ao,&ArpackOperations::OpBiA,&ao,&ArpackOperations::OpB,
			"LM",0,0.0,0,initVec);
			
		// Trace problem
		ProbSymGenEigAB_->Trace();
				
		// Solve eigenvalue problem
		nconvGenAB_ = ProbSymGenEigAB_->FindEigenvectors();
		*/
		
		FEI_COUT << " Done... " << FEI_ENDL;
		
	}
	
	
	// Compute solver time
	double solve_time = fei::utils::cpu_time()-start_solve_time;

	// Print solver time
	if (verbose_) {
		FEI_COUT << "DONE! " << FEI_ENDL;
		FEI_COUT << "CPU time for solving: " << solve_time << FEI_ENDL;
	}
	
	// Check positive definiteness
	CHK_ERR(MPI_Bcast(&isPosDef,1,MPI_INT,0,MPI_COMM_WORLD));
	if(!isPosDef){
		if(verbose_)
			FEI_COUT << "ERROR: THE EFFECTIVE STIFFNESS MATRIX FOR BUCKLING ANALYSIS IS NOT POSITIVE DEFINITE" << FEI_ENDL;
		return 1;
	}
	
	
	
	
	return 0;
}

// WriteOutput for LinearBuckling system
int LinearBucklingSolver::WriteOutput( bool writeLoadCaseInfo){
	/*
	if(verbose_)
	FEI_COUT << "---- LINEAR BUCKLING SOLVE - WRITE OUTPUT ---" << FEI_ENDL;
	
	if(verbose_){
		int nDOF = DomainPtr_->NumTotalNodes * DomainConstants::DISP_FIELD_SIZE;
		
		std::cout << std::endl << "------------- OUTPUT (BiA-lambda*I)x=0 --------------"<< std::endl<< std::endl;
		// Print eigenvalues
		for(std::size_t i=0;i<nconvStdB_;i++){
			std::cout << "Eigenvalue [" << i+1<<"] = ";
			std::cout << ProbSymStdEigB_->Eigenvalue(i) << std::endl;
			std::cout << "Eigenvector[" << i+1<<"] = ";
			for(std::size_t j=0;j<nDOF;j++){
				std::cout << ProbSymStdEigB_->Eigenvector(i,j) << "    ";
			}
			std::cout<<std::endl<<std::endl;
		}
				
		std::cout << std::endl << "-------------- OUTPUT (A-lambda*B)x=0 ---------------"<< std::endl<< std::endl;
		// Print eigenvalues
		for(std::size_t i=0;i<nconvGenAB_;i++){
			std::cout << "Eigenvalue [" << i+1<<"] = ";
			std::cout << ProbSymGenEigAB_->Eigenvalue(i) << std::endl;
			std::cout << "Eigenvector[" << i+1<<"] = ";
			for(std::size_t j=0;j<nDOF;j++){
				std::cout << ProbSymGenEigAB_->Eigenvector(i,j) << "    ";
			}
			std::cout<<std::endl<<std::endl;
		}
		
	}
	*/
	
	// Write output for Solve and Prestress Linear Static (if applicable)
	// has already been performed.
	
	// Write output for this Linear Buckling solve
	if(verbose_){ 
		CHK_ERR(DomainPtr_->WriteLinearBucklingOutput(this,OutputFile_));
	}
	
	return 0;
}


// Create random initial vector for buckling analysis
void LinearBucklingSolver::InitialVector(int n, double* initVec){
	
	// Create initial vector using random number generator
	std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    
	double number = 0, maxVal = 0;
	for(std::size_t i=0; i<n; i++){
		number = distribution(generator);
		initVec[i] = number;
		if(fabs(number)>maxVal) maxVal = fabs(number);
	}
	// Normalize wrt maxVal
	for(std::size_t i=0; i<n; i++)initVec[i]/=maxVal;
}
