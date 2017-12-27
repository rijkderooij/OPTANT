#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "fei_base.hpp"
#include <fei_Factory_Trilinos.hpp>
#include "fei_ErrMacros.hpp"


// Constructors
Solver::Solver(Domain* domPtr, int globalLoadCaseID, const char* outputFile, bool verbose) 
	: DomainPtr_(domPtr), GlobalLoadCaseID_(globalLoadCaseID), verbose_(verbose),
	OutputFile_(outputFile) {
}

// Destructor
Solver::~Solver(){
	Factory_.reset();
}


// Initialize Solver
int Solver::Initialize(MPI_Comm &comm){
	
	// ======================================================== //
	// 			I N I T I A L I Z E  F A C T O R Y				//
	// -------------------------------------------------------- //
	
	// Declare a shared pointer: The shared_ptr class template stores a 
	// pointer to a dynamically allocated object, typically with a C++ 
	// new-expression. The object pointed to is guaranteed to be deleted 
	// when the last shared_ptr pointing to it is destroyed or reset.
	// The fei::Factory is an interface for creating fei::instances
	// fei::SharedPtr<fei::Factory> factory;
	
	feiFactoryPtr factory(new Factory_Trilinos(comm));
	fei::ParameterSet paramset;
	paramset.add(fei::Param("Trilinos_Solver", "Amesos"));
	factory->parameters(paramset);
	factory->parameters(DomainPtr_->SolverParams);
	
	// Check whether the factory allocation was succesful
	if (factory.get() == NULL) {
		FEI_COUT << "fei::Factory allocation failed." << FEI_ENDL;

		#ifndef FEI_SER
		MPI_Finalize();
		#endif

		return(-1);
	}
	
	// Set factory_ as factory
	Factory_ = factory;
	
	
	// ======================================================== //
	// 		I N I T I A L I Z E  M A T R I X G R A P H			//
	// -------------------------------------------------------- //
	
	if(verbose_) FEI_COUT<< FEI_ENDL << "INITIALIZING MATRIX GRAPH..." << FEI_ENDL;
	
	// Start time for initializing computations
	double start_init_time = fei::utils::cpu_time();
	
	// Create rowspace, which is equalt to column space
	NodeSpace_ = Factory_->createVectorSpace(comm,NULL); 	//row space
	feiVectorSpacePtr colSpace ; // == NULL, hence == row space
	
	// Create matrixgraph from these space. Note: if the matrix is 
	// symmetric, then the columns space is not required,i.e. NULL
	MatrixGraph_ = Factory_->createMatrixGraph(NodeSpace_, colSpace,"StiffnessMatrix"); 
		
	//load some solver-control parameters.
	NodeSpace_->setParameters(DomainPtr_->SolverParams);
	MatrixGraph_->setParameters(DomainPtr_->SolverParams);
		
	// Define displacement/rotation field, with 6 dof per node. In case
	// only displacement should be considered (i.e. for solids), the 
	// rotational dof should be deactivated using slave constraints.
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int dispFieldSize = DomainConstants::DISP_FIELD_SIZE;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int constraintTypeID = 2;
	
	NodeSpace_->defineFields( 1, &dispFieldID, &dispFieldSize );
	NodeSpace_->defineIDTypes(1, &nodeTypeID );
	NodeSpace_->defineIDTypes(1, &constraintTypeID );
	
	// Initialize element connectivities
	if(verbose_)
		CHK_ERR(DomainPtr_->InitializeElementConnectivities(MatrixGraph_.get()));
	
	
	// Initialize constraints for first load case: id=0
	if(verbose_){
		CHK_ERR(DomainPtr_->InitializeConstraints(MatrixGraph_.get(),GlobalLoadCaseID_));
	}
	
	// Initialization is complete
	CHK_ERR(MatrixGraph_->initComplete());

	// Compute cpu time used for initialization
	double fei_init_time = fei::utils::cpu_time() - start_init_time;

	if (verbose_) 
		FEI_COUT << "DONE!" << FEI_ENDL << "CPU time for matrixgraph initialization:   " << fei_init_time << FEI_ENDL;
	
	return 0;
}

// Write output
int Solver::WriteOutput(bool writeLoadCaseInfo){
	
	if(verbose_){
		
		// Write general Loadcase info
		LoadCase* loadcase = DomainPtr_->LoadCaseList[GlobalLoadCaseID_];
		
		// Create output stream
		std::ofstream fout;
		
		// Check whether the output file has been opened already for a previous
		// load case. If not, write new file; if yes, append to file;
		if(!DomainPtr_->OutputReq->OutputFileAlreadyOpened_){
			fout.open(OutputFile_, std::ofstream::out);
			DomainPtr_->OutputReq->OutputFileAlreadyOpened_ = true;
		}else fout.open(OutputFile_, std::ofstream::app);
		
		// Write general load case info
		fout << "------ LOADCASE " << std::setw(3)
			 << GlobalLoadCaseID_ << " ------" << FEI_ENDL
			 << "SPC set:   " << loadcase->SPCSetID()  << FEI_ENDL
			 << "MPC set:   " << loadcase->MPCSetID()  << FEI_ENDL
			 << "LOAD set:  " << loadcase->LOADSetID() << FEI_ENDL
			 << "PLOAD set: " << loadcase->PLOADSetID() << FEI_ENDL
			 << "TEMP set:  " << loadcase->TEMPSetID() << FEI_ENDL << FEI_ENDL;
			 
		// Close outputFile
		fout.close();
	}
		
	
	
	return 0;
	
}
