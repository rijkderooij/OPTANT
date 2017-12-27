//
// This is a simple program to exercise FEI classes for the
// purposes of testing, code tuning and scaling studies.
//

#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "fei_iostream.hpp" 
#include "fei_base.hpp"
#include <fei_Factory_Trilinos.hpp>
#include "fei_ErrMacros.hpp"


// Include own classes
#include "include/PCH_OPTANT.h"

// ================================================ //
// 				M A I N  F U N C T I O N 			//
// ------------------------------------------------ //

int main(int argc, char** argv ) {
	using std::cout;
	using std::endl;
	
	// ======================================================== //
	// 			I N I T I A L I Z E  M P I  D A T A 			//
	// -------------------------------------------------------- //
		
	int numProcs, localProc;
	#ifndef FEI_SER
	  CHK_ERR(MPI_Init(&argc, &argv));
	  CHK_ERR(MPI_Comm_rank(MPI_COMM_WORLD, &localProc));
	  CHK_ERR(MPI_Comm_size(MPI_COMM_WORLD, &numProcs));
	#else
	  localProc = 0;
	  numProcs = 1;
	#endif
	
	MPI_Comm comm = MPI_COMM_WORLD;
	bool verbose = (localProc==0);
	
	// ======================================================== //
	// 			R E A D  I N P U T  P A R A M E T E R S			//
	// -------------------------------------------------------- //
	
	// Start timer
	double start_time = fei::utils::cpu_time();
	
	
	// Print processes
	if(verbose)
		FEI_COUT << endl << endl << "Total number of processes:  " 
		     << numProcs << FEI_ENDL;
	
	// Read names Inputfile, solverParameters, and OutputRequest name
	char InputFile[DomainConstants::MAXLINESIZE];
	char SolverParamFile[DomainConstants::MAXLINESIZE];
	char OutputReqFile[DomainConstants::MAXLINESIZE];
	char OutputFile[DomainConstants::MAXLINESIZE];
	
	CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-inp",InputFile));
	CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-sol",SolverParamFile));
	CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-req",OutputReqFile));
	//CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-out",OutputFile));
	
	// Outputfile has same name as input file, but different file format:
	CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-inp",OutputFile));
	char * pch = strstr (OutputFile,".");
	strncpy (pch+1,"c06",3);
	
	// Create domain on heap memory from input file
	Domain domain(InputFile);
	
	if(verbose){
		FEI_COUT << "InputFile: 		" << InputFile << FEI_ENDL;
		FEI_COUT << "SolverParamFile: 	" << SolverParamFile << FEI_ENDL;
		FEI_COUT << "OutputReqFile is: 	" << OutputReqFile << FEI_ENDL;
		FEI_COUT << "OutputFile: 		" << OutputFile << FEI_ENDL;

		// Read input file
		CHK_ERR(domain.ReadInputFile());
	}

	// Read solver parameters and outputRequests
	std::vector<std::string> solParamString;
	std::vector<std::string> outReqString;
	CHK_ERR(UtilityFunctions::read_param_file(SolverParamFile,
												comm,solParamString));
	CHK_ERR(UtilityFunctions::read_param_file(OutputReqFile,
												comm,outReqString));
	
	if(verbose){
		cout <<endl << "Solver Parameters: " << endl;
		for(int i=0;i<solParamString.size();++i){
			cout << "  " << i << ":\t " << solParamString.at(i) << endl;			
		}
		
		cout <<endl << "Output Requests: " << endl;
		for(int i=0;i<outReqString.size();++i){
			cout << "  " << i << ":\t " << outReqString.at(i) << endl;			
		}
	}
	
	// Put the input in a parameterset. Note the parse_strings() function
	// assumes each entry in stdstrings is a key-value combination 
	// seperated by a " " 
	fei::ParameterSet solverParams, outputParams;
    fei::utils::parse_strings(solParamString, " ", domain.SolverParams);
    fei::utils::parse_strings(outReqString, " ", outputParams);
    
    // Create an OutputRequest object and set parameters
	//OutputRequest OutputReq;
	domain.OutputReq = new OutputRequest();
	domain.OutputReq->SetParameters(outputParams);
		
	// ======================================================== //
	// 			 S T A R T  L O A D C A S E  L O O P			//
	// -------------------------------------------------------- //
	
	// Broadcast the number of loadcases and solutionType to all processes
	int nLoadCases, solutionType;
	if(verbose) nLoadCases = domain.LoadCaseList.size();
	if(verbose) solutionType = domain.SolutionType;
	
	CHK_ERR(MPI_Bcast(&nLoadCases,1,MPI_INT,0,MPI_COMM_WORLD));
	CHK_ERR(MPI_Bcast(&solutionType,1,MPI_INT,0,MPI_COMM_WORLD));
		
	// Loop through load cases. NOTE: the matrixgraph is the same for 
	// each load case. This cannot be the case if we would initialize
	// SPC or MPC with initSlaveConstraint(). In that case, the
	// matrixgraphs needs to be computed inside the loop
	for(size_t lcLID=0; lcLID<nLoadCases;lcLID++){
		
		// Compute GlobalLoadCaseID and share among processes
		int GlobalLoadCaseID = -1;
		if(verbose){
			std::map<const int, LoadCase*>::iterator itrLC = domain.LoadCaseList.begin();
			std::advance(itrLC,lcLID);
			GlobalLoadCaseID = (*itrLC).first;
		}
		CHK_ERR(MPI_Bcast(&GlobalLoadCaseID,1,MPI_INT,0,MPI_COMM_WORLD));
		
		if(verbose) FEI_COUT<< FEI_ENDL << "--------------  LOAD CASE " 
							<< std::setw(3) << GlobalLoadCaseID
							<< "  ---------------" << FEI_ENDL;
							
		// ======================================================== //
		// 					C R E A T E  S O L V E R				//
		// -------------------------------------------------------- //
		Solver* theSolver = NULL; 
		
		switch(solutionType){
			case Domain::LinearStatic:
			
				if(verbose) FEI_COUT << "   -------------  LINEAR STATIC  ------------   " << FEI_ENDL;
	
				theSolver = new LinearStaticSolver(&domain,GlobalLoadCaseID,OutputFile,verbose);
				break;
			case Domain::LinearBuckling:
				{
				// Check if a prestress loadcase is applied
				int GlobalPSLoadCaseID = 0;
				if(verbose) domain.SolverParams.getIntParamValue("preStress", GlobalPSLoadCaseID); 
				CHK_ERR(MPI_Bcast(&GlobalPSLoadCaseID,1,MPI_INT,0,MPI_COMM_WORLD));
				
				// Only initiate solver when prestress loadcase is not the same as buckling force
				if(GlobalPSLoadCaseID!=GlobalLoadCaseID){
					theSolver = new LinearBucklingSolver(&domain,GlobalLoadCaseID,GlobalPSLoadCaseID,OutputFile,verbose);
				
					if(verbose && GlobalPSLoadCaseID==0) 
						FEI_COUT << "-------------  LINEAR BUCKLING  --------------" << FEI_ENDL;
							
					if(verbose && GlobalPSLoadCaseID!=0) 
						FEI_COUT << "------  LINEAR BUCKLING WITH PRESTRESS  ------" << FEI_ENDL;
				}else{
					if(verbose) 
						FEI_COUT << "---------  PRESTRESS - NOT COMPUTED  ---------" << FEI_ENDL;
				}
				break;
				}
			default:
				break;
		}
		
		// Continue to next load case if theSolver has not been initialized
		if(theSolver==NULL) continue;
		
		// ======================================================== //
		// 				I N I T I A L I Z E  S O L V E R			//
		// -------------------------------------------------------- //
		
		// Initialize solver (this includes creating the fei::factory 
		// and fei::matrixgraph)
		if(theSolver->Initialize(comm)!=0){
			delete theSolver;
			return(EXIT_FAILURE);
		}
		
		// ======================================================== //
		// 					P R E P A R E  S O L V E R				//
		// -------------------------------------------------------- //
		
		// Prepare the solve (this includes e.g. setting up the linear 
		// system to be solved
		if(theSolver->Prepare()!=0){
			delete theSolver;
			return(EXIT_FAILURE);
		}
		
		// ======================================================== //
		// 				    E X E C U T E  S O L V E				//
		// -------------------------------------------------------- //
		
		// Solve the prepared system
		if(theSolver->Solve()!=0){
			delete theSolver;
			return(EXIT_FAILURE);
		}
		
		// ======================================================== //
		// 					W R I T E  O U T P U T					//
		// -------------------------------------------------------- //
		
		// Write output for the solved system
		if(theSolver->WriteOutput(OutputFile)!=0){
			delete theSolver;
			return(EXIT_FAILURE);
		}
		
		// ======================================================== //
		// 				 F I N I S H  L O A D C A S E				//
		// -------------------------------------------------------- //
		
		// Reset all element temperatures in the domain
		if(verbose)
			domain.ResetElementTemperatures();
				
		// Delete the solver
		delete theSolver;
	
    }
    
    // ======================================================== //
	// 			 E N D  O F  L O A D C A S E  L O O P			//
	// -------------------------------------------------------- //
    
	
	// Finalize
	#ifndef FEI_SER
		MPI_Finalize();
	#endif

	return(EXIT_SUCCESS);

}
