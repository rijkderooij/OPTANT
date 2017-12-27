//
// This is a simple program to exercise FEI classes for the
// purposes of testing, code tuning and scaling studies.
//
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"

#include "fei_iostream.hpp" 

//Including the header fei_base.hpp pulls in the declaration for
//various classes in the 'fei' namespace.

#include "fei_base.hpp"


//Make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.{hpp,cpp}.

#include <fei_Factory_Trilinos.hpp>

//
//Include definitions of macros like 'CHK_ERR' to call functions and check
//the return code.
//
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

    
	
	// Finalize
	#ifndef FEI_SER
		MPI_Finalize();
	#endif

	return(EXIT_SUCCESS);

}
