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
		cout << endl << endl << "The total number of processes is:  " 
		     << numProcs << endl;
	cout << " My process is : " << localProc << endl;
	
	// Read names Inputfile, solverParameters, and OutputRequest name
	char InputFile[DomainConstants::MAXLINESIZE];
	
	CHK_ERR(UtilityFunctions::getFileName(argc, argv,"-inp",InputFile));
	
	// Create domain on heap memory from input file
	Domain domain(InputFile);
	
	if(verbose){
		FEI_COUT << " The InputFile is: " << InputFile << FEI_ENDL;

		// Read input file
		CHK_ERR(domain.ReadInputFile());
	}
	
	
	
	// ======================================================== //
	// 		C O M P U T E  E L E M E N T  M A T R I C E S		//
	// -------------------------------------------------------- //
	
	if(verbose){
		// We will test the first element of each block
		std::vector<Element*>::iterator itrE = domain.ElementList.begin(); // Iterator for the element to compute
		
		// Allocate some space to store the element matrix
		int maxDOFPerElem =(*(std::max_element(domain.NodesPerElementBlock.begin(), 
							domain.NodesPerElementBlock.end())))
							*DomainConstants::DISP_FIELD_SIZE;
		double* elementMat = new double[maxDOFPerElem*maxDOFPerElem];
		double* elementMatAccess = new double[maxDOFPerElem*maxDOFPerElem]; // Used to access the elementMatrix
		int elementSize;
		
		FEI_COUT << FEI_ENDL << FEI_ENDL << "The maximum number of dof per element is: " << 
				  maxDOFPerElem << FEI_ENDL;
				  
		// We will try to compute the number of elements generated per second, hence introduce some variables
		int num_elem_to_build;				// Number of elements to build
		double factor_to_build;				// Factor to increase num_of_elem_to_build in next iteration
		double start_time;					// Start time
		double build_time;		 			// Time required to build
		int iter;							// Iteration counter
		
		// Loop through the element blocks
		for (std::vector<int>::iterator itrB = domain.NumElementsBlock.begin();
			itrB != domain.NumElementsBlock.end(); ++itrB){
			
			// Set pointer of element matrix
			(*itrE)->SetElementMatrixLocation(elementMat);
			
			// Print info of this element
			int GlobalElID = domain.ElementID_LG[(*itrE)->LocalElementID()];
			FEI_COUT << FEI_ENDL << "Test build time of element " 
			         << GlobalElID << ":" << FEI_ENDL;
			
			
			// Set variables
			num_elem_to_build=10;		// Start by building 100 elements
			build_time = 0.0;
			iter = 0;	
			
			while(build_time<1.0 && iter<50){		
				// Increment iter
				iter++;
				
				// Compute the element matrix num_of_elem_to_build times
				start_time = fei::utils::cpu_time(); 	// Set start time
				for(std::size_t nBuild=0; nBuild<num_elem_to_build; nBuild++){
					CHK_ERR((*itrE)->BuildElementMatrix());
				}
				build_time = fei::utils::cpu_time() - start_time;
				
				// Print time
				FEI_COUT << " #" << iter <<"\t - Build Time: " << build_time <<
							"\t - Num of Elem: " << num_elem_to_build <<
							"\t - Speed: " << num_elem_to_build/build_time 
						 << "\telem/s" << FEI_ENDL;
						 
				// Update number of elements to build by aiming at a total build time of 0.8 seconds (to prevent overshoot over 1 sec)
				factor_to_build = 0.8/build_time;
				if(factor_to_build>5) factor_to_build=5;
				if(factor_to_build<1.1) factor_to_build = 1.1;
				num_elem_to_build *= factor_to_build;
			}
			
			
			// Access the element matrix
			//(*itrE)->ElementMatrix(elementMatAccess,elementSize);
			//FEI_COUT<< "ElementMatrix:" << FEI_ENDL; UtilityFunctions::PrintMatrix(elementMatAccess,elementSize,elementSize);
	
			
			// Increase element iterator by number of elements in this block
			itrE += *itrB;
		}
	
	
	}
	
	
	

	#ifndef FEI_SER
		MPI_Finalize();
	#endif

	return(EXIT_SUCCESS);

}
