#ifndef OUTPUTREQUEST_H
#define OUTPUTREQUEST_H

#include "fei_base.hpp"
#include "fei_fstream.hpp"

// Forward declarations
class Domain;

class OutputRequest{
	
	public:
		OutputRequest(): OutputFileAlreadyOpened_(false),
						 PrintModelInfo_(false),
						 PrintNodalOutput_(false),
						 PrintElementOutput_(false),
						 PrintElementStrains_(false),
						 PrintElementForces_(false),
						 PrintElementStrainEnergy_(false),
						 PrintNodalBucklingOutput_(false),
						 numNodalOutputs_(0),
						 numElementOutputs_(0){};
		
		// Functions
		void SetParameters(fei::ParameterSet& params);
		
		
		// General write functions which do not belong to particular class
		void PrintModelInfo(Domain &domain, std::ofstream& fout);
		
		// Flag to check whether output file was opened already 
		bool OutputFileAlreadyOpened_;
		
		// Variables
		bool PrintModelInfo_;
		bool PrintNodalOutput_;
		bool PrintElementOutput_;
		bool PrintElementStrains_;
		bool PrintElementForces_;
		bool PrintElementStrainEnergy_;
		
		bool PrintNodalBucklingOutput_;
		
		int numNodalOutputs_;		// Used to determine the number of temporary files generated when writing output
		int numElementOutputs_;
		
		
	private:
	
	
};
#endif
