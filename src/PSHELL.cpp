#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "fei_base.hpp"


// Constructors:
PSHELL::PSHELL(): Property(Shell){
	for(int i=0;i<6;i++) A_[i]=0;
	for(int i=0;i<6;i++) B_[i]=0;
	for(int i=0;i<6;i++) D_[i]=0;
	for(int i=0;i<3;i++) a_[i]=0;
	for(int i=0;i<3;i++) b_[i]=0;
}


// Read PSHELL from file
int PSHELL::ReadFromFile(std::ifstream& fin, int PropOption, 
						int NumLines, Domain& domain){
	
						
	char line[DomainConstants::MAXLINESIZE];
	char temp[DomainConstants::MAXNAMESIZE];
	int ierr=0;
	
	if(PropOption==0) ; // SECTION INPUT
	
	else if(PropOption==1){
		
		// Check that numLines should be 3 or 4
		if(NumLines!=3 && NumLines!=4){
			FEI_COUT << "ERROR: NUMLINES IN PROPERTY NOT CORRECT" << FEI_ENDL;
			return 1;
		} 		
		
		// Next loop through the following numLines lines
		for (int i=0; i<NumLines; i++){
			// Get next line and parse it
			fin.getline(line,DomainConstants::MAXLINESIZE);	
			std::stringstream parse(line);
			
			// First line
			if(i==0){
				// This lines fills A_ 
				if(parse>>A_[0]>>A_[1]>>A_[2]>>A_[3]>>A_[4]>>A_[5]);
				else ierr=1;
				
				// Line should be complete
				if(parse>>temp) ierr=1;	  
			}
			// Second line
			else if(i==1){
				// This lines fills B_ 
				if(parse>>B_[0]>>B_[1]>>B_[2]>>B_[3]>>B_[4]>>B_[5]);
				else ierr=1;
				
				// Line should be complete
				if(parse>>temp) ierr=1;	  
			}
			// Third line
			else if(i==2){
				// This lines fills D_ 
				if(parse>>D_[0]>>D_[1]>>D_[2]>>D_[3]>>D_[4]>>D_[5]);
				else ierr=1;
				
				// Line should be complete
				if(parse>>temp) ierr=1;	  
			}
			// Fourth line (if present)
			else if(i==3){
				
				// This lines fills a_ and b_
				if(parse>>a_[0]>>a_[1]>>a_[2]>>b_[0]>>b_[1]>>b_[2]);
				else ierr=1;
				
				// Line should be complete
				if(parse>>temp) ierr=1;	  
			}
			
		}
	}
	
	return ierr;
}



