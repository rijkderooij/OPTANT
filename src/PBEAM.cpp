#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "fei_base.hpp"

// Constructors:
PBEAM::PBEAM(): Property(Beam){
	for(std::size_t i=0;i<16;i++) BendingStiffness_[i]=0;
	for(std::size_t i=0;i<4;i++) ShearStiffness_[i]=0;
	for(std::size_t i=0;i<8;i++) CouplingStiffness_[i]=0;
	for(std::size_t i=0;i<4;i++) BendingExpansion_[i]=0;
	for(std::size_t i=0;i<2;i++) ShearExpansion_[i]=0;
}


// Read PBEAM from file
int PBEAM::ReadFromFile(std::ifstream& fin, int PropOption, 
						int NumLines, Domain& domain){
							
	char line[DomainConstants::MAXLINESIZE];
	char temp[DomainConstants::MAXNAMESIZE];
	int ierr=0;
	
	if(PropOption==0) ; // SECTION INPUT
	
	else if(PropOption==1){
				
		Material* mat;
		double e11;
		double g12,g13,g23;
		double alpha;
		double area;
		double kshearInv;
		
		// Next loop through the following numLines lines
		for (int i=0; i<NumLines; i++){
			// Get next line and parse it
			fin.getline(line,DomainConstants::MAXLINESIZE);	
			std::stringstream parse(line);
			
			// First line
			if(i==0){
				// First entry gives material id:
				int matID; parse>>matID;
				mat = domain.MaterialList[matID];
				
				// Material parameters
				e11 = (mat->EngineeringMaterialProperties())[0];
				g23 = (mat->EngineeringMaterialProperties())[7];
				g13 = (mat->EngineeringMaterialProperties())[7];
				g12 = (mat->EngineeringMaterialProperties())[8];
				alpha = (mat->ThermalExpansion())[0];
				
				// Second entry gives area, third 1/kShear
				if(parse>>area>>kshearInv){
					ShearStiffness_[0] = kshearInv/(area*g12);
					ShearStiffness_[3] = kshearInv/(area*g13);
					
					BendingExpansion_[0] = e11*area*alpha;}
				else ierr=1;
								
				// Check whether this line was completed correctly
				if(parse>>temp) ierr=1;
			}
			// Second line
			else if(i==1){
				// First three entries Ixx,Iyy,Izz are used
				double Jxx,Iyy,Izz;
				if(parse>>Jxx>>Iyy>>Izz){
					BendingStiffness_[0] = e11*area;
					BendingStiffness_[5] = g23*Jxx;
					BendingStiffness_[10] = e11*Iyy;
					BendingStiffness_[15] = e11*Izz;
				} else ierr=1;
				
				// Last three entries not used:
				if(parse>>temp>>temp>>temp); else  ierr=1;
				
				// Line should be complete
				if(parse>>temp) ierr=1;	  
			}
			
			
		}
	// Manual input of matrices:
	}else if(PropOption==2){
		double temp;
		// First two lines consist of BendingStiffnessMatrix
		
		// Get first line and parse it
		{
		fin.getline(line,DomainConstants::MAXLINESIZE);	
		std::stringstream parse(line);
		if(parse>>BendingStiffness_[0]>>BendingStiffness_[1]>>
		          BendingStiffness_[2]>>BendingStiffness_[3]>>
		          BendingStiffness_[5]>>BendingStiffness_[6]);
		else ierr=1;}
		
		// Get second line and parse it
		{
		fin.getline(line,DomainConstants::MAXLINESIZE);	
		std::stringstream parse(line);
		if(parse>>BendingStiffness_[7]>>BendingStiffness_[10]>>
		          BendingStiffness_[11]>>BendingStiffness_[15]);
		else ierr=1;
		
		if(parse>>temp)ierr=1; // Check if no extra info is given
		}
		// Third and fourth line consist of CouplingStiffness and 
		// ShearStiffness. Note that in C++ the matrix are given columns 
		// by column, whereas the input file has row by row
		
		// Get third line and parse it
		{
		fin.getline(line,DomainConstants::MAXLINESIZE);	
		std::stringstream parse(line);
		if(parse>>CouplingStiffness_[0]>>CouplingStiffness_[2]>>
		          CouplingStiffness_[4]>>CouplingStiffness_[6]>>
		          CouplingStiffness_[1]>>CouplingStiffness_[3]);
		else ierr=1;
		}
		// Get fourth line and parse it
		{
		fin.getline(line,DomainConstants::MAXLINESIZE);	
		std::stringstream parse(line);
		if(parse>>CouplingStiffness_[5]>>CouplingStiffness_[7]>>
		          ShearStiffness_[0]>>ShearStiffness_[1]>>
		          ShearStiffness_[3]);
		else ierr=1;
		
		if(parse>>temp)ierr=1; // Check if no extra info is given
		}
		// Get fifth line (if available) and parse it
		if(NumLines==5){
			fin.getline(line,DomainConstants::MAXLINESIZE);	
			std::stringstream parse(line);
			if(parse>>BendingExpansion_[0]>>BendingExpansion_[1]>>
					  BendingExpansion_[2]>>BendingExpansion_[3]>>
					  ShearExpansion_[0]>>ShearExpansion_[1]);
			else ierr=1;
			
			if(parse>>temp)ierr=1; // Check if no extra info is given
		}
		
		// Next fill BendingStiffness_ and ShearStiffness_ to be symmetric
		BendingStiffness_[4] = BendingStiffness_[1];
		BendingStiffness_[8] = BendingStiffness_[2];
		BendingStiffness_[9] = BendingStiffness_[6];
		BendingStiffness_[12] = BendingStiffness_[3];
		BendingStiffness_[13] = BendingStiffness_[7];
		BendingStiffness_[14] = BendingStiffness_[11];
		
		ShearStiffness_[2] = ShearStiffness_[1];
		
		
	}else ierr=1;
	
	
	return ierr;
}



