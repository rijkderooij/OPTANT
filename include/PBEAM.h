#ifndef PBEAM_H
#define PBEAM_H

#include "Property.h"

#include "fei_iostream.hpp" 

//Forward declarations
class Domain;
class Material;

// NOTE: THIS CLASS CONTAINS GENERAL STIFFNESS MATRICES CTILDE, GTILDE 
// AND STILDE, WHICH MAY REPRESENT THE EULER BERNOULLI OR THE TIMOSHENKO
// BEAM, BUT MAY ALSO CONTAIN A DIFFERENT TYPE

class PBEAM : public Property{
	
	public:		
		// Constructor
		PBEAM();
		
		// Functions
		int ReadFromFile(std::ifstream& fin, int PropOption,
						   int NumLines, Domain& domain);
		
		// Getters		
		double* BendingStiffness() {return BendingStiffness_;};
		double* ShearStiffness() {return ShearStiffness_;};
		double* CouplingStiffness() {return CouplingStiffness_;};
		double* BendingExpansion() {return BendingExpansion_;};
		double* ShearExpansion() {return ShearExpansion_;};
		
		
		
		
	private:
		//Variables
		Material *material_;
		double BendingStiffness_[16]; // 4x4 Ctilde matrix (bending and normal)
		double ShearStiffness_[4];    // 2x2 Stilde matrix (shear)
		double CouplingStiffness_[8]; // 2x4 coupling between Ctilde and Stilde
		
		double BendingExpansion_[4];  // 4x1 ctilde vector 
		double ShearExpansion_[2];    // 2x1 gtilde vector
		
		
	
};
#endif
