#ifndef MATERIAL_H
#define MATERIAL_H

#include "fei_iostream.hpp"
#include <string.h>

#include "DomainConstants.h"

class Material{
	
	public:
		enum Type{ 
				Isotropic=0,
				Orthotropic2D=1,
				Anisotropic2D=2,
				Orthotropic3D=3
			};
			
		// Constructors
		Material();
		Material(const char* name, int type, double* propArray);
		
		// Getters
		const char* MaterialName(){return MaterialName_;};		
		Material::Type MaterialType(){return MaterialType_;};
		double Density(){return Density_;};
		double ReferenceTemperature() {return ReferenceTemperature_;};
		double* EngineeringMaterialProperties() {return EngineeringMaterialProperties_;};
		double* MaterialStiffness() {return MaterialStiffness_;};
		double* PlaneStressStiffness() {return PlaneStressStiffness_;};
		double* PlaneStrainStiffness() {return PlaneStrainStiffness_;};
		double* ThermalExpansion() {return ThermalExpansion_;};
		double* ThermalConductivity() {return ThermalConductivity_;};
		
	private:
		
			
		char MaterialName_[DomainConstants::MAXNAMESIZE];		
		
		Material::Type MaterialType_;
		double Density_;
		double ReferenceTemperature_;
		double EngineeringMaterialProperties_[9];
		double MaterialStiffness_[36];
		double PlaneStressStiffness_[9];
		double PlaneStrainStiffness_[9];
		double ThermalExpansion_[6];
		double ThermalConductivity_[9];
		
		//Functions
		int IsotropicMaterial(double* propArray); //Handle isotropric
		
	
	
};
#endif
