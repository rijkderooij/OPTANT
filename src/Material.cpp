#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 



Material::Material(){
}

Material::Material(const char* name, int type, double* propArray)
				 : MaterialType_((Type) type){
	 
	 // Copy name into MaterialName_
	 strcpy(MaterialName_, name);
	 
	 // Switch for materialType
	 switch(MaterialType_){
		 case Isotropic:
			IsotropicMaterial(propArray);
			break;
		 default:
			FEI_COUT << " My type is not included " << FEI_ENDL;
			break;
		}
 }
 
 
 int Material::IsotropicMaterial(double* propArray){
	 Density_ = propArray[0];
	 ReferenceTemperature_ = propArray[3];
	 
	 // EngineeringMatProperties
	 double E = propArray[1];
	 double nu = propArray[2];
	 double G = E/2/(1+nu);
	 double emp[] = {E,E,E,nu,nu,nu,G,G,G};
	 for(int i=0;i<9;i++) EngineeringMaterialProperties_[i]=emp[i];
	 
	 // MaterialStiffness
	 double C1 = E/(1+nu)/(1-2*nu);
	 double cii = C1*(1-nu), cij = C1*nu, gii = C1*(1-2*nu)/2;
	 double ms[] = {cii,cij, cij,0,0,0, cij,cii,cij,0,0,0, cij,cij,cii,0,0,0,
				   0,0,0,gii,0,0, 0,0,0,0,gii,0, 0,0,0,0,0,gii};
	 for(int i=0;i<9;i++) MaterialStiffness_[i]=ms[i];
	 
	 // PlaneStressStiffness
	 double Cpss = E/(1-nu*nu);
	 double pss[] = {Cpss,Cpss*nu,0, Cpss*nu,Cpss,0, 0,0,Cpss*(1-nu)/2};
	 for(int i=0;i<9;i++) PlaneStressStiffness_[i]=pss[i];
	 
	 // PlaneStrainStiffness
	 double psn[] = {cii,cij,0, cij,cii,0, 0,0,gii};
	 for(int i=0;i<9;i++) PlaneStrainStiffness_[i]=psn[i];
	 
	 // ThermalExpansion
	 double alpha = propArray[4];
	 double te[] = {alpha,alpha,alpha,0,0,0};
	 for(int i=0;i<6;i++) ThermalExpansion_[i]=te[i];
	 
	 // ThermalConductivity (general: anisotropic diffusion equation)
	 double lambda = propArray[5];
	 double tc[] = {lambda,0,0,0,lambda,0,0,0,lambda};
	 for(int i=0;i<9;i++) ThermalConductivity_[i]=tc[i];
 }
