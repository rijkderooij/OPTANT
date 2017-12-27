#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 



// Constructor with type initialization
Element::Element(Element::Type et): ElementType_(et),
									NumDeltaTemperatures_(0),
									DeltaTemperatures_(NULL){}

Element::~Element(){delete [] Properties_; }


void Element::SetNumProperties(int np){
	NumProperties_=np;
	Properties_ = new Property*[np];
}

void Element::SetNumDeltaTemperatures(int nt){
	NumDeltaTemperatures_=nt;
	DeltaTemperatures_ = new double[nt];
}

void Element::ResetDeltaTemperatures(){
	NumDeltaTemperatures_=0;
	if(DeltaTemperatures_!=NULL){
		delete [] DeltaTemperatures_;
		DeltaTemperatures_ = NULL;
	}
}

int Element::AddProperty(int index, int GlobalPropID, Domain& domain){
	if(index<0 || index>NumProperties_){
		FEI_COUT << " ERROR: INDEX FOR PROPERTY " << GlobalPropID 
				 << " OUT OF BOUNDS" <<FEI_ENDL;
		return 1;
	}
	
	// Compute localPropID
	int locPropID = domain.PropertyID_GL[GlobalPropID];
	
	// Set property
	if(locPropID>-1 && locPropID<domain.PropertyList.size())
		Properties_[index] = domain.PropertyList[locPropID];
	else{
		FEI_COUT << " ERROR: COULD NOT SET PROPERTY " << GlobalPropID
				 << " TO ELEMENT " <<FEI_ENDL;
		return 1;};
		
	return 0;
}

int Element::AddDeltaTemperature(int index, double dTemp){
	if(index<0 || index>NumDeltaTemperatures_){
		FEI_COUT << " ERROR: INDEX FOR DELTA TEMPERATURE " 
				 << " OUT OF BOUNDS" <<FEI_ENDL;
		return 1;
	}
	
	// Set Delta Temperature
	DeltaTemperatures_[index]=dTemp;
		
	return 0;
}



int Element::ReadPropertyLines(std::ifstream& fin, Domain& domain, Property::Type propType){
	char line[DomainConstants::MAXLINESIZE];
	// Read line containing properties
	int numProp,propID;
	fin.getline(line,DomainConstants::MAXLINESIZE);
	std::stringstream parse(line);
	if(parse>>numProp); 
	else {FEI_COUT << " ERROR: COULD NOT READ # OF PROPERTIES " << FEI_ENDL;
		return 1;}
	
	SetNumProperties(numProp);
	
	// Read properties and add to property array
	for (int i=0;i<numProp;i++){
		if(parse>>propID){
			if(AddProperty(i,propID,domain)!=0) return 1;
			// Check property type
			if(Properties_[i]->PropertyType()!=propType){
				FEI_COUT << " ERROR: PROPERTY " << propID <<
				" INCONSISTENT WITH ELEMENT TYPE" << FEI_ENDL;
				return 1;
			}
		}
		else {FEI_COUT << " ERROR: COULD NOT READ PROPERTY ID " << FEI_ENDL;
			return 1;}
	}
	
	return 0;
}
