#ifndef PROPERTY_H
#define PROPERTY_H

#include "fei_iostream.hpp" 

//Forward declarations
class Domain;
class Section;

class Property{
	
	public:	
		enum Type{ 
				Rod=0,
				Beam=1,
				Shell=2,
				Solid=3
			};	
			
		Property(Property::Type pt);
			
		
		
		// General Functions
		virtual Section* GetSection() {return Section_;};
		virtual int SetSection(Section *s) {Section_ = s;};
		
		virtual int ReadFromFile(std::ifstream& fin, 
		            int PropOption, int NumLines, Domain& domain) = 0;
		
		Property::Type PropertyType(){return PropertyType_;};
		
		
	protected:
		Section* Section_;		
		Property::Type PropertyType_;
		
		//Functions
	
	
};
#endif
