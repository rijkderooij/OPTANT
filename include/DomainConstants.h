#ifndef DOMAINCONSTANTS_H
#define DOMAINCONSTANTS_H


#include <cmath>

namespace DomainConstants{
	
	const int MAXNAMESIZE = 20;
	const int MAXLINESIZE = 100;
	
	// Constants regarding the allowed angle between the longitudinal and
	// normal direction of a beam element. This angle should not be too
	// small. These checks are made when reading CBEAM element
	const double WARNINGANGLE = 30*atan(1)*4/180; 	// < 30deg -> warning
	const double ERRORANGLE   = 1*atan(1)*4/180;	// < 1deg -> error
	
	
	// Displacement/rotation field id and size
	const int DISP_FIELD_ID = 0;
	const int DISP_FIELD_SIZE = 6;
	
	// Type ID's
	const int NODE_TYPE_ID = 0;
	
	
	// OUTPUTREQUEST
	const int OUTPUT_PRECISION = 9;
	const int OUTPUT_DOUBLE_FIELD_WIDTH = 17;
	const int OUTPUT_INT_FIELD_WIDTH = 6;
	
}

#endif
