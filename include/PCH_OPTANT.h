// PCH_OPTANT.H 
#ifndef PCH_OPTANT
#define PCH_OPTANT

// ASCM_CPP classes
#include "../include/Domain.h"
#include "../include/Material.h"
#include "../include/Section.h"
#include "../include/Property.h"
#include "../include/PBEAM.h"
#include "../include/PSHELL.h"
#include "../include/Node.h"
#include "../include/Element.h"
#include "../include/CBEAM.h"
#include "../include/CTRIA.h"
#include "../include/CQUAD.h"
#include "../include/LoadCase.h"
#include "../include/LOAD.h"
#include "../include/PLOAD.h"
#include "../include/SPC.h"
#include "../include/MPC.h"
#include "../include/TEMP.h"
#include "../include/MatrixOperations.h"
#include "../include/UtilityFunctions.h"
#include "../include/OutputRequest.h"
#include "../include/ArpackOperations.h"
#include "../include/Solver.h"
#include "../include/LinearStaticSolver.h"
#include "../include/LinearBucklingSolver.h"

// Type definitions
typedef fei::SharedPtr<fei::Factory> feiFactoryPtr;
typedef fei::SharedPtr<fei::MatrixGraph> feiMatrixGraphPtr;
typedef fei::SharedPtr<fei::VectorSpace> feiVectorSpacePtr;

typedef fei::SharedPtr<fei::Vector> feiVectorPtr;
typedef fei::SharedPtr<fei::Matrix> feiMatrixPtr;
typedef fei::SharedPtr<fei::LinearSystem> feiLinearSystemPtr;

typedef fei::SharedPtr<fei::Solver> feiSolverPtr;


#endif // PCH_ASCM_H
