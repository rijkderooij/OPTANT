#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 

// Constructors
LoadCase::LoadCase(){
	std::cout << " I created a LoadCase" << std::endl;
}

LoadCase::LoadCase(int spc, int mpc, int load, int pload, int temp): 
    SPCSetID_(spc), MPCSetID_(mpc), LOADSetID_(load), 
    PLOADSetID_(pload), TEMPSetID_(temp){};
	
// Initialize loadcase
int LoadCase::Initialize(Domain& domain){return 0;}

// Apply loadcase
int LoadCase::Apply(Domain& domain){return 0;}

