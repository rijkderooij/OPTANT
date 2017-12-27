#ifndef LOADCASE_H
#define LOADCASE_H

class Domain;

class LoadCase{
	
	public:
		LoadCase();
		LoadCase(int spc, int mpc, int load, int pload, int temp);
		
		int Initialize(Domain& domain);
		int Apply(Domain& domain);
		
		
		// Getters
		int SPCSetID(){return SPCSetID_;};
		int MPCSetID(){return MPCSetID_;};
		int LOADSetID(){return LOADSetID_;};
		int PLOADSetID(){return PLOADSetID_;};
		int TEMPSetID(){return TEMPSetID_;};
	
	private:
		// Variables
		int SPCSetID_;  // Global id of activated spc set
		int MPCSetID_;  // Global id of activated mpc set
		int LOADSetID_; // Global id of activated load set
		int PLOADSetID_;// Global id of activated pload set
		int TEMPSetID_; // Global id of activated temp set
	
};
#endif
