// Include own classes
#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 
#include "fei_fstream.hpp"
#include "fei_base.hpp"


// Functions
void OutputRequest::SetParameters(fei::ParameterSet& params){
	params.getBoolParamValue("PrintModelInfo", PrintModelInfo_);
	params.getBoolParamValue("PrintNodalOutput", PrintNodalOutput_);                 numNodalOutputs_++;
	params.getBoolParamValue("PrintElementOutput", PrintElementOutput_);
	params.getBoolParamValue("PrintElementStrains", PrintElementStrains_);           numElementOutputs_++;
	params.getBoolParamValue("PrintElementForces", PrintElementForces_);             numElementOutputs_++;
	params.getBoolParamValue("PrintElementStrainEnergy", PrintElementStrainEnergy_); numElementOutputs_++;
	
	params.getBoolParamValue("PrintNodalBucklingOutput", PrintNodalBucklingOutput_);
}

// General write functions which do not belong to particular class
void OutputRequest::PrintModelInfo(Domain &domain, std::ofstream& fout){
	
	if(PrintModelInfo_){
		
	// Write content of MaterialList
	fout << FEI_ENDL << FEI_ENDL << "Elements in MaterialList " << FEI_ENDL;
	for (std::map<const int, Material*>::iterator it = domain.MaterialList.begin();
		it != domain.MaterialList.end();
		++it)
		fout << "  [" << (*it).first << ", " << (*it).second << "]" << FEI_ENDL;
	
	std::map<const int, Material*>::iterator itr = domain.MaterialList.begin();
	fout <<  FEI_ENDL << "Properties of first material: " << FEI_ENDL; 
	fout << "Name: " << itr->second->MaterialName() << FEI_ENDL;
	fout << "PlaneStressStiffness: "; 
	UtilityFunctions::PrintArray(itr->second->PlaneStressStiffness(),9,fout);
	fout << "ThermalExpansion:     "; 
	UtilityFunctions::PrintArray(itr->second->ThermalExpansion(),6,fout);
	
	if(domain.MaterialList.size()>1){
		++itr;
		fout <<  FEI_ENDL << "Properties of second material: " << FEI_ENDL; 
		fout << "Name: " << itr->second->MaterialName() << FEI_ENDL;
		fout << "PlaneStressStiffness: "; 
		UtilityFunctions::PrintArray(itr->second->PlaneStressStiffness(),9,fout);
		fout << "ThermalExpansion:     "; 
		UtilityFunctions::PrintArray(itr->second->ThermalExpansion(),6,fout);
	}
	
	
	// Write content of PropertyList
	fout << FEI_ENDL << FEI_ENDL << "Elements in PropertyList " << FEI_ENDL;
	for (std::vector<Property*>::iterator it = domain.PropertyList.begin();
		it != domain.PropertyList.end();
		++it)
		fout <<  *it << FEI_ENDL;
	
	std::vector<Property*>::iterator itrP = domain.PropertyList.begin();
	fout <<  FEI_ENDL << "First property: " << FEI_ENDL; 
	fout << "Global ID: " << domain.PropertyID_LG[0] << FEI_ENDL;
	fout << "Local ID:  " << domain.PropertyID_GL[domain.PropertyID_LG[0]] << FEI_ENDL;
	fout << "A (packed):        " << FEI_ENDL;
	UtilityFunctions::PrintMatrix(static_cast<PSHELL*>(*itrP)->A(),1,6,fout);
	fout << "B (packed):        " << FEI_ENDL;
	UtilityFunctions::PrintMatrix(static_cast<PSHELL*>(*itrP)->B(),1,6,fout);
	fout << "D (packed):        " << FEI_ENDL;
	UtilityFunctions::PrintMatrix(static_cast<PSHELL*>(*itrP)->D(),1,6,fout);
	

	
	// Write content of NodeList
	fout << FEI_ENDL << FEI_ENDL << "Elements in NodeList " << FEI_ENDL;
	for (std::vector<Node*>::iterator it = domain.NodeList.begin();
		it != domain.NodeList.end();
		++it)
		fout <<  *it << FEI_ENDL;
		 
	fout << "NumNodes from vectorsize:      " << domain.NodeList.size() << FEI_ENDL;
	fout << "NumNodes from domain variable: " << domain.NumTotalNodes << FEI_ENDL;
	
	std::vector<Node*>::iterator itrN = domain.NodeList.begin(); ++itrN;
	fout <<  FEI_ENDL << "Second node: " << FEI_ENDL; 
	fout << "Global ID:             " << domain.NodeID_LG[1] << FEI_ENDL;
	fout << "Local ID (NodeID_GL):  " << domain.NodeID_GL[domain.NodeID_LG[1]] << FEI_ENDL;
	fout << "Local ID (getter):     " << domain.NodeList[1]->LocalNodeID() << FEI_ENDL;
	fout << "Number of DOF:     	" << domain.NodalNumDOF[1] << FEI_ENDL;
	fout << "GlobalCoordinates:   "; 
	UtilityFunctions::PrintArray((*itrN)->GlobalCoordinates(),3,fout);
	
	
	// Write content of ElementList
	fout << FEI_ENDL << FEI_ENDL << "Elements in ElementList " << FEI_ENDL;
	for (std::vector<Element*>::iterator it = domain.ElementList.begin();
		it != domain.ElementList.end();
		++it)
		fout <<  *it << FEI_ENDL;
	 
	fout << "NumElements from vectorsize:      " << domain.ElementList.size() << FEI_ENDL;
	fout << "NumElements from domain variable: " << domain.NumTotalElements << FEI_ENDL;
	
	std::vector<Element*>::iterator itrE = domain.ElementList.begin(); ++itrE;
	fout <<  FEI_ENDL << "Second element: " << FEI_ENDL; 
	fout << "Global ID:  " << domain.ElementID_LG[1] << FEI_ENDL;
	fout << "Local ID (ElementID_GL):  " << domain.ElementID_GL[domain.ElementID_LG[1]] << FEI_ENDL;
	fout << "Local ID (getter):        " << (*itrE)->LocalElementID() << FEI_ENDL;
	fout << "MaterialDirection:        " << FEI_ENDL;
	UtilityFunctions::PrintMatrix(static_cast<CTRIA*>(*itrE)->MaterialDirection(),1,3,fout);
	fout << "NumProperties:   " << (*itrE)->NumProperties() << FEI_ENDL;
	for(int i=0; i<(*itrE)->NumProperties(); i++){
		fout << " Prop " << i << ": " << (*itrE)->Properties()[i] <<FEI_ENDL;
	}
	fout << "Nodes:   "<< (*itrE)->Nodes()[0] << "  " 
		 << (*itrE)->Nodes()[1]<< "  "  << (*itrE)->Nodes()[2] << FEI_ENDL;
	
	// Write contents of NodesPerElementBlock, NumElementsBlock, 
	// and ElementBlockID
	fout << FEI_ENDL << FEI_ENDL << "NodesPerElementBlock: " << FEI_ENDL;
	for (std::vector<int>::iterator it = domain.NodesPerElementBlock.begin();
		it != domain.NodesPerElementBlock.end();
		++it)
		fout <<  *it << "  ";
	fout << FEI_ENDL << FEI_ENDL << "NumElementsBlock: " << FEI_ENDL;
	for (std::vector<int>::iterator it = domain.NumElementsBlock.begin();
		it != domain.NumElementsBlock.end();
		++it)
		fout <<  *it << "  ";
	fout << FEI_ENDL << FEI_ENDL << "ElementBlockID: " << FEI_ENDL;
	for (std::vector<int>::iterator it = domain.ElementBlockID.begin();
		it != domain.ElementBlockID.end();
		++it)
		fout <<  *it << "  ";
	
	// Write content of SPCList
	fout << FEI_ENDL << FEI_ENDL << "Elements in SPCList " << FEI_ENDL;
	for (std::multimap<const int, SPC*>::iterator it = domain.SPCList.begin();
		it != domain.SPCList.end();
		++it)
		fout << "  [" << (*it).first << ", " << (*it).second << "]" << FEI_ENDL;
	
	std::multimap<const int, SPC*>::iterator itrS = domain.SPCList.begin();
	fout <<  FEI_ENDL << "Properties of first SPC: " << FEI_ENDL; 
	fout << "SPCNode: " << itrS->second->GetNode() << FEI_ENDL;
	fout << "PrescribedValue: " << itrS->second->PrescribedValue() << FEI_ENDL;
	fout << "PrescribedDOF: "; 
	UtilityFunctions::PrintArray(itrS->second->PrescribedDOF(),6,fout);
	
	fout <<  FEI_ENDL << "The SPC in set " << itrS->first << " are: "; 	
	std::multimap<const int, SPC*>::iterator itrS2;
	for (itrS2 = domain.SPCList.equal_range(itrS->first).first; 
			itrS2!=domain.SPCList.equal_range(itrS->first).second; ++itrS2)
		fout <<  " " << (*itrS2).second ;
	fout<< FEI_ENDL;
	
	// Write content of MPCList
	fout << FEI_ENDL << FEI_ENDL << "Elements in MPCList " << FEI_ENDL;
	for (std::multimap<const int, MPC*>::iterator it = domain.MPCList.begin();
		it != domain.MPCList.end();
		++it)
		fout << "  [" << (*it).first << ", " << (*it).second << "]" << FEI_ENDL;
	
	std::multimap<const int, MPC*>::iterator itrM = domain.MPCList.begin();
	fout <<  FEI_ENDL << "Properties of first MPC: " << FEI_ENDL;
	fout << "numNodes: " << itrM->second->NumNodes() << FEI_ENDL; 
	fout << "PrescribedDOF: "; UtilityFunctions::PrintArray(itrM->second->PrescribedDOF(),itrM->second->NumNodes(),fout);
	fout << "Weights:       "; UtilityFunctions::PrintArray(itrM->second->Weights(),itrM->second->NumNodes(),fout);
	fout << "RHSConstant:   " << itrM->second->RHSConstant() << FEI_ENDL;
	fout << "MPCNodes:      ";
			
	for(int i=0; i<itrM->second->NumNodes(); i++){
		fout << domain.NodeID_LG[itrM->second->MPCNodes()[i]->LocalNodeID()] << "  ";
	}
	fout << FEI_ENDL;
	
	fout <<  FEI_ENDL << "The MPC in set " << itrM->first << " are: "; 	
	std::multimap<const int, MPC*>::iterator itrM2;
	for (itrM2 = domain.MPCList.equal_range(itrM->first).first; 
			itrM2!=domain.MPCList.equal_range(itrM->first).second; ++itrM2)
		fout <<  " " << (*itrM2).second ;
	fout<< FEI_ENDL;
	
	// Write content of LOADList
	fout << FEI_ENDL << FEI_ENDL << "Elements in LOADList " << FEI_ENDL;
	for (std::multimap<const int, LOAD*>::iterator it = domain.LOADList.begin();
		it != domain.LOADList.end();
		++it)
		fout << "  [" << (*it).first << ", " << (*it).second << "]" << FEI_ENDL;
	
	std::multimap<const int, LOAD*>::iterator itrL = domain.LOADList.begin();
	fout <<  FEI_ENDL << "Properties of first LOAD: " << FEI_ENDL; 
	fout << "Node: " << itrL->second->GetNode() << FEI_ENDL;
	fout << "LoadVector: ";
	UtilityFunctions::PrintArray(itrL->second->LoadVector(),6,fout);
	
	fout <<  FEI_ENDL << "The LOADs in set " << itrL->first << " are: "; 	
	std::multimap<const int, LOAD*>::iterator itrL2;
	for (itrL2 = domain.LOADList.equal_range(itrL->first).first; 
			itrL2!=domain.LOADList.equal_range(itrL->first).second; ++itrL2)
		fout <<  " " << (*itrL2).second ;
	fout<< FEI_ENDL;
	
	// Write content of LoadCaseList
	fout << FEI_ENDL << FEI_ENDL << "Elements in LoadCaseList " << FEI_ENDL;
	for (std::map<const int, LoadCase*>::iterator it = domain.LoadCaseList.begin();
		it != domain.LoadCaseList.end();
		++it)
		fout << (*it).first << "    " <<  (*it).second << FEI_ENDL;
		 
	fout << "NumLoadCases from vectorsize:      " << domain.LoadCaseList.size() << FEI_ENDL;
	
	std::map<const int, LoadCase*>::iterator itrLC = domain.LoadCaseList.begin();
	fout <<  FEI_ENDL << "First LoadCase: " << FEI_ENDL; 
	fout << "Global ID:  " << (*itrLC).first << FEI_ENDL;
	fout << "SPCset:     " << (*itrLC).second->SPCSetID() << FEI_ENDL;
	fout << "MPCset:     " << (*itrLC).second->MPCSetID() << FEI_ENDL;
	fout << "LOADset:    " << (*itrLC).second->LOADSetID() << FEI_ENDL;
	fout << "PLOADset:   " << (*itrLC).second->PLOADSetID() << FEI_ENDL;
	fout << "TEMPset:    " << (*itrLC).second->TEMPSetID() << FEI_ENDL;

	
	}
}



