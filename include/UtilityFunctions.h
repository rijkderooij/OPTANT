#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_fwd.hpp>
#include "fei_iostream.hpp"
#include "fei_fstream.hpp"

#include <string>

//Forward declarations
class Domain;

namespace UtilityFunctions {
	
	// Get filename given in execute command
	int getFileName(int argc, char** argv, const char* argTag, char* destination);
	
	// Set of functions to read input file
	int read_input_file(const char* filename, Domain &domain);
	int ReadSolutionType(std::ifstream& fin, Domain &domain);
	int ReadMaterials(std::ifstream& fin, Domain &domain);
	int ReadProperties(std::ifstream& fin, Domain &domain);
	int ReadNodes(std::ifstream& fin, Domain &domain);
	int ReadElements(std::ifstream& fin, Domain &domain);
	int ReadSPC(std::ifstream& fin, Domain &domain);
	int ReadMPC(std::ifstream& fin, Domain &domain);
	int ReadLOAD(std::ifstream& fin, Domain &domain);
	int ReadPLOAD(std::ifstream& fin, Domain &domain);
	int ReadTEMP(std::ifstream& fin, Domain &domain);
	int ReadLoadCases(std::ifstream& fin, Domain &domain);
	
	// Read file with solver parameter
	int read_param_file(const char* filename, MPI_Comm comm,
							std::vector<std::string>& file_contents);

	
	// Get id of argument in argc, argv
	int whichArg(int argc, const char*const* argv, const char* findarg);
	
	// Read file lines, containing a key and value, into a string vector;
	void read_file_lines_into_strings(const char* filename,
							std::vector<std::string>& file_contents);
	
	// Print array
	void PrintArray(double* array, int size);
	void PrintArray(int* array, int size);
	void PrintArray(double* array, int size,std::ofstream& fout);
	void PrintArray(int* array, int size,std::ofstream& fout);
	
	void PrintMatrix(double* a, int nrow,int ncol);
	void PrintMatrix(double* a, int nrow,int ncol,std::ofstream& fout);
	
	
	
	
};
#endif
