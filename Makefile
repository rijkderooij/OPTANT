# CMAKE File for "MyApp" application building against an installed Trilinos

#This file is an adaptation of the CMakeLists.txt file that was converted from 
#the buildAgainstTrilinos example. This Makefile was designed to be used in a
#flat directory structure. If you would like to run this example you will need
#put this file and src_file.cpp, src_file.hpp, main_file.cpp from
#buildAgainstTrilinos into a new directory. You will then need to set the
#environment variable MYAPP_TRILINOS_DIR to point to your base installation of
#Trilinos. Note that this example assumes that the installation of Trilinos that
#you point to has Epetra enabled.

# Get Trilinos as one entity
include ~/Libraries/trilinos-11.4.1-Source/Amesos_MPI_MUMPS/include/Makefile.export.Trilinos

# Set Arpack and Cholmod entities
ARPACKPP_INC      = -I/usr/include/arpack++
CHOLMOD_INC   	  = -I/usr/include/suitesparse
ARPACK_LIB   	  = /usr/lib/libarpack.a /usr/lib/libarpack++.a
CHOLMOD_LIB       = /usr/lib/libcholmod.a /usr/lib/libcolamd.a /usr/lib/libamd.a
ARPACK_CXXFLAG    =  -Wno-write-strings

# Additional directories to include (later used in variable INCLUDE_DIRS)

# Make sure to use same compilers and flags as Trilinos
CXX=$(Trilinos_CXX_COMPILER)
CC=$(Trilinos_C_COMPILER)
FORT=$(Trilinos_Fortran_COMPILER)

CXX_FLAGS= -std=c++0x $(Trilinos_CXX_COMPILER_FLAGS) $(USER_CXX_FLAGS) $(ARPACK_CXXFLAG)
C_FLAGS=$(Trilinos_C_COMPILER_FLAGS) $(USERC_FLAGS)
FORT_FLAGS=$(Trilinos_Fortran_COMPILER_FLAGS) $(USER_FORT_FLAGS)

INCLUDE_DIRS_T=$(Trilinos_INCLUDE_DIRS) $(ARPACKPP_INC) $(CHOLMOD_INC)
LIBRARY_DIRS_T=$(Trilinos_LIBRARY_DIRS)
LIBRARIES_T=$(Trilinos_LIBRARIES)

INCLUDE_DIRS=$(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(ARPACKPP_INC) $(CHOLMOD_INC)
LIBRARY_DIRS=$(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS)
LIBRARIES=$(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(ARPACK_LIB) $(CHOLMOD_LIB)

LINK_FLAGS=$(Trilinos_EXTRA_LD_FLAGS)

#just assuming that epetra is turned on.
DEFINES=-DMYAPP_EPETRA

# Specify number of processors
NPROC = 2

#Define file names (make sure there are no empty spaces behind the names)
MAINFILE = OPTANT
TESTFILE = OPTANT_test
ELEMTESTFILE = OPTANT_test_elem
INPUTFILE = models/PlateBucklingPSModel.dat
SOLVERFILE = SolverParameters.dat
REQUESTFILE = OutputRequests.dat

CPPFILES    = \
src/Domain.cpp src/Material.cpp src/Section.cpp src/Property.cpp \
src/PBEAM.cpp src/PSHELL.cpp \
src/Node.cpp src/Element.cpp \
src/CBEAM.cpp src/CTRIA.cpp src/CQUAD.cpp\
src/LoadCase.cpp src/LOAD.cpp src/PLOAD.cpp src/SPC.cpp src/MPC.cpp src/TEMP.cpp \
src/UtilityFunctions.cpp src/MatrixOperations.cpp src/ArpackOperations.cpp \
src/Solver.cpp src/LinearStaticSolver.cpp src/LinearBucklingSolver.cpp \
src/OutputRequest.cpp

PCH_OPTANT_FILENAME = include/PCH_OPTANT.h
PCH_OPTANT = $(PCH_OPTANT_FILENAME).gch


default: print_info compile 	  exec
etest:   print_info compile_etest exec_etest
test:    print_info compile_test  exec_test

# Echo trilinos build info just for fun
print_info:
	@echo "\nFound Trilinos!  Here are the details: "
	@echo "   Trilinos_VERSION = $(Trilinos_VERSION)"
	@echo "   Trilinos_PACKAGE_LIST = $(Trilinos_PACKAGE_LIST)"
	@echo "   Trilinos_LIBRARIES = $(Trilinos_LIBRARIES)"
	@echo "   Trilinos_INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS)"
	@echo "   Trilinos_LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS)"
	@echo "   Trilinos_TPL_LIST = $(Trilinos_TPL_LIST)"
	@echo "   Trilinos_TPL_INCLUDE_DIRS = $(Trilinos_TPL_INCLUDE_DIRS)"
	@echo "   Trilinos_TPL_LIBRARIES = $(Trilinos_TPL_LIBRARIES)"
	@echo "   Trilinos_TPL_LIBRARY_DIRS = $(Trilinos_TPL_LIBRARY_DIRS)"
	@echo "   Trilinos_BUILD_SHARED_LIBS = $(Trilinos_BUILD_SHARED_LIBS)"
	@echo "End of Trilinos details\n"

# ----------------------------------------------------------------------------------#
# -----------------------------------  MAIN FILE -----------------------------------#
# ----------------------------------------------------------------------------------#

# Compile main file	
compile: $(MAINFILE).o 
	$(CXX) $(CXX_FLAGS) -c $(PCH_OPTANT_FILENAME) -o $(PCH_OPTANT) \
	$(INCLUDE_DIRS_T) $(DEFINES) $(LIBRARY_DIRS_T) $(LIBRARIES_T)	
	$(CXX) $(CXX_FLAGS) $(MAINFILE).cpp $(CPPFILES) -o \
	$(MAINFILE) $(LINK_FLAGS) $(INCLUDE_DIRS) $(DEFINES) \
	$(LIBRARY_DIRS) $(LIBRARIES)

# build the main file
$(MAINFILE).o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) $(DEFINES) $(MAINFILE).cpp

# execute main file
exec: 
	mpirun -np $(NPROC) ./$(MAINFILE) -inp $(INPUTFILE) -sol $(SOLVERFILE) -req $(REQUESTFILE)


# ----------------------------------------------------------------------------------#
# -------------------------------- CUSTOM TEST FILE --------------------------------#
# ----------------------------------------------------------------------------------#

# Compile and execute test file
compile_test: $(TESTFILE).o 
	$(CXX) $(CXX_FLAGS) -c $(PCH_OPTANT_FILENAME) -o $(PCH_OPTANT) \
	$(INCLUDE_DIRS_T) $(DEFINES) $(LIBRARY_DIRS_T) $(LIBRARIES_T)
	$(CXX) $(CXX_FLAGS) $(TESTFILE).cpp $(CPPFILES) -o \
	$(TESTFILE) $(LINK_FLAGS) $(INCLUDE_DIRS) $(DEFINES) \
	$(LIBRARY_DIRS) $(LIBRARIES)

# build the test file
$(TESTFILE).o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) $(DEFINES) $(TESTFILE).cpp

# execute test file
exec_test: 
	mpirun -np $(NPROC) ./$(TESTFILE) -inp $(INPUTFILE) -sol $(SOLVERFILE) -req $(REQUESTFILE)


# ----------------------------------------------------------------------------------#
# -------------------------------- ELEMENT TEST FILE -------------------------------#
# ----------------------------------------------------------------------------------#

# Compile and execute element_test_file
compile_etest: $(ELEMTESTFILE).o 
	$(CXX) $(CXX_FLAGS) -c $(PCH_OPTANT_FILENAME) -o $(PCH_OPTANT) \
	$(INCLUDE_DIRS_T) $(DEFINES) $(LIBRARY_DIRS_T) $(LIBRARIES_T)	
	$(CXX) $(CXX_FLAGS) $(ELEMTESTFILE).cpp $(CPPFILES) -o \
	$(ELEMTESTFILE) $(LINK_FLAGS) $(INCLUDE_DIRS) $(DEFINES) \
	$(LIBRARY_DIRS) $(LIBRARIES)

# build the element_test file
$(ELEMTESTFILE).o:
	$(CXX) -c $(CXX_FLAGS) $(INCLUDE_DIRS) $(DEFINES) $(ELEMTESTFILE).cpp

# execute element_test file
exec_etest: 
	mpirun -np $(NPROC) ./$(ELEMTESTFILE) -inp $(INPUTFILE)



.PHONY: clean
clean:
	rm -f *.o *.a *.exe 
