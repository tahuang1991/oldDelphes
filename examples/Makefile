# -*- Makefile -*- 

#
# This Makefile is intended for compiling cc under delphes
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#Tao Huang, email: taohuang@physics.tamu.edu
include Makefile.arch
#LDLIBS=`root-config --glibs`
#ROOTLIBS='-lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory'
DelphesDir=/home/taohuang/Herwig++/Delphes-3.2.0
ROOTINCLUDE=$(shell root-config --glibs 2> /dev/null)
#DelphesINCLUDE=-I$(DelphesDir)/classes -I$(DelphesDir)/modules -I$(DelphesDir)/readers
DelphesINCLUDE=-I$(DelphesDir) -I$(DelphesDir)/external 
DelphesLIB=-L$(DelphesDir) -lDelphes
#Other_INCLUDE_if_necessary= 
LDFLAGS=-Wl,--no-as-needed
LDFLAGS+=`root-config --ldflags`
INCLUDE = $(DelphesINCLUDE) $(ROOTINCLUDE)
LIB=$(DelphesLIB)
#
# C++ flags
# 
CXX=g++
CXXFLAGS=-pedantic -ansi -Wall -Wno-long-long -Wno-format -Werror=uninitialized -Werror=delete-non-virtual-dtor -O2
CXXFLAGS+=`root-config --cflags`
#src file
Example_SRC = Example1.cpp Example1.C
Example_EXE = Example1.exe
# test print out
#var=TEST
#PHONY:
#	@echo $(INCLUDE)
	#echo $(LIB)
#default : DiHiggs
Example1:
	$(CXX) $(CXXFLAGS) -g $(Example_SRC) $(LDFLAGS) $(INCLUDE) $(LIB) -o $(Example_EXE)	
#project : DiHiggsWWBBAna.cc 
#	$(CXX) $(CXXFLAGS) -c $<
#	$(CXX) $(LDFLAGS) DiHiggsWWBBAna.o $(CXXOBJS) $(LIBS) -o DiHiggsWWBBAna

clean:
	rm -f *.exe
