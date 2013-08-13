ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)
#ROOTCFLAGS    = -pthread -m64 -I/Users/lucamartini/root/include
#ROOTGLIBS     = -L/Users/lucamartini/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/Users/lucamartini/root/lib -lm -ldl

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -O2 -Wall -Wextra 
LDFLAGS       = -g3 -lRooFit -lRooFitCore -lMinuit -lFoam -lGraf -lRooStats -lTMVA -lPostscript

SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)
#ARCH         := macosx64
#PLATFORM     := macosx

CXXFLAGS      += $(ROOTCFLAGS)
NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../../../../../lib/

###########
# TARGETS #
###########

all: main_make_pdf make_bdt_uml_inputs main_fitData main_toyMC

mass_vs_bdt: mass_vs_bdt.cc
		$(CXX) $(CXXFLAGS) -o mass_vs_bdt.o $(GLIBS)  $(LDFLAGS)  $<

main_toyMC: main_toyMC.cc pdf_analysis.o pdf_fitData.o pdf_toyMC.o
		$(CXX) $(CXXFLAGS) -o main_toyMC.o pdf_analysis.o pdf_fitData.o pdf_toyMC.o $(GLIBS)  $(LDFLAGS)  $<

main_fitData: main_fitData.cc pdf_analysis.o pdf_fitData.o
		$(CXX) $(CXXFLAGS) -o main_fitData.o pdf_analysis.o pdf_fitData.o $(GLIBS)  $(LDFLAGS)  $<

plot_arc: plot_arc.cc plot_arc.o 
		$(CXX) $(CXXFLAGS) -o plot_arc.o $(GLIBS)  $(LDFLAGS)  $<

make_bdt_uml_inputs: make_bdt_uml_inputs.cc pdf_analysis.o 
		$(CXX) $(CXXFLAGS) -o make_bdt_uml_inputs.o pdf_analysis.o $(GLIBS)  $(LDFLAGS)  $<

main_make_pdf: main_make_pdf.cc pdf_analysis.o
		$(CXX) $(CXXFLAGS) -o main_make_pdf.o pdf_analysis.o $(GLIBS)  $(LDFLAGS)  $<

pdf_analysis.o: pdf_analysis.cpp
	$(CXX) $(CXXFLAGS) -c pdf_analysis.cpp

pdf_toyMC.o: pdf_toyMC.cpp
	$(CXX) $(CXXFLAGS) -c pdf_toyMC.cpp

pdf_fitData.o: pdf_fitData.cpp
	$(CXX) $(CXXFLAGS) -c pdf_fitData.cpp

clean:
	rm -f *.o *~


