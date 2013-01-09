ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

CXXFLAGS      += $(ROOTCFLAGS)

NGLIBS         = $(ROOTGLIBS) 
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/cms/cmssw/CMSSW_3_6_3/external/slc5_ia32_gcc434/lib -lMinuit -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)
# Hack here, because RooFit is removed from ROOT:
CXX	         += -I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms6/include
OUTLIB	         = ../../../../lib/slc4_ia32_gcc345/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../../../../lib/slc4_ia32_gcc345/

MergeDataSets:  $(INCLUDEDIR)MergeDataSets.C  
	$(CXX) $(CXXFLAGS) -o MergeDataSets $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f $(OUTLIB)*.so
	rm -f MergeDataSets

all: MergeDataSets
