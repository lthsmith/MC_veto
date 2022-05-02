# Makefile for dShower

# Library for gfortran
LDLIBS=-lgfortran

# Compiler flags:
CXXFLAGS= -Wall -W
CCFLAGS= -Wall

# Directory containing the .o files
ODIR = ../tmp

# External flags, if necessary:
# Sherpa
ifdef DSHOWER_WITH_SHERPA
LIBEXT+=-L$(SHERPA)/lib/SHERPA-MC -lSherpaMain -lMEProcess -lToolsPhys -lToolsMath -lToolsOrg -lPhasicProcess -lSherpaPerturbativePhysics -lPhasicMain -lModelMain
EXTFLAGS+=-I$(SHERPA)/include/SHERPA-MC 
DEFINENAME+=-DDSHOWER_WITH_SHERPA
endif

# OpenLoops
ifdef DSHOWER_WITH_OL
LIBEXT+=-L$(OPENLOOPS)/lib -lopenloops -lcollier -lcuttools -lolcommon -loneloop -lrambo -ltrred
EXTFLAGS+=-I$(OPENLOOPS)/include
DEFINENAME+=-DDSHOWER_WITH_OL
endif

# ChiliPDF
ifdef DSHOWER_WITH_CHILIPDF
LIBEXT+=-L$(CHILIPDF)/build/lib -lchilipdf -lxmath
EXTFLAGS+=-I$(CHILIPDF)/include -I$(CHILIPDF)/xmath/include -I$(CHILIPDF)/build/include
DEFINENAME+=-DDSHOWER_WITH_CHILIPDF
endif

# ChiliPDF requires c++14
ifdef DSHOWER_WITH_CHILIPDF
CXXFLAGS += -std=c++14
else
CXXFLAGS += -std=c++11
endif

EXEC=main main_check main_check_luminosity MC_MT_veto Integrand_scan

# Only include main_chilipdf if ChiliPDF is used
ifdef DSHOWER_WITH_CHILIPDF
EXEC += main_chilipdf
endif

all: $(EXEC)

#####################################
# Fortran:

$(ODIR)/gsalps.o : gsalps.f
	$(FC) -c $(FFLAGS) $< -o $@

$(ODIR)/gsdpdf.o : gsdpdf.f
	$(FC) -c $(FFLAGS) $< -o $@

# Optional fortran files

$(ODIR)/gsdpdfYinter.o : gsdpdfYinter.f
	$(FC) -c $(FFLAGS) $< -o $@

$(ODIR)/intrisplit.o : intrisplit.f
	$(FC) -c $(FFLAGS) $< -o $@

$(ODIR)/gsdpdfYinter_c.o : gsdpdfYinter_c.f
	$(FC) -c $(FFLAGS) $< -o $@

#####################################
# C:

$(ODIR)/exponential_integral_Ei.o : exponential_integral_Ei.c exponential_integral_Ei.h
	$(CC) -c $(CCFLAGS) $< -o $@


#####################################
# C++:

# Note: gsdpdfYinter_c.o for y-dependent DGS set and gsdpdf.o for GS09 set. Remove .o before compiling.

objects = $(ODIR)/gsalps.o $(ODIR)/gsdpdf.o 
objects += $(ODIR)/Basics.o $(ODIR)/ColourFlow.o $(ODIR)/DGS.o $(ODIR)/dShower.o $(ODIR)/Event.o $(ODIR)/General.o $(ODIR)/HardProcess.o $(ODIR)/Kinematics.o $(ODIR)/MECorr.o $(ODIR)/mstwpdf.o $(ODIR)/PDF.o $(ODIR)/ProtoDGS.o $(ODIR)/Shower.o $(ODIR)/XSection.o
objects += $(ODIR)/exponential_integral_Ei.o

ifdef DSHOWER_WITH_OL
objects += $(ODIR)/External.o
endif

ifdef DSHOWER_WITH_SHERPA
objects += $(ODIR)/External.o
endif

main : $(ODIR)/main.o $(objects)
	$(CXX) $^ -o dShower $(LDLIBS) $(LIBEXT)

main_check : $(ODIR)/main_check.o $(objects)
	$(CXX) $^ -o dShower_check $(LDLIBS) $(LIBEXT)
	
main_check_luminosity : $(ODIR)/main_check_luminosity.o $(objects)
	$(CXX) $^ -o dShower_check_luminosity $(LDLIBS) $(LIBEXT)
	
MC_MT_veto : $(ODIR)/MC_MT_veto.o $(objects)
	$(CXX) $^ -o MC_MT_veto $(LDLIBS) $(LIBEXT)
	
Integrand_scan : $(ODIR)/Integrand_scan.o $(objects)
	$(CXX) $^ -o Integrand_scan $(LDLIBS) $(LIBEXT)

main_chilipdf : $(ODIR)/main_chilipdf.o $(objects)
ifdef DSHOWER_WITH_CHILIPDF
	$(CXX) $^ -o dShower_chilipdf $(LDLIBS) $(LIBEXT)
else
	@echo "Error: $@ requires ChiliPDF"
endif

# Create .o files.

$(ODIR)/main.o : main.cc dShower.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/main_check.o : main_check.cc dShower.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/main_check_luminosity.o : main_check_luminosity.cc dShower.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/MC_MT_veto.o : MC_MT_veto.cc dShower.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@
	
$(ODIR)/Integrand_scan.o : Integrand_scan.cc dShower.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/main_chilipdf.o : main_chilipdf.cc dShower.h
ifdef DSHOWER_WITH_CHILIPDF
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $(DEFINENAME) $< -o $@
else
	@echo "Error: $@ requires ChiliPDF"
endif

$(ODIR)/Basics.o : Basics.cc Basics.h dShowerStdlib.h 
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/ColourFlow.o : ColourFlow.cc ColourFlow.h dShowerStdlib.h Event.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/DGS.o : DGS.cc DGS.h dShowerStdlib.h 
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/dShower.o : dShower.cc dShower.h dShowerStdlib.h cxxopts.hpp INIReader.h Basics.h General.h Event.h Shower.h HardProcess.h XSection.h MECorr.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/Event.o : Event.cc Event.h dShowerStdlib.h General.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/General.o : General.cc General.h dShowerStdlib.h LHEF.h Basics.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/HardProcess.o : HardProcess.cc HardProcess.h dShowerStdlib.h Basics.h General.h XSection.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $< -o $@

$(ODIR)/Shower.o : Shower.cc Shower.h dShowerStdlib.h Basics.h General.h Event.h PDF.h Kinematics.h ColourFlow.h MECorr.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/Kinematics.o : Kinematics.cc Kinematics.h dShowerStdlib.h Basics.h General.h Event.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/MECorr.o : MECorr.cc MECorr.h dShowerStdlib.h Basics.h General.h PDF.h Event.h ColourFlow.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/mstwpdf.o : mstwpdf.cc mstwpdf.h dShowerStdlib.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/PDF.o : PDF.cc PDF.h dShowerStdlib.h General.h mstwpdf.h DGS.h ProtoDGS.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $(DEFINENAME) $< -o $@

$(ODIR)/ProtoDGS.o : ProtoDGS.cc ProtoDGS.h dShowerStdlib.h General.h mstwpdf.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/XSection.o : XSection.cc XSection.h dShowerStdlib.h Basics.h General.h PDF.h ColourFlow.h External.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $(DEFINENAME) $< -o $@

$(ODIR)/External.o : External.cc External.h dShowerStdlib.h Basics.h General.h
	$(CXX) -c $(CXXFLAGS) $(EXTFLAGS) $(DEFINENAME) $< -o $@

#####################################
# Install:

install :
	mv dShower ../bin/
	mv dShower_check ../bin/
ifdef DSHOWER_WITH_CHILIPDF
	mv dShower_chilipdf ../bin/
endif

#####################################

clean :
	-rm $(ODIR)/*.o
