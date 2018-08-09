#
# General Makefile for the OfflineUser package
#
#

# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
TOPDIR= $(shell pwd)

BINDIR  = $(TOPDIR)/bin
LIBDIR  = $(TOPDIR)/lib
SRCDIR  = $(TOPDIR)/src
INCDIR  = $(TOPDIR)/include
OBJDIR  = $(TOPDIR)/obj

USER_SRCS = $(wildcard $(SRCDIR)/*.cc)
files := $(foreach USER_SRCS,$(SRCDIR),$(USER_SRCS))

HEADERS_DICT = $(INCDIR)/DataReader.h $(INCDIR)/PropagationMCReader.h $(INCDIR)/SourceFitter.h $(INCDIR)/DrawFitResults.h $(INCDIR)/ConfigParser.h 

OBJS = $(USER_SRCS:.cc=.o)

MAIN_CRSOURCEFITTER = CRSourceFitter.cc
MAIN_OBJ = $(MAIN_CRSOURCEFITTER:.cc=.o)

## Get platform 32/64 bit
LBITS   = $(shell getconf LONG_BIT)

# Set executable a name
SHARED_LIB = libCRSourceFitter.so
EXE = CRSourceFitter
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean


######################################
###  CPPFLAGS & CXXFLAGS  ############
######################################

CPPFLAGS = -I$(INCDIR)

ifeq ($(LBITS),64)
  # do 64 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m64
	CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64 
	SYSLIBDIR = /usr/lib/x86_64-linux-gnu
else
  # do 32 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m32
  CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m32 
	SYSLIBDIR = /usr/lib
endif

SOFLAGS = -fPIC -ggdb3 -Wall -shared

CPPFLAGS_ROOT= -I$(ROOTSYS)/include
CPPFLAGS += $(CPPFLAGS_ROOT) 

###########################
###  LDFLAGS   ############
###########################
LDFLAGS_ROOT = $(shell root-config --libs) -L$(ROOTSYS)/lib -lSpectrum -lMathMore -lMinuit
LDFLAGS_SYSTEM = -L$(SYSLIBDIR) -lrt 
LDFLAGS = $(LDFLAGS_ROOT) 

################################################################


#all: GETOBJS PRINTINFO $(SHARED_LIB) PUTOBJS $(EXE) PUTBINARIES PUTLIBRARIES
all: GETOBJS PRINTINFO $(SHARED_LIB) $(EXE) PUTOBJS

PRINTINFO: 
	@echo 'Compiling $(EXE) on a $(LBITS) bit machine' \

GETOBJS:
	@echo "Put object and lib files again to $(SRCDIR) dir"
	- mv -f $(OBJDIR)/*.o $(SRCDIR) 2>/dev/null
	- mv -f $(LIBDIR)/$(SHARED_LIB) $(TOPDIR) 2>/dev/null
	- mv -f $(LIBDIR)/ClassDictionary_rdict.pcm $(TOPDIR) 2>/dev/null
	- mv -f $(BINDIR)/$(EXE) $(TOPDIR) 2>/dev/null

PUTOBJS:
	@echo "Moving object files to $(OBJDIR) dir, binariers to $(BINDIR) and libraries to $(LIB) dir"
	- mv -f $(SRCDIR)/*.o $(OBJDIR) 2>/dev/null
	- mv -f *.o $(OBJDIR) 2>/dev/null
	- mv -f $(SHARED_LIB) $(LIBDIR) 2>/dev/null
	- mv -f ClassDictionary_rdict.pcm $(LIBDIR) 2>/dev/null
	- mv -f $(EXE) $(BINDIR) 2>/dev/null

PUTBINARIES:
	@echo "Moving binary files to $(BINDIR) dir"
	- mv -f $(EXE) $(BINDIR) 2>/dev/null

PUTLIBRARIES:
	@echo "Moving library files to $(LIBDIR) dir"
	- mv -f $(SHARED_LIB) $(LIBDIR) 2>/dev/null

ClassDictionary.cc: $(HEADERS_DICT) LinkDef.h
	rootcint -f $@ -c -p $(CPPFLAGS) $^
	
$(SHARED_LIB): $(OBJS) ClassDictionary.o
	@echo "Compiler is $(CXX) or $(CC), options are $(CXXFLAGS), generating $@ shared library..."
	@$(CXX) $(CXXFLAGS) $(SOFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) 
	

$(EXE): $(MAIN_CRSOURCEFITTER_OBJ)
	@echo "Building CRSourceFitter ..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(MAIN_CRSOURCEFITTER) $(LDFLAGS) -L$(TOPDIR) -lCRSourceFitter
#@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(MAIN_CRSOURCEFITTER) $(LDFLAGS) -L$(LIBDIR) -lCRSourceFitter


#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $^ > $@

clean:
	- rm -f *.o $(OBJDIR)/*.o $(SRCDIR)/*.o
	- rm -f $(EXE) $(BINDIR)/$(EXE)
	- rm -f $(SHARED_LIB) $(LIBDIR)/$(SHARED_LIB)
	- rm -f ClassDictionary.cc ClassDictionary.h ClassDictionary_rdict.pcm $(LIBDIR)/ClassDictionary_rdict.pcm
	- rm -f core Make-depend

#############################################################
# 'make run' will run the thing



#############################################################
# the lines below are for running with debugger 'make run_gdb'

.INTERMEDIATE: gdb.cmdl

# batch mode gdb needs a file with commands
gdb.cmdl:
	echo "r -b bootstrap.xml" > $@

run_gdb: gdb.cmdl $(EXE)
	gdb -batch -x $< ./$(EXE) && touch $@

#include Make-depend
