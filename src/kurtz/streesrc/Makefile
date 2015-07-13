#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS


INCLUDEDIR=../shared

# STREEHUGE  means that the maximal sequence length is 500 MB
SIZEFLAG=-DSTREEHUGE

CXXFLAGS+=-I$(INCLUDEDIR) $(SIZEFLAG) -c
LDFLAGS=-I$(INCLUDEDIR) $(SIZEFLAG)

CXX=g++
LD=g++

LIBBASE=${LIBBASEDIR}/libbase.a

SOURCE=access.cpp construct.cpp scanpref.cpp linkloc.cpp depthtab.cpp ex2leav.cpp \
						dfs.cpp overmax.cpp oversucc.cpp addleafcount.cpp iterator.cpp 

OBJECTS=access.o construct.o scanpref.o linkloc.o depthtab.o ex2leav.o \
									dfs.o overmax.o oversucc.o addleafcount.o iterator.o
									
##objects and source that conflict with each other or are not suitable  
##for inclusion into a Libarary
ESOURCE=stree.spp

EOBJECTS=stree.o

EXECS=stree loc

OUTPUTLIBS=libstree.a

all: $(OUTPUTLIBS) $(EXECS)

$(EXECS): $(OUTPUTLIBS) $(OBJECTS) $(EOBJECTS)
	$(LD) $(OBJECTS) $(EOBJECTS) $(OUTPUTLIBS) $(LDFLAGS) -o $@

$(OUTPUTLIBS) : $(OBJECTS)
	ar ruv $@ $(OBJECTS)

%.o : %.cpp $(EXTERNLIBS)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECS) $(OUTPUTLIBS)
