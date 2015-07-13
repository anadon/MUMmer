#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

##include ../Makedef

LIBBASEDIR=../libbasedir
INCLUDEDIR=../libbasedir

# the following flags determines the maximal allowed input size of
# the suffix tree construction
# STREESMALL means that the maximal sequence length is 2 MB
# STREELARGE means that the maximal sequence length is 134 MB
# STREEHUGE  means that the maximal sequence length is 500 MB

#SIZEFLAG=-DSTREESMALL
#SIZEFLAG=-DSTREELARGE
SIZEFLAG=-DSTREEHUGE

CFLAGS+=-I$(INCLUDEDIR) $(SIZEFLAG)

LD=${CC}

LIBBASE=${LIBBASEDIR}/libbase.a
LIBBASEDBG=${LIBBASEDIR}/libbase.dbg.a

PROTOFILES=access.cpp scanpref.cpp linkloc.cpp depthtab.cpp\
           ex2leav.cpp dfs.cpp overmax.cpp overs.cpp\
           addleafcount.cpp iterator.cpp 

OBJECTS=construct.o access.o scanpref.o linkloc.o depthtab.o\
          ex2leav.o dfs.o overmax.o overs.cpp.o\
          addleafcount.o iterator.o


all: libs progs

libs: libstree.a

progs: stree.x loc.x

libstree.a: $(OBJECTS)
	ar sruv $@ $(OBJECTS)

stree.x: stree.o libstree.a
	$(LD) $(LDFLAGS) stree.o libstree.a $(LIBBASE) -o $@

loc.x: loc.o libstree.a
	$(LD) $(LDFLAGS) loc.o libstree.a $(LIBBASE) -o $@

clean:
	rm -f *.[aox] *~
