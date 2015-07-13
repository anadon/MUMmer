#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS


EXECS=mummer maxmat3.x

LIBBASEDIR=../shared/
LIBSTREEDIR=../streesrc/
INCLUDEDIR= -I${LIBBASEDIR} -I${LIBSTREEDIR}

CFLAGS= -c -O2 -Wall

KURTZLIBS=$(LIBSTREEDIR)/libstree.a $(LIBBASEDIR)/libbase.a

MUM3OBJECTS=maxmat3.o maxmatopt.o maxmatinp.o findmumcand.o \
																							findmaxmat.o procmaxmat.o
																							
SOURCE=findmaxmat.cpp maxmat3.cpp maxmatopt.cpp findmumcand.cpp \
																						maxmatinp.cpp procmaxmat.cpp


all : $(MUM3OBJECTS) $(KURTZLIBS)
	$(CXX) $(INCLUDEDIR) $(LDFLAGS) $(MUM3OBJECTS) $(KURTZLIBS) -o $@
	cp $@ maxmat3.x

%.o : %.cpp $(SOURCE)
	$(CXX) $(CFLAGS) $(INCLUDEDIR)  $< -o $@

clean:
	rm -f *.[oa] *~
	cd $(BIN_DIR); rm -f $(ALL) mummer

