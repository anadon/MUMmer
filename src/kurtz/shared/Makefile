#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS


SOURCE=clock.cpp safescpy.cpp cleanMUMcand.cpp mapfile.cpp \
				multiseq.cpp procopt.cpp  seterror.cpp

OBJECTS=clock.o safescpy.o cleanMUMcand.o mapfile.o \
				multiseq.o procopt.o  seterror.o



CXXFLAGS=-Wall

LIBBASE=libbase.a


all: $(LIBBASE)

*.o: $(SOURCE)
	$(CXX) $(CXXFLAGS) -c -o $@


$(LIBBASE): $(OBJECTS)
	ar ruv $@ $(LIBOBJECTS)


clean:
	rm -f *.[oa]
