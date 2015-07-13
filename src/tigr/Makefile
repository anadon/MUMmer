#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif
ifndef AUX_BIN_DIR
AUX_BIN_DIR := $(CURDIR)
endif

OBJ_RULE = $(CXX) $(CXXFLAGS) $< -c -o $@
BIN_RULE = $(CXX) $(CXXFLAGS) $^ -o $(BIN_DIR)/$@; \
           chmod 755 $(BIN_DIR)/$@
AUX_BIN_RULE = $(CXX) $(CXXFLAGS) $^ -o $(AUX_BIN_DIR)/$@; \
           chmod 755 $(AUX_BIN_DIR)/$@
VPATH := $(AUX_BIN_DIR):$(BIN_DIR)

ALL := annotate combineMUMs delta-filter gaps mgaps \
     postnuc postpro prenuc prepro repeat-match \
     show-aligns show-coords show-tiling show-snps \
     show-diff show-variants


#-- PHONY rules --#
.PHONY: all clean


all: $(ALL)


clean:
	rm -f *.o *~
	cd $(BIN_DIR); rm -f $(ALL)
	cd $(AUX_BIN_DIR); rm -f $(ALL)




#-- not so PHONY rules --#
delta.o: delta.cpp delta.h
	$(OBJ_RULE)

tigrinc.o: tigrinc.cpp tigrinc.h
	$(OBJ_RULE)

sw_align.o: sw_align.cpp sw_align.h tigrinc.h
	$(OBJ_RULE)

translate.o: translate.cpp translate.h
	$(OBJ_RULE)


annotate: annotate.cpp tigrinc.o
	$(BIN_RULE)

combineMUMs: combineMUMs.cpp tigrinc.o
	$(BIN_RULE)

delta-filter: delta-filter.cpp tigrinc.o delta.o
	$(BIN_RULE)

gaps: gaps.cpp tigrinc.o
	$(BIN_RULE)

mgaps: mgaps.cpp tigrinc.o
	$(BIN_RULE)

postnuc: postnuc.cpp tigrinc.o sw_align.o
	$(AUX_BIN_RULE)

postpro: postpro.cpp tigrinc.o sw_align.o translate.o
	$(AUX_BIN_RULE)

prenuc: prenuc.cpp tigrinc.o
	$(AUX_BIN_RULE)

prepro: prepro.cpp tigrinc.o translate.o
	$(AUX_BIN_RULE)

repeat-match: repeat-match.cpp tigrinc.o
	$(BIN_RULE)

show-aligns: show-aligns.cpp tigrinc.o translate.o delta.o
	$(BIN_RULE)

show-coords: show-coords.cpp tigrinc.o delta.o
	$(BIN_RULE)

show-tiling: show-tiling.cpp tigrinc.o delta.o
	$(BIN_RULE)

show-snps: show-snps.cpp tigrinc.o translate.o delta.o
	$(BIN_RULE)

show-variants: show-variants.cpp tigrinc.o translate.o delta.o
	$(BIN_RULE)

show-diff: show-diff.cpp tigrinc.o delta.o
	$(BIN_RULE)
