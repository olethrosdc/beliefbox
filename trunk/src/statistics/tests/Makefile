SMPL_DIR := $(shell cd ../..; pwd)

include $(SMPL_DIR)/Make-default.mk

%: %.cc
	cd $(SMPL_DIR); ${MAKE}
	mkdir -p bin/
	$(CXX) $(CFLAGS) $(INCS) -DMAKE_MAIN -o bin/$@ $< ${MAIN_LIB} $(LIBS)

