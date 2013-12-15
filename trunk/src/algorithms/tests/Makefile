SMPL_DIR := $(shell cd ../..; pwd)

include $(SMPL_DIR)/Make-default.mk

sample_agent_rl_glue: sample_agent_rl_glue.cc
	cd $(SMPL_DIR); ${MAKE}
	mkdir -p bin/
	$(CXX) $(CFLAGS) $(INCS) -DMAKE_MAIN -o bin/$@ $< ${MAIN_LIB} $(LIBS) -lrlutils -lrlgluenetdev -lrlagent 


%: %.cc
	cd $(SMPL_DIR); ${MAKE}
	mkdir -p bin/
	$(CXX) $(CFLAGS) $(INCS) -DMAKE_MAIN -o bin/$@ $< ${MAIN_LIB} $(LIBS)

