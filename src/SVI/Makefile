# get user and architecture specific options
SMPL_DIR := $(shell cd ..; pwd)
include ../Make-default.mk

CC_FILES := $(wildcard *.cc)
OBJS := $(foreach f,$(CC_FILES),$(OBJS_DIR)/$(patsubst %.cc,%.o,$(f)))

all:	$(LIBSMPL)

$(LIBSMPL): $(OBJS)
	@echo "Archiving..." $(LIBS_DIR)
	@$(AR) $(LIBSMPL) $(OBJS)

$(OBJS_DIR)/%.o: %.cc .depend
	@mkdir -p $(OBJS_DIR)
	@echo $<
	@$(CXX) $(CXXFLAGS) $(INCS) -o $@ -c $<

distclean:
	@\rm -f .depend

clean:
	@echo "Remove objects file and dependencies..."
	@\rm -Rf $(OBJS) $(LIBSMPL)
	@\rm -f .depend
	@\rm -f *~

.depend:
	@echo "Tracking dependencies..."
	@\rm -f .depend
	@for file in *.cc ; do printf "$(OBJS_DIR)/" >> .depend; $(DEP) $(CFLAGS) $(INCS) $$file >> .depend; done


ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),.depend)
include .depend
endif
endif
endif

# DO NOT DELETE
