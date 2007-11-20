SMPL_DIR := $(shell cd ..; pwd)
SUBDIRS = ranlib
EXPINCDIRS = $(SUBDIRS)
include ../Make-default.mk

all: 
	@for subdir in ${SUBDIRS}; do ( cd $$subdir ; ${MAKE} $@) || exit 10; done
	@echo "Done."

clean:
	@\rm -f *~
	@for subdir in ${SUBDIRS} ; do ( cd $$subdir ; ${MAKE} $@ ) || exit 10 ; done	

distclean:
	@\rm -f .deps*

depend:
	@echo ">>> Dependencies <<<"
	@echo $(INCS)
	@echo "--------------------"
	@for subdir in ${SUBDIRS} ; do ( cd $$subdir ; ${MAKE} $@ ) || exit 10 ; done
	@\mkdir -p $(OBJS_DIR)
	@\mkdir -p $(LIBS_DIR)
	@echo ">> !!! Ok !!! <<<"
