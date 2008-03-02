SMPL_DIR := $(shell pwd)
include Make-default.mk
SUBDIRS := libs
SUBDIRS += $(PACKAGES)

all: 
	@echo $SMPL_DIR
	@mkdir -p ${LIBS_DIR}
	@mkdir -p ${OBJS_DIR}
	@mkdir -p ${LIBS_EXPORT}
	@mkdir -p ${INCS_EXPORT}
	@for subdir in ${SUBDIRS}; do ( cd $$subdir ; ${MAKE} $@) || exit 10; done
	@echo "Done."

clean:
	@echo "Making clean"
	@\rm -f .deps*
	@\rm -f *~
	@\rm -rf ${LIBS_EXPORT}
	@\rm -rf ${INCS_EXPORT}
	@for subdir in ${SUBDIRS} ; do ( cd $$subdir ; ${MAKE} $@ ) || exit 10 ; done

distclean:
	@\rm -f .deps*

.depend: depend

depend:
	@echo $(SMPL_DIR)
	@echo ${LIBS_EXPORT}
	@echo ${INCS_EXPORT}
	@echo $(OBJS_DIR)
	@echo $(LIBS_DIR)
	@mkdir -p ${LIBS_EXPORT}
	@mkdir -p ${INCS_EXPORT}
	@echo ">>> Dependencies <<<"
	@echo $(INCS)
	@echo "--------------------"
	@for subdir in ${SUBDIRS} ; do ( cd $$subdir ; ${MAKE} $@ ) || exit 10 ; done
	@mkdir -p $(OBJS_DIR)
	@mkdir -p $(LIBS_DIR)
	@echo ">> !!! Ok !!! <<<"




