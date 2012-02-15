SMPL_DIR := $(shell cd ../..; pwd)
include ../../Make-default.mk
SOLIBRARY = libranlib.so
SOLIBDIR = .
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)

all: ${SOLIBRARY}

${SOLIBRARY}: ${OBJECTS}
	@cp ranlib.h $(INCS_EXPORT)
	${CC} -shared -o ${SOLIBRARY} ${OBJECTS} ${LDFALGS} ${LIBSPATH} ${DEBUG_LIBS}
	cp ${SOLIBRARY} $(LIBS_EXPORT)
	@$(AR) $(LIBSMPL) $(OBJECTS)

ifeq (.depend,$(wildcard .depend))
include .depend
endif

depend : .depend
	@cp ranlib.h ${INCS_EXPORT}

.depend:
	$(DEP) $(SOURCES) >.depend

clean:
	@\rm -f .depend
	@\rm -Rf $(OBJECTS)
	@\rm -f *~
	@\rm -f $(SOLIBRARY)
# DO NOT DELETE
