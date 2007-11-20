TOPDIRS=src

all:
	@for subdir in ${TOPDIRS}; do ( cd $$subdir ; ${MAKE} $@) || exit 10; done
	@echo "Done."

clean:
	@for subdir in ${TOPDIRS}; do ( cd $$subdir ; ${MAKE} $@) || exit 10; done
	@echo "Done."



