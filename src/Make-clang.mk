#	$Id: Makefile,v 1.3 2006/11/06 23:42:32 olethros Exp olethros $	
# My packages

PACKAGES = core algorithms geometry models statistics environments #florian
# Compiler and Linker
CC = clang++ # - some people may have to use clang instead
CXX = clang++ 
LD = clang # - you can link with 'ld' too
AR = ar -rus
DEP = clang -MM -D__DEPEND__


# Libraries

LIBS_EXPORT=$(SMPL_DIR)/export/lib
INCS_EXPORT=$(SMPL_DIR)/export/inc
MYLIBS = -L$(LIBS_EXPORT)
MYINCS = -I$(INCS_EXPORT)

# Flags

## Use DBG to compile a debug version.

#DBG_OPT=DBG
DBG_OPT=OPT

# Add -pg flag for profiling
CFLAGS_DBG = -fPIC -g -std=c++14 -Wall -DUSE_DOUBLE -Wno-overloaded-virtual
CFLAGS_OPT = -fPIC -g -std=c++14 -O3 -Wall -DUSE_DOUBLE -DNDEBUG -Wno-overloaded-virtual
#CFLAGS_DBG = -fPIC -g -Wall -pipe -pg
#CFLAGS_OPT = -fPIC -g -O3 -Wall -DNDEBUG -pipe -pg
CFLAGS=$(CFLAGS_$(DBG_OPT))
CXXFLAGS=$(CFLAGS)

# DIRECTORIES - You might want to change those

LIB_DIR_DBG=lib_dbg
LIB_DIR_OPT=lib
LIB_DIR_NAME=$(LIB_DIR_$(DBG_OPT))
OBJ_DIR_DBG=objs_dbg
OBJ_DIR_OPT=objs
OBJ_DIR_NAME=$(OBJ_DIR_$(DBG_OPT))



##SMPL_DIR := $(shell pwd)
BINDIR = $(SMPL_DIR)/bin
BIN_INSTALL_DIR = $(HOME)/bin
LIBS_DIR = $(SMPL_DIR)/$(LIB_DIR_NAME)
OBJS_DIR = $(SMPL_DIR)/$(OBJ_DIR_NAME)
LIBSMPL = $(LIBS_DIR)/libsmpl.a
LIBSMPLXX = $(LIBS_DIR)/libsmpl++.a
LIBS = -L$(LIBS_DIR) $(MYLIBS) -latlas -lcblas -lgsl
EXPORTED_LIBS = -lranlib
MAIN_LIB = -lsmpl
INCS := -I$(SMPL_DIR)/core $(MYINCS)
INCS += $(foreach f,$(PACKAGES),-I$(SMPL_DIR)/$(f))

# -lnjamd -lefence - add either for debugging

C_FILES := $(wildcard *.c)
OBJS := $(foreach f,$(C_FILES),$(OBJS_DIR)/$(patsubst %.c,%.o,$(f)))


