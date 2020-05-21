# CCFD Makefile

include config.mk

### Build options:
TARGET = ccfd
GCC    = cc
ICC    = icc
BINDIR = bin
OBJDIR = obj
SRCDIR = src
LIBDIR = lib

### Library options:
LIBS	 = -lm
CGNS_DIR = $(LIBDIR)/CGNS-$(CGNS_VERSION)
CGNS_LIB = $(CGNS_DIR)/BUILD/src/libcgns.a
INCDIR   = -I $(CGNS_DIR)/BUILD/include
LIBS    += -L $(CGNS_DIR)/BUILD/lib -lcgns

### Compile- and linkflags:
ifeq ($(COMPILER), gnu)
  CC    = $(GCC)
  FLAGS = -std=c99 -Wall -Wextra -pedantic -Wno-unknown-pragmas
  ifeq ($(DEBUG), on)
    FLAGS += -g -O0 -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined
  else
    FLAGS += -Ofast -flto -march=native
  endif
  ifeq ($(PARALLEL), on)
    FLAGS += -fopenmp
  endif
  ifeq ($(PROF), on)
    FLAGS += -pg
  endif
  CFLAGS = $(FLAGS) $(INCDIR) -D$(EQNSYS)
  LFLAGS = $(FLAGS)
endif
ifeq ($(COMPILER), intel)
  CC    = $(ICC)
  FLAGS = -std=c99 -Wall -Wno-unknown-pragmas
  ifeq ($(DEBUG), on)
    FLAGS += -g -O0
  else
    FLAGS += -Ofast -ipo -xHost -mtune=native
  endif
  ifeq ($(PARALLEL), on)
    FLAGS += -qopenmp
  endif
  ifeq ($(PROF), on)
    FLAGS += -pg
  endif
  CFLAGS = $(FLAGS) $(INCDIR) -D$(EQNSYS)
  LFLAGS = $(FLAGS)
endif

### Build directions:
.PHONY: clean allclean check cleancheck

SRC = $(wildcard $(SRCDIR)/*.c)
OBJ = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRC))
TGT = $(BINDIR)/$(TARGET)

all: libs $(OBJDIR) $(BINDIR) $(TGT)
libs: $(CGNS_LIB)

$(OBJDIR):
	-mkdir -p $(OBJDIR)

$(BINDIR):
	-mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/%.h Makefile $(CGNS_LIB) config.mk
	$(CC) $(CFLAGS) -c $< -o $@

$(TGT): $(OBJ)
	$(CC) $(LFLAGS) $^ -o $@ $(LIBS)

$(CGNS_LIB): $(CGNS_DIR)
	-@mkdir $(CGNS_DIR)/BUILD && \
	cd $(CGNS_DIR)/BUILD && \
	cmake \
		-Wno-dev \
		-DCMAKE_C_COMPILER=$(CC) \
		-DCMAKE_INSTALL_PREFIX=. \
		-DCGNS_BUILD_CGNSTOOLS=OFF \
		-DCGNS_BUILD_SHARED=OFF \
		-DCGNS_BUILD_TESTING=OFF \
		-DCGNS_USE_SHARED=ON \
		-DCGNS_ENABLE_64BIT=ON \
		-DCGNS_ENABLE_BASE_SCOPE=OFF \
		-DCGNS_ENABLE_FORTRAN=OFF \
		-DCGNS_ENABLE_HDF5=OFF \
		-DCGNS_ENABLE_LEGACY=OFF \
		-DCGNS_ENABLE_MEM_DEBUG=OFF \
		-DCGNS_ENABLE_SCOPING=OFF \
		-DCGNS_ENABLE_TESTS=OFF \
		.. && \
	$(MAKE) install

$(CGNS_DIR):
	-tar -xzf $(CGNS_DIR).tar.gz -C $(LIBDIR)

clean:
	-rm -rf $(OBJDIR)
	-rm -rf $(BINDIR)

cleancheck:
	-rm -f check/*.csv check/*.log

allclean: clean cleancheck
	-rm -rf $(LIBDIR)/CGNS-*/

check:
	@cd check && python3 check.py ../$(TGT) $(EQNSYS)
