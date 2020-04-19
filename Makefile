#
# CCFD Makefile
#
# 'make' 	  : build libraries and executable
# 'make libs' 	  : build only libraries
# 'make clean' 	  : remove executable
# 'make allclean' : remove executable and libraries
# 'make doc'      : create the documentation with doxygen

### Equation system:
EQNSYS = EULER
#EQNSYS = NAVIERSTOKES

### Build options:
TARGET = ccfd
CC     = gcc
BINDIR = bin
OBJDIR = obj
SRCDIR = src
LIBDIR = lib

### Library options:
# Math
LIBS	     = -lm
# CGNS
CGNS_VERSION = 3.1.4
#CGNS_VERSION = 4.1.1
CGNS_DIR     = $(LIBDIR)/CGNS-$(CGNS_VERSION)
CGNS_LIB     = $(CGNS_DIR)/BUILD/src/libcgns.a
INCDIR       = -I $(CGNS_DIR)/BUILD/include
LIBS        += -L $(CGNS_DIR)/BUILD/lib -lcgns

### Compile- and linkflags:
FLAGS  = -std=c99 -Wall -pedantic -Wno-unknown-pragmas -march=native -g -O3 -fopenmp
CFLAGS = $(FLAGS) $(INCDIR) -D $(EQNSYS)
LFLAGS = $(FLAGS)

### Build directions:
.PHONY: clean allclean doc latex docclean

SRC = $(wildcard $(SRCDIR)/*.c)
OBJ = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRC))
TGT = $(BINDIR)/$(TARGET)

default: libs $(OBJDIR) $(BINDIR) $(TGT)
all: default
libs: $(CGNS_LIB)

$(OBJDIR):
	-mkdir -p $(OBJDIR)

$(BINDIR):
	-mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/%.h Makefile $(CGNS_LIB)
	$(CC) $(CFLAGS) -c $< -o $@

$(TGT): $(OBJ)
	$(CC) $(LFLAGS) $^ -o $@ $(LIBS)

$(CGNS_LIB) : $(CGNS_DIR)
	-@mkdir $(CGNS_DIR)/BUILD && \
	cd $(CGNS_DIR)/BUILD && \
	cmake \
		-Wno-dev \
		-DCMAKE_INSTALL_PREFIX=. \
		-DBUILD_CGNSTOOLS=OFF \
		-DCGNS_BUILD_SHARED=OFF \
		-DCGNS_USE_SHARED=ON \
		-DENABLE_64BIT=ON \
		-DENABLE_FORTRAN=OFF \
		-DENABLE_HDF5=OFF \
		-DENABLE_LEGACY=OFF \
		-DENABLE_SCOPING=OFF \
		-DENABLE_TESTS=OFF \
		-DCGNS_BUILD_CGNSTOOLS=OFF \
		-DCGNS_BUILD_SHARED=OFF \
		-DCGNS_BUILD_TESTING=OFF \
		-DCGNS_ENABLE_64BIT=ON \
		-DCGNS_ENABLE_BASE_SCOPE=OFF \
		-DCGNS_ENABLE_FORTRAN=OFF \
		-DCGNS_ENABLE_HDF5=OFF \
		-DCGNS_ENABLE_LEGACY=OFF \
		-DCGNS_ENABLE_MEM_DEBUG=OFF \
		-DCGNS_ENABLE_SCOPING=OFF \
		-DCGNS_ENABLE_TESTS=OFF \
		-DCGNS_USE_SHARED=ON \
		.. && \
	$(MAKE) install

$(CGNS_DIR):
	-tar -xzf $(CGNS_DIR).tar.gz -C $(LIBDIR)

doc:
	-@doxygen Doxyfile

latex: doc
	-@cd docs/latex && $(MAKE) && mv refman.pdf ../../ccfd.pdf

docclean:
	-rm -rf docs/*

clean:
	-rm -rf $(OBJDIR)
	-rm -rf $(BINDIR)

allclean: clean
	-rm -rf $(LIBDIR)/CGNS-*/
