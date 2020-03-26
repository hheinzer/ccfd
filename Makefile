#
# CCFD Makefile
#
# 'make' 	  : build libraries and executable
# 'make libs' 	  : build only libraries
# 'make clean' 	  : remove executable
# 'make allclean' : remove executable and libraries

### Equation system:
EQNSYS = euler
#EQNSYS = navierStokes

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
ARCH         = LINUX64
CGNS_VERSION = 4.1.0
CGNS_DIR     = CGNS-$(CGNS_VERSION)
CGNS_LIB_DIR = $(ARCH)
INCDIR       = -I $(CGNS_DIR)/$(CGNS_LIB_DIR)/include
LIBS        += -L $(CGNS_DIR)/$(CGNS_LIB_DIR)/lib -lcgns

### Compile- and linkflags:
FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -O0 -g
#FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -march=native -O3
#FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -march=native -O3 -fopenmp
#FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -march=native -O0 -g
#FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -march=native -O3 -pg
#FLAGS = -std=c99 -pedantic -Wall -Wno-unknown-pragmas -march=native -O3 -fopenmp -pg
INCDIR += -I $(SRCDIR) -I $(SRCDIR)/$(EQNSYS)
CFLAGS = $(FLAGS) $(INCDIR) -D $(EQNSYS)
LFLAGS = $(FLAGS)

### Build directions:
.PHONY: default all clean allclean

SRC_CORE   = $(wildcard $(SRCDIR)/*.c)
OBJ_CORE   = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRC_CORE))
SRC_EQNSYS = $(wildcard $(SRCDIR)/$(EQNSYS)/*.c)
OBJ_EQNSYS = $(patsubst $(SRCDIR)/$(EQNSYS)/%.c, $(OBJDIR)/%.o, $(SRC_EQNSYS))
SRC        = $(SRC_CORE) $(SRC_EQNSYS)
OBJ        = $(OBJ_CORE) $(OBJ_EQNSYS)
TGT        = $(BINDIR)/$(TARGET)

default: cgns $(OBJDIR) $(BINDIR) $(TGT)
all: default
libs: cgns

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/%.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/$(EQNSYS)/%.c $(SRCDIR)/$(EQNSYS)/%.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

$(TGT): $(OBJ)
	$(CC) -Wall $^ -o $@ $(LIBS)

cgns:
	@cd $(LIBDIR) && \
	if [ ! -f $(CGNS_DIR)/build/src/libcgns.a ]; then \
		printf "Building CGNS library\n"; \
		if [ ! -d $(CGNS_DIR) ]; then \
			tar -xzf $(CGNS_DIR).tar.gz >/dev/null; \
		fi; \
		mkdir $(CGNS_DIR)/build && \
		cd $(CGNS_DIR)/build && \
		cmake \
			-Wno-dev \
			-DCMAKE_INSTALL_PREFIX=$(shell pwd)/$(LIBDIR)/$(CGNS_DIR)/$(CGNS_LIB_DIR) \
			-DCGNS_BUILD_CGNSTOOLS=OFF \
			-DCGNS_BUILD_SHARED=OFF \
			-DCGNS_BUILD_TESTING=OFF \
			-DCGNS_ENABLE_64BIT=ON \
			-DCGNS_ENABLE_BASE_SCOPE=OFF \
			-DCGNS_ENABLE_FORTRAN=OFF \
			-DCGNS_ENABLE_HDF5=OFF \
			-DCGNS_ENABLE_MEM_DEBUG=OFF \
			-DCGNS_ENABLE_SCOPING=OFF \
			-DCGNS_ENABLE_TESTS=OFF \
			-DCGNS_ENABLE_LEGACY=OFF \
			.. && \
		$(MAKE) && \
		$(MAKE) install; \
	fi;

clean:
	-rm -rf $(OBJDIR)
	-rm -rf $(BINDIR)

allclean: clean
	-rm -rf $(LIBDIR)/$(CGNS_DIR)

run:
	@./$(TGT) calc/sod.ini

debug:
	@gdb -q -tui --args $(TGT) calc/sod.ini
