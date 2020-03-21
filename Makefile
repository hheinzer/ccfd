#
# CCFD Makefile (copy)
#
# 'make' 	  : build libraries and executable
# 'make libs' 	  : build only libraries
# 'make clean' 	  : remove executable
# 'make allclean' : remove executable and libraries

.PHONY: default all clean allclean

default:
	cd src && $(MAKE) default

all: default

libs:
	cd src && $(MAKE) libs

clean:
	cd src && $(MAKE) clean

allclean:
	cd src && $(MAKE) allclean
