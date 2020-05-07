#
# CCFD Documentation Makefile
#
# 'make'          : create the documentation with doxygen

.PHONY: doc clean cleandoc

all: doc

ccfd:
	-@git clone https://github.com/hhh95/ccfd.git

doc: ccfd
	-@doxygen Doxyfile
	-@cd latex && $(MAKE) && mv refman.pdf ../ccfd.pdf

clean:
	-rm -rf ccfd

cleandoc:
	-rm -rf html/*
	-rm -rf latex/*
