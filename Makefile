#
# CCFD Documentation Makefile
#
# 'make'      : create the documentation with doxygen
# 'make latex'    : create latex document with doxygen

.PHONY: doc latex docclean

all: doc latex

ccfd:
	-@git clone https://github.com/hhh95/ccfd.git

doc: ccfd
	-@cd ccfd && git pull
	-@doxygen Doxyfile

latex: doc
	-@cd docs/latex && $(MAKE) && mv refman.pdf ../../refman.pdf

docclean:
	-rm -rf html/*
