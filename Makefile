#
# CCFD Documentation Makefile
#
# 'make'      : create the documentation with doxygen
# 'make latex'    : create latex document with doxygen

.PHONY: doc latex docclean

all: doc

doc:
	-@cd ccfd && git pull
	-@doxygen Doxyfile

latex: doc
	-@cd latex && $(MAKE) && mv refman.pdf ../refman.pdf

docclean:
	-rm -rf html/*
	-rm -rf latex/*
