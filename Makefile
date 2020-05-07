#
# CCFD Documentation Makefile
#
# 'make'          : create the documentation with doxygen
# 'make latex'    : create latex document with doxygen

.PHONY: doc latex clean

all: doc

ccfd:
	-@git clone https://github.com/hhh95/ccfd.git

doc: ccfd
	-@doxygen Doxyfile

latex: doc
	-@cd latex && $(MAKE) && mv refman.pdf ../refman.pdf

clean:
	-rm -rf html/*
	-rm -rf latex/*
	-rm -rf ccfd
