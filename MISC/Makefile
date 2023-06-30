#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
TOP := $(CURDIR)/..

include ../Makefile.compiler
include ../Makefile.include

LIBRARY = misc.a
FILES = locate.o sorti.o memusg.o binsrch.o binsrch_decr.o warn.o lcase.o specint.o

.PHONY: all clean

all: $(LIBRARY)
Makefile.depend:
	$(TOP)/depend.sh $(FILES) >Makefile.depend
clean:
	rm -f *.o *.mod *.MOD Makefile.depend $(LIBRARY)

$(LIBRARY): $(FILES)

include Makefile.depend
