#!/bin/bash

if [ 4a_simple.instr -nt 4a_simple.out ] || [ ! -f 4a_simple.out ]; then
	rm 4a_simple.c 4a_simple.out
	mcstas -o 4a_simple.c 4a_simple.instr $*
	mpicc -O2 -o 4a_simple.out 4a_simple.c -lm -DUSE_MPI \
	   -Wno-unused-result -Wno-format-truncation -Wno-format-overflow -Wno-stringop-overflow
	   # -DUSE_NEXUS -lNeXus \
fi
