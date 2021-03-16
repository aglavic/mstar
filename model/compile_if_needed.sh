#!/bin/bash

if [ mstar.instr -nt mstar.out ] || [ ! -f mstar.out ]; then
	rm mstar.c mstar.out
	mcstas -o mstar.c mstar.instr $*
	mpicc -O2 -o mstar.out mstar.c -lm -DUSE_MPI -lmcpl \
	   -Wno-unused-result -Wno-format-truncation -Wno-format-overflow -Wno-stringop-overflow
	   # -DUSE_NEXUS -lNeXus \
fi
