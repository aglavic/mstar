#!/usr/bin/bash

./compile_if_needed.sh


# beamline transmission with normal chopper
rm -r results/transmission_normal
mpirun --use-hwthread-cpus ./mstar.out -d results/transmission_normal -n 1e9 \
       omega=0.0 reference=0 skipping_chopp=0 lambda_min=2.5 sample_height=0.005 sample_length=0.001

# beamline transmission with hybrid pulse-skipping chopper
rm -r results/transmission_skip
mpirun --use-hwthread-cpus ./mstar.out -d results/transmission_skip -n 1e9 \
       omega=0.0 reference=0 skipping_chopp=1 lambda_min=2.5 sample_height=0.005 sample_length=0.001
