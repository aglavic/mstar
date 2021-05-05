#!/usr/bin/sh
./compile_if_needed.sh


mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_reference -n 2e8 omega=0.0 d_omega=0.005

mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_0p25deg -n 2e8 omega=0.25 d_omega=0.0025
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_0p45deg -n 2e8 omega=0.45 d_omega=0.0045
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_0p80deg -n 1e8 omega=0.80 d_omega=0.0080
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_1p46deg -n 1e8 omega=1.46 d_omega=0.0146
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_2p62deg -n 1e8 omega=2.62 d_omega=0.0262
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_4p72deg -n 1e8 omega=4.72 d_omega=0.0472
mpirun --use-hwthread-cpus ./4a_simple.out -d ../results/4a_ni_8p50deg -n 1e8 omega=8.50 d_omega=0.0850
