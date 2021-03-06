#!/usr/bin/bash

./compile_if_needed.sh


# # 10x10 sample, two angles to get >0.3 Å^-1
# rm -r results/sample_10 results/reference_10 results/sample_34 results/reference_34
# mpirun --use-hwthread-cpus ./mstar.out -d results/sample_10 -n 5e9 --format=NeXuS \
       # omega=1.0 reference=0 skipping_chopp=0
# nohup analysis/compress_h5.sh results/sample_10/mccode.h5&
# mpirun --use-hwthread-cpus ./mstar.out -d results/reference_10 -n 5e9 --format=NeXuS \
       # omega=1.0 reference=1 skipping_chopp=0
# nohup analysis/compress_h5.sh results/reference_10/mccode.h5&

# mpirun --use-hwthread-cpus ./mstar.out -d results/sample_34 -n 5e8 --format=NeXuS \
       # omega=3.4 reference=0 skipping_chopp=0
# nohup analysis/compress_h5.sh results/sample_34/mccode.h5&
# mpirun --use-hwthread-cpus ./mstar.out -d results/reference_34 -n 5e8 --format=NeXuS \
       # omega=3.4 reference=1 skipping_chopp=0
# nohup analysis/compress_h5.sh results/reference_34/mccode.h5&

# # 10x10 sample, single angle with partial pulse skipping
# rm -r results/sample_20 results/reference_20 
# mpirun --use-hwthread-cpus ./mstar.out -d results/sample_20 -n 1e10 --format=NeXuS \
       # omega=2.0 reference=0 skipping_chopp=1
# analysis/compress_h5.sh results/sample_20/mccode.h5&
# mpirun --use-hwthread-cpus ./mstar.out -d results/reference_20 -n 1e10 --format=NeXuS \
       # omega=2.0 reference=1 skipping_chopp=1
# analysis/compress_h5.sh results/reference_20/mccode.h5&

# 2x2 sample, two angles to get >0.3 Å^-1
rm -r results/sample_2x2_10 results/sample_2x2_34
mpirun --use-hwthread-cpus ./mstar.out -d results/sample_2x2_10 -n 1e10 --format=NeXuS \
       omega=1.0 reference=0 skipping_chopp=0 \
	   sample_length=0.002 sample_height=0.002 over_illumination=0.0002
nohup analysis/compress_h5.sh results/sample_2x2_10/mccode.h5&
mpirun --use-hwthread-cpus ./mstar.out -d results/sample_2x2_34 -n 2e10 --format=NeXuS \
       omega=3.4 reference=0 skipping_chopp=0 \
	   sample_length=0.002 sample_height=0.002 over_illumination=0.0002
nohup analysis/compress_h5.sh results/sample_2x2_34/mccode.h5&

# 10x10 sample, semi-conventional with 2% resolution
rm -r results/reference_coll results/sample_coll_04 results/sample_coll_15 results/sample_coll_60
mpirun --use-hwthread-cpus ./mstar.out -d results/reference_coll -n 1e9 --format=NeXuS \
       omega=0.0 reference=0 skipping_chopp=0 lambda_min=3.9 \
	   divergence_slit_H=0.08
nohup analysis/compress_h5.sh results/reference_coll/mccode.h5&
mpirun --use-hwthread-cpus ./mstar.out -d results/sample_coll_04 -n 1e10 --format=NeXuS \
       omega=0.4 reference=0 skipping_chopp=0 lambda_min=3.9 \
	   divergence_slit_H=0.008
nohup analysis/compress_h5.sh results/sample_coll_04/mccode.h5&
mpirun --use-hwthread-cpus ./mstar.out -d results/sample_coll_15 -n 2e9 --format=NeXuS \
       omega=1.5 reference=0 skipping_chopp=0 lambda_min=3.9 \
	   divergence_slit_H=0.03
nohup analysis/compress_h5.sh results/sample_coll_15/mccode.h5&
mpirun --use-hwthread-cpus ./mstar.out -d results/sample_coll_60 -n 2e9 --format=NeXuS \
       omega=6.0 reference=0 skipping_chopp=0 lambda_min=3.9 \
	   divergence_slit_H=0.12
nohup analysis/compress_h5.sh results/sample_coll_60/mccode.h5&
