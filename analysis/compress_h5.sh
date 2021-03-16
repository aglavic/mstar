#!/bin/bash
# compress McStas hdf5 files to increase performance and save disk space

for fi in $*
do
  echo $fi
  h5repack -f /entry1/data/tof_detector_list_p_x_y_t_L/events:GZIP=5 \
           -l /entry1/data/tof_detector_list_p_x_y_t_L/events:CHUNK=3072x5 \
           $fi $fi.c
  mv $fi.c $fi
done
