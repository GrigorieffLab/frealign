#!/bin/csh -f
limit coredumpsize 0

set working_directory	= `pwd -L`
set SCRATCH		= `grep scratch_dir mparameters | awk '{print $2}'`
if ( $status || $SCRATCH == "" ) then
  set SCRATCH		= ${working_directory}/scratch
endif
if ( ! -d $SCRATCH ) then
  mkdir $SCRATCH
endif

set bin_dir		= `grep frealign_bin_dir $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $bin_dir == "" ) then
  set bin_dir		= `which frealign_v9.exe`
  set bin_dir		= ${bin_dir:h}
endif

time ${bin_dir}/apply_mask.exe < $SCRATCH/apply_mask_r${1}.in >& $SCRATCH/apply_mask_r${1}.log

echo 'apply_mask.com finished' >> $SCRATCH/apply_mask_r${1}.log
date
#
