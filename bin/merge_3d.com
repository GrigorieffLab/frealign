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

set mp_cpus		= `grep mp_cpus $SCRATCH/mparameters_run | awk '{print $2}'`
#
set m = ""
if ( ${mp_cpus} > 1 ) then
  set m = "_mp"
  setenv NCPUS ${mp_cpus}
endif

time ${bin_dir}/merge_3d${m}.exe < $SCRATCH/merge_3d_r${1}.in >& $SCRATCH/merge_3d_r${1}.log

echo 'merge_3d.com finished' >> $SCRATCH/merge_3d_r${1}.log
date
#
