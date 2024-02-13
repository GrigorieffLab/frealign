#!/bin/csh -f
#
#   Control script to submit multiple jobs on a cluster using
#   the Sun Grid Engine. Each job processes a sampled set of
#   particles.
#
set working_directory	= `pwd -L`
set SCRATCH		= `grep scratch_dir mparameters | awk '{print $2}'`
if ( $status || $SCRATCH == "" ) then
  set SCRATCH		= ${working_directory}/scratch
endif
if ( ! -d $SCRATCH ) then
  mkdir $SCRATCH
endif

cp mparameters $SCRATCH/mparameters_run

if ( ${1} != "" ) then
  set start = ${1}
else
  set start		= `grep start_process $SCRATCH/mparameters_run | awk '{print $2}'`
endif
set first		= `grep first_particle $SCRATCH/mparameters_run | awk '{print $2}'`
set last		= `grep last_particle $SCRATCH/mparameters_run | awk '{print $2}'`
set data_input		= `grep data_input $SCRATCH/mparameters_run | awk '{print $2}'`
set pix			= `grep pix_size $SCRATCH/mparameters_run | awk '{print $2}'`
set nclass		= `grep nclasses $SCRATCH/mparameters_run | awk '{print $2}'`
set bin_dir		= `grep frealign_bin_dir $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $bin_dir == "" ) then
  set bin_dir		= `which frealign_v9.exe`
  set bin_dir		= ${bin_dir:h}
endif
set cluster_type	= `grep cluster_type $SCRATCH/mparameters_run | awk '{print $2}'`

set raw_images		= `grep raw_images_ref $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $raw_images == "" ) set raw_images = `grep raw_images $SCRATCH/mparameters_run | awk '{print $2}'`
set raw_images = `echo ${raw_images:r}`
set extension = `ls $raw_images.* | head -1`
if ( $extension == "" ) then
  set raw_images = ${working_directory}/${raw_images}
  set extension = `ls $raw_images.* | head -1`
endif
set extension = `echo ${extension:e}`

if ( $first == "" || $last == "" ) then
  set first = 1
  set last = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text NX | awk '{print $6}'`
endif

@ prev = $start - 1

mv ${data_input}_${prev}_r1.par ${data_input}_${prev}.par
mv ${data_input}_${prev}_r1.${extension} ${data_input}_${prev}.${extension}

time ${bin_dir}/rsample.exe << eot
${data_input}_${prev}.par
${pix}
${nclass}
${data_input}_${prev}_r.par
eot

# ${working_directory}/mult_reconstruct.com $first $last $prev

date
