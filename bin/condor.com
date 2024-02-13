#!/bin/csh -f
#
universe = parallel
machine_count = 1
request_cpus = 1
request_memory = 6 GB
executable = /opt/condor/etc/examples/openmpiscript 
arguments = -np 320 /mnt/lustre1/users/luwei/em/relion-1.3/bin/relion_refine_mpi --o Class2D/run1 --i particles_autopick_sort_2d.star --particle_diameter 170 --angpix 1.06 --ctf  --iter 30 --tau2_fudge 2 --K 100 --flatten_solvent  --zero_mask  --oversampling 1 
--psi_step 6 --offset_range 5 --offset_step 2 --norm --scale  --j 1 --memory_per_thread 2
output = log/class2d_output.$(Node)
error = log/class2d_error.$(Node)
log = log/class2d.log
queue
