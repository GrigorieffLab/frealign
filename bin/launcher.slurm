#!/bin/csh
#
# Simple SGE script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#      <------------ Account String ------------>
# <--- (Use this ONLY if you have MULTIPLE accounts) --->
##SBATCH -A 
#------------------------------------------------------

set working_directory   = `pwd`
set SCRATCH             = `grep scratch_dir mparameters | awk '{print $2}'`
if ( $status || $SCRATCH == "" ) then
  set SCRATCH           = ${working_directory}/scratch
  setenv CONTROL_FILE   scratch/slurm.com
else
 set back = `echo $working_directory | fgrep -o / | awk '{printf "%3s","../"}'`
 setenv CONTROL_FILE   ${back}${SCRATCH}/slurm.com
endif

module load launcher
setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
# setenv CONTROL_FILE   ${SCRATCH}/slurm.com
# setenv CONTROL_FILE   ../../../..${SCRATCH}/slurm.com
# setenv CONTROL_FILE   scratch/slurm.com
setenv WORKDIR        $working_directory
# 
# Variable description:
#
#  EXECUTABLE     = full path to the job launcher executable
#  CONTROL_FILE   = text input file which specifies
#                   executable for each process
#                   (should be located in WORKDIR)
#  WORKDIR        = location of working directory
#
#      <------ End Setup Parameters ------>
#--------------------------------------------------------
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
