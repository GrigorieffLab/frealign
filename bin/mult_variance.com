#!/bin/csh -f
#
#   Control script to submit multiple jobs on a cluster using
#   the Sun Grid Engine, or on a multi-core machine.
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

set first		= ${1}
set last		= ${2}
set start		= ${3}
set data_input		= `grep data_input $SCRATCH/mparameters_run | awk '{print $2}'`
set pix			= `grep pix_size $SCRATCH/mparameters_run | awk '{print $2}'`
set nclass		= `grep nbootstrap $SCRATCH/mparameters_run | awk '{print $2}'`
set npar = `ls ${data_input}_${start}_r*.par | wc -l`
set bin_dir		= `grep frealign_bin_dir $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $bin_dir == "" ) then
  set bin_dir		= `which frealign_v9.exe`
  set bin_dir		= ${bin_dir:h}
endif
set cluster_type	= `cat $SCRATCH/cluster_type.log`
if ( $status || $cluster_type == "" ) then
  set cluster_type	= `grep cluster_type $SCRATCH/mparameters_run | awk '{print $2}'`
endif
set stn			= `grep qsub_string_rec mparameters | awk -F\" '{print $2}'`
set no_delete           = `grep delete_scratch $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $no_delete == "F" ) then
  set no_delete = 1
else
  set no_delete = 0
endif

set raw_images		= `grep raw_images_rec $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $raw_images == "" ) set raw_images = `grep raw_images $SCRATCH/mparameters_run | awk '{print $2}'`
set raw_images = `echo ${raw_images:r}`
set extension = `ls $raw_images.* | head -1`
if ( $extension == "" ) then
  set raw_images = ${working_directory}/${raw_images}
  set extension = `ls $raw_images.* | head -1`
endif
set extension = `echo ${extension:e}`

set nx = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text NX | awk '{print $4}'`
set mem_small = `echo $nx | awk '{print int(10 * $1^3 * 4 * 3 /1024^3 + 1)/10}'`  
set mem_small = `echo $mem_small | awk '{if ($1 < 1) {print 1} else {print $1} }'`

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then
 set mem_small = `echo $mem_small | awk '{print int($1 * 1024)}'`
endif

cp mparameters $SCRATCH/mparameters_run
set nproc		= `grep nprocessor_rec $SCRATCH/mparameters_run | awk '{print $2}'`
set sym			= `grep Symmetry $SCRATCH/mparameters_run | awk '{print $2}'`

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
  set hosts = hosts
  \rm $SCRATCH/hosts >& /dev/null
  foreach h (`cat $hosts | awk '{print $1}'`)
    set n = `grep $h $hosts | awk '{print $2}'`
    awk 'BEGIN { while (i++ < '$n') print "'$h'" }' >> $SCRATCH/hosts
  end
  set nhosts = `wc -l $SCRATCH/hosts | awk '{print $1}'`
  set nh = 1
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set nproc = `echo ${nproc} | awk '{print int($1/16+0.5)*16}'`
endif

# echo "Calculating bootstrap volumes...." >> frealign.log

${bin_dir}/rsample.exe << eot
${data_input}_${start}_r${4}.par
${pix}
${nclass}
$SCRATCH/${data_input}_${start}_r.par
eot

\rm $SCRATCH/pid_temp.log >& /dev/null
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  if ( $sym == "H" || $sym == "HP" ) then
    echo "executable = ${bin_dir}/mult_hreconstruct_v.com" > $SCRATCH/condor.com
  else
    echo "executable = ${bin_dir}/mult_reconstruct_v.com"  > $SCRATCH/condor.com
  endif
  echo "universe = vanilla"                               >> $SCRATCH/condor.com
  echo "request_memory = ${mem_small} GB"                 >> $SCRATCH/condor.com
endif

  set nc = 1
  set firstn = $first
  set lastn = $last

  while ( $nc <= $nclass )

    echo "Cycle "$start": reconstructing bootstrap volume $nc on "`date` >> frealign.log

    grep kill $SCRATCH/pid.log >& /dev/null
    if ( ! $status ) exit

    \rm $SCRATCH/${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >& /dev/null
    \rm $SCRATCH/${data_input}_${start}_r${nc}.shft_* >& /dev/null
    \rm $SCRATCH/${data_input}_${start}_r${nc}.${extension} >& /dev/null
    \rm $SCRATCH/${data_input}_mult_refine_n_r${nc}.log_* >& /dev/null
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
      if ( $sym == "H" || $sym == "HP" ) then
        set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N r${nc}_${firstn} ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc | awk '{print $3}'`
      else
        set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N r${nc}_${firstn} ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc | awk '{print $3}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
      if ( $sym == "H" || $sym == "HP" ) then
        set pid = `sbatch $stn -D $working_directory -J r${nc}_${firstn} -o /dev/null ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc | awk '{print $4}'`
      else
        set pid = `sbatch $stn -D $working_directory -J r${nc}_${firstn} -o /dev/null ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc | awk '{print $4}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
      if ( $sym == "H" || $sym == "HP" ) then
#        set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
        set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
      else
#        set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
        set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
      if ( $sym == "H" || $sym == "HP" ) then
        echo ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc \& >> $SCRATCH/slurm.com
      else
        echo ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc \& >> $SCRATCH/slurm.com
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
      echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${firstn}_${nc}.com
      echo cd $working_directory >> $SCRATCH/pbs_${firstn}_${nc}.com
      if ( $sym == "H" || $sym == "HP" ) then
        echo ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc >> $SCRATCH/pbs_${firstn}_${nc}.com
      else
        echo ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc >> $SCRATCH/pbs_${firstn}_${nc}.com
      endif
      set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${firstn}_${nc}.com`
      set pid = `echo $pid | awk -F. '{print $1}'`
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
      echo "arguments = $firstn $lastn $start $nc" >> $SCRATCH/condor.com
      echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
      echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
      echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
      echo "queue" >> $SCRATCH/condor.com
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
      if ( $nh > $nhosts ) set nh = 1
      set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
      @ nh++
      if ( $sym == "H" || $sym == "HP" ) then
        ssh $rh "cd $working_directory; ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc" &
      else
        ssh $rh "cd $working_directory; ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc" &
      endif
      set pid = $!
    else
      if ( $sym == "H" || $sym == "HP" ) then
        ${bin_dir}/mult_hreconstruct_v.com $firstn $lastn $start $nc &
      else
        ${bin_dir}/mult_reconstruct_v.com $firstn $lastn $start $nc &
      endif
      set pid = $!
    endif

    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" ) then
      echo $pid ${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >> $SCRATCH/pid_temp.log
    endif
    @ nc++
  end

#  sleep 2

submission_done:

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then          
  chmod +x $SCRATCH/slurm.com
  set incs = `grep com $SCRATCH/slurm.com | wc -l`
  set nnode = `echo $incs | awk '{print int(($1+15)/16)}'`
  set pid = `sbatch  -D $working_directory -J recon -N $nnode -n $incs -p normal -o $SCRATCH/frealign.o%j -t 12:00:00 ${bin_dir}/launcher.slurm`
  set pid = `echo $pid | awk '{ print $NF }'`
  echo $pid ${data_input}_mult_reconstruct_${first}_${last}.log >> $SCRATCH/pid_temp.log
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  set pid = `condor_submit $SCRATCH/condor.com | tail -1`
  set pid = `echo $pid | awk '{print $6}' | awk -F. '{print $1}'`
  set i = 0
  foreach fst (`grep arguments $SCRATCH/condor.com | awk '{print $3}'`)
    set nc = `grep arguments scratch/condor.com | cat -n | awk '{if ($1 - 1 == '$i') print $7}'`
    echo ${pid}.${i} ${data_input}_mult_reconstruct_r${nc}_n${fst}.log >> $SCRATCH/pid_temp.log
    @ i++
  end
endif

sleep 10

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
\rm $SCRATCH/pid_temp.log >& /dev/null

set firstn = $first
set lastn = $last

set nc = 1
while ( $nc <= $nclass )

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

checkdone_rec:

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

#  sleep 5

    grep --binary-files=text finished $SCRATCH/${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >& /dev/null

    if ($status) then
      sleep 5
      goto checkdone_rec
    endif

    if ($status) goto checkdone_rec

    echo "Cycle "$start": reconstruction of bootstrap volume $nc finished "`date` >> frealign.log
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" ) then
      set pid = `grep "${data_input}_mult_reconstruct_r${nc}_n${firstn}.log" $SCRATCH/pid.log | awk '{print $1}'`
      cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
      if ($pid != "") grep -v $pid $SCRATCH/pid_temp.log > $SCRATCH/pid.log
      \rm $SCRATCH/pid_temp.log >& /dev/null
    endif

  @ nc++

end

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set pid = `grep "${data_input}_mult_reconstruct_${first}_${last}" $SCRATCH/pid.log | awk '{print $1}'`
  cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
  if ($pid != "") grep -v $pid $SCRATCH/pid_temp.log > $SCRATCH/pid.log
  \rm $SCRATCH/pid_temp.log >& /dev/null
endif

\rm $SCRATCH/pid_temp.log >& /dev/null

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  echo "Cycle "$start": calculating 3D variance map "`date` >> frealign.log

  @ incr = 1 + $last - $first
  echo $SCRATCH/${data_input}_${start}_r.${extension} > $SCRATCH/bootvar.in
  echo 1,$nclass >> $SCRATCH/bootvar.in
  echo $incr >> $SCRATCH/bootvar.in
  echo 0 >> $SCRATCH/bootvar.in
  echo $SCRATCH/${data_input}_${start}_ave.${extension} >> $SCRATCH/bootvar.in
  echo $SCRATCH/${data_input}_${start}_var.${extension} >> $SCRATCH/bootvar.in

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
    set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N bootvar ${bin_dir}/bootvar.com | awk '{print $3}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
    set pid = `sbatch $stn -D $working_directory -J bootvar -o /dev/null ${bin_dir}/bootvar.com | awk '{print $4}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
#    set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/bootvar.com | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
    set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/bootvar.com | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then 
    echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
    echo ${bin_dir}/bootvar.com \& >> $SCRATCH/slurm.com
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
    echo \#\!\/bin\/csh -f > $SCRATCH/pbs.com
    echo cd $working_directory > $SCRATCH/pbs.com
    echo ${bin_dir}/bootvar.com >> $SCRATCH/pbs.com
    set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs.com`
    set pid = `echo $pid | awk -F. '{print $1}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
    echo "executable = ${bin_dir}/boorvar.com"   > $SCRATCH/condor.com
    echo "universe = vanilla"                   >> $SCRATCH/condor.com
    echo "request_memory = ${mem_small} GB"     >> $SCRATCH/condor.com
    echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
    echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
    echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
    echo "queue" >> $SCRATCH/condor.com
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
    if ( $nh > $nhosts ) set nh = 1
    set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
    @ nh++
    ssh $rh "cd $working_directory; ${bin_dir}/bootvar.com" &
    set pid = $!
  else
    ${bin_dir}/bootvar.com &
    set pid = $!
  endif

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" ) then
    echo $pid bootvar.log >> $SCRATCH/pid_temp.log
  endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  chmod +x $SCRATCH/slurm.com
  set incs = `grep com $SCRATCH/slurm.com | wc -l`
  set nnode = `echo $incs | awk '{print int(($1+15)/16)}'`
  set pid = `sbatch  -D $working_directory -J merge -N $nnode -n $incs -p normal -o $SCRATCH/frealign.o%j -t 12:00:00 ${bin_dir}/launcher.slurm`
  set pid = `echo $pid | awk '{ print $NF }'`
  echo $pid bootvar.log >> $SCRATCH/pid_temp.log
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  set pid = `condor_submit $SCRATCH/condor.com | tail -1`
  set pid = `echo $pid | awk '{print $6}' | awk -F. '{print $1}'`
  echo ${pid}.1 bootvar.log >> $SCRATCH/pid_temp.log
endif

sleep 10

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
\rm $SCRATCH/pid_temp.log >& /dev/null

checkdonem:

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  sleep 2

  grep --binary-files=text "bootvar.com finished" $SCRATCH/bootvar.log >& /dev/null
  if ($status) goto checkdonem

  ls $SCRATCH/${data_input}_${start}_ave.${extension} >& /dev/null
  if ($status) goto checkdonem

  ls $SCRATCH/${data_input}_${start}_var.${extension} >& /dev/null
  if ($status) goto checkdonem

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE") then
    set pid = `grep "bootvar.log" $SCRATCH/pid.log | awk '{print $1}'`
    cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
    if ($pid != "") grep -v $pid $SCRATCH/pid_temp.log > $SCRATCH/pid.log
    \rm $SCRATCH/pid_temp.log >& /dev/null
  endif

  sleep 1

  echo "Cycle "$start": variance calculation for class ${4} finished "`date` >> frealign.log
  mv $SCRATCH/${data_input}_${start}_ave.${extension} ${working_directory}/${data_input}_${start}_r${4}_ave.${extension}
  mv $SCRATCH/${data_input}_${start}_var.${extension} ${working_directory}/${data_input}_${start}_r${4}_var.${extension}

  set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
  if ( $form == "IMAGIC" ) then
    mv $SCRATCH/${data_input}_${start}_ave.img ${working_directory}/${data_input}_${start}_r${4}_ave.img
    mv $SCRATCH/${data_input}_${start}_var.img ${working_directory}/${data_input}_${start}_r${4}_var.img
  endif

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set pid = `grep "bootvar" $SCRATCH/pid.log | awk '{print $1}'`
  cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
  if ($pid != "") grep -v $pid $SCRATCH/pid_temp.log > $SCRATCH/pid.log
  \rm $SCRATCH/pid_temp.log >& /dev/null
endif

# Clean up

@ prev = $start - 1

if ( ! $no_delete ) then
  set nc = 1
  while ( $nc <= $nclass )

    \rm $SCRATCH/${data_input}_${start}_r*.* >& /dev/null
    @ nc++

  end
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  \rm $SCRATCH/slurm.com >& /dev/null
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  \rm $SCRATCH/condor.* >& /dev/null
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
  \rm $SCRATCH/pbs*.com* >& /dev/null
endif

echo "Cycle "$start": variance calculation done." >> frealign.log

grep reconstruct $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) then
    echo "monitor term" >> $SCRATCH/pid.log
  endif
endif

date
