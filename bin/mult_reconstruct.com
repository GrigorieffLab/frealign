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
set nclass		= `grep nclasses $SCRATCH/mparameters_run | awk '{print $2}'`
set npar = `ls ${data_input}_${start}_r*.par | wc -l`
if ( $nclass > $npar ) set nclass = $npar
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

mainloop:

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
@ incr = $last + 1 - $first
set incr = `echo ${incr} ${nclass} ${nproc} | awk '{print int($1*$2/$3+1)}'`
# set incr = `expr ${incr} \* ${nclass} / ${nproc}`
# @ incr++

reconstruct:
# echo "Calculating 3D structure...." >> frealign.log

\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  if ( $sym == "H" || $sym == "HP" ) then
    echo "executable = ${bin_dir}/mult_hreconstruct_n.com" > $SCRATCH/condor.com
  else
    echo "executable = ${bin_dir}/mult_reconstruct_n.com"  > $SCRATCH/condor.com
  endif
  echo "universe = vanilla"                               >> $SCRATCH/condor.com
  echo "request_memory = ${mem_small} GB"                 >> $SCRATCH/condor.com
endif

set npart = 0
set firstn = $first
@ lastn = $first + $incr - 1
if ( $lastn >= $last ) set lastn = $last

while ( $lastn <= $last )
  echo "Cycle "$start": reconstructing particles $firstn to $lastn on "`date` >> frealign.log

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  @ npart++

  set nc = 1
  while ( $nc <= $nclass )
    \rm $SCRATCH/${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >& /dev/null
#    \rm $SCRATCH/${data_input}_${start}_r${nc}.shft_* >& /dev/null
    \rm $SCRATCH/${data_input}_mult_refine_n_r${nc}.log_* >& /dev/null
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
      if ( $sym == "H" || $sym == "HP" ) then
        set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N r${nc}_${firstn} ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc | awk '{print $3}'`
      else
        set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N r${nc}_${firstn} ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc | awk '{print $3}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
      if ( $sym == "H" || $sym == "HP" ) then
        set pid = `sbatch $stn -D $working_directory -J r${nc}_${firstn} -o /dev/null ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc | awk '{print $4}'`
      else
        set pid = `sbatch $stn -D $working_directory -J r${nc}_${firstn} -o /dev/null ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc | awk '{print $4}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
      if ( $sym == "H" || $sym == "HP" ) then
#        set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
        set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
      else
#        set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
        set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
      if ( $sym == "H" || $sym == "HP" ) then
        echo ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc \& >> $SCRATCH/slurm.com
      else
        echo ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc \& >> $SCRATCH/slurm.com
      endif
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
      echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${firstn}_${nc}.com
      echo cd $working_directory >> $SCRATCH/pbs_${firstn}_${nc}.com
      if ( $sym == "H" || $sym == "HP" ) then
        echo ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc >> $SCRATCH/pbs_${firstn}_${nc}.com
      else
        echo ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc >> $SCRATCH/pbs_${firstn}_${nc}.com
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
        ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc > /dev/null &" >> $SCRATCH/pid_ssh.log
      else
        ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc > /dev/null &" >> $SCRATCH/pid_ssh.log
      endif
      set pid = $rh
    else
      if ( $sym == "H" || $sym == "HP" ) then
        ${bin_dir}/mult_hreconstruct_n.com $firstn $lastn $start $nc &
      else
        ${bin_dir}/mult_reconstruct_n.com $firstn $lastn $start $nc &
      endif
      set pid = $!
    endif

    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "SSH" ) then
      echo $pid ${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >> $SCRATCH/pid_temp.log
    endif
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
      echo ${data_input}_mult_reconstruct_r${nc}_n${firstn}.log $pid >> $SCRATCH/pid_temp.log
    endif
    @ nc++
  end

#  sleep 2

  if ( $lastn == $last ) then
    goto submission_done
  endif
  @ firstn = $firstn + $incr
  @ lastn = $lastn + $incr
  if ( $firstn >= $last ) set firstn = $last
  if ( $lastn >= $last ) set lastn = $last
end

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

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
  cat $SCRATCH/pid_ssh.log | awk '{print " "$2}' > $SCRATCH/pid_ssh_temp.log
  paste -d " " $SCRATCH/pid_ssh_temp.log $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
else
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
endif
\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
\rm $SCRATCH/pid_ssh_temp.log >& /dev/null

set nc = 1
while ( $nc <= $nclass )

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  set firstn = $first
  @ lastn = $first + $incr - 1
  if ( $lastn >= $last ) set lastn = $last

checkdone_rec:

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

#  sleep 5
  while ( $firstn <= $last )

    grep --binary-files=text finished $SCRATCH/${data_input}_mult_reconstruct_r${nc}_n${firstn}.log >& /dev/null

    if ($status) then
      sleep 5
      goto checkdone_rec
    endif

    if ($status) goto checkdone_rec

    echo "Cycle "$start": reconstruction for particles $firstn to $lastn, ref $nc, finished "`date` >> frealign.log
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" ) then
      set pid = `grep "${data_input}_mult_reconstruct_r${nc}_n${firstn}.log" $SCRATCH/pid.log | awk '{print $1}'`
      cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
      if ($pid != "") grep -v "${data_input}_mult_reconstruct_r${nc}_n${firstn}.log" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
      \rm $SCRATCH/pid_temp.log >& /dev/null
    endif

    @ firstn = $firstn + $incr
    @ lastn = $lastn + $incr
    if ( $lastn >= $last ) set lastn = $last
  end

  @ nc++

end

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set pid = `grep "${data_input}_mult_reconstruct_${first}_${last}" $SCRATCH/pid.log | awk '{print $1}'`
  cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
  if ($pid != "") grep -v "${data_input}_mult_reconstruct_${first}_${last}" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
  \rm $SCRATCH/pid_temp.log >& /dev/null
endif

\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
set nc = 1
while ( $nc <= $nclass )

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  echo "Cycle "$start": merging 3D dump files for ref $nc "`date` >> frealign.log
  echo $npart > $SCRATCH/merge_3d_r${nc}.in
  set firstn = $first
  @ lastn = $first + $incr - 1
  while ( $firstn <= $last )

    echo $SCRATCH/${data_input}_${start}_r${nc}_n${firstn}.${extension} >> $SCRATCH/merge_3d_r${nc}.in

    @ firstn = $firstn + $incr
    @ lastn = $lastn + $incr
    if ( $lastn >= $last ) set lastn = $last

  end

  echo $SCRATCH/${data_input}_${start}_r${nc}.res >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}.${extension} >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}_weights >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}_map1.${extension} >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}_map2.${extension} >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}_phasediffs >> $SCRATCH/merge_3d_r${nc}.in
  echo $SCRATCH/${data_input}_${start}_r${nc}_pointspread >> $SCRATCH/merge_3d_r${nc}.in

  cp $SCRATCH/${data_input}_${start}_r${nc}_n${first}.res $SCRATCH/${data_input}_${start}_r${nc}.res

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
    set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N m${nc} ${bin_dir}/merge_3d.com ${nc} | awk '{print $3}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
    set pid = `sbatch $stn -D $working_directory -J m${nc} -o /dev/null ${bin_dir}/merge_3d.com ${nc} | awk '{print $4}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
#    set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/merge_3d.com ${nc} | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
    set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/merge_3d.com ${nc} | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then 
    echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
    echo ${bin_dir}/merge_3d.com ${nc} \& >> $SCRATCH/slurm.com
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
    echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${nc}.com
    echo cd $working_directory >> $SCRATCH/pbs_${nc}.com
    echo ${bin_dir}/merge_3d.com ${nc} >> $SCRATCH/pbs_${nc}.com
    set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${nc}.com`
    set pid = `echo $pid | awk -F. '{print $1}'`
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
    echo "executable = ${bin_dir}/merge_3d.com"  > $SCRATCH/condor.com
    echo "universe = vanilla"                   >> $SCRATCH/condor.com
    echo "request_memory = ${mem_small} GB"     >> $SCRATCH/condor.com
    echo "arguments = $nc" >> $SCRATCH/condor.com
    echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
    echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
    echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
    echo "queue" >> $SCRATCH/condor.com
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
    if ( $nh > $nhosts ) set nh = 1
    set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
    @ nh++
    ssh -f $rh "cd $working_directory; nohup ${bin_dir}/merge_3d.com ${nc} > /dev/null &" >> $SCRATCH/pid_ssh.log
    set pid = $rh
  else
    ${bin_dir}/merge_3d.com ${nc} &
    set pid = $!
  endif

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "SSH" ) then
    echo $pid merge_3d_r${nc}.log >> $SCRATCH/pid_temp.log
  endif
  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
    echo merge_3d_r${nc}.log $pid >> $SCRATCH/pid_temp.log
  endif
  @ nc++

end

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  chmod +x $SCRATCH/slurm.com
  set incs = `grep com $SCRATCH/slurm.com | wc -l`
  set nnode = `echo $incs | awk '{print int(($1+15)/16)}'`
  set pid = `sbatch  -D $working_directory -J merge -N $nnode -n $incs -p normal -o $SCRATCH/frealign.o%j -t 12:00:00 ${bin_dir}/launcher.slurm`
  set pid = `echo $pid | awk '{ print $NF }'`
  echo $pid merge_3d.log >> $SCRATCH/pid_temp.log
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  set pid = `condor_submit $SCRATCH/condor.com | tail -1`
  set pid = `echo $pid | awk '{print $6}' | awk -F. '{print $1}'`
  set i = 0
  foreach nc (`grep arguments $SCRATCH/condor.com | awk '{print $3}'`)
    echo ${pid}.${i} merge_3d_r${nc}.log >> $SCRATCH/pid_temp.log
    @ i++
  end
endif

sleep 10

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
  cat $SCRATCH/pid_ssh.log | awk '{print " "$2}' > $SCRATCH/pid_ssh_temp.log
  paste -d " " $SCRATCH/pid_ssh_temp.log $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
else
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
endif
\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
\rm $SCRATCH/pid_ssh_temp.log >& /dev/null

set nc = 1
while ( $nc <= $nclass )
  
checkdonem:

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  sleep 2

  grep --binary-files=text "merge_3d.com finished" $SCRATCH/merge_3d_r${nc}.log >& /dev/null
  if ($status) goto checkdonem

  ls $SCRATCH/${data_input}_${start}_r${nc}.${extension} >& /dev/null
  if ($status) goto checkdonem

  ls $SCRATCH/${data_input}_${start}_r${nc}.res >& /dev/null
  if ($status) goto checkdonem

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE") then
    set pid = `grep "merge_3d_r${nc}.log" $SCRATCH/pid.log | awk '{print $1}'`
    cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
    if ($pid != "") grep -v "merge_3d_r${nc}.log" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
    \rm $SCRATCH/pid_temp.log >& /dev/null
  endif

  sleep 1

  echo "Cycle "$start": merging for ref $nc finished "`date` >> frealign.log
  cat $SCRATCH/${data_input}_${start}_r${nc}.res >> ${working_directory}/${data_input}_${start}_r${nc}.par
  mv $SCRATCH/${data_input}_${start}_r${nc}.${extension} ${working_directory}/${data_input}_${start}_r${nc}.${extension}

  set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
  if ( $form == "IMAGIC" ) mv $SCRATCH/${data_input}_${start}_r${nc}.img ${working_directory}/${data_input}_${start}_r${nc}.img

  @ nc++

end

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set pid = `grep "merge_3d" $SCRATCH/pid.log | awk '{print $1}'`
  cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
  if ($pid != "") grep -v "merge_3d" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
  \rm $SCRATCH/pid_temp.log >& /dev/null
endif

# Clean up

@ prev = $start - 1

if ( ! $no_delete ) then
  set nc = 1
  while ( $nc <= $nclass )

    \rm $SCRATCH/${data_input}_${start}_r${nc}.res >& /dev/null
    \rm $SCRATCH/${data_input}_${start}_r${nc}.par >& /dev/null
    set firstn = $first
    @ lastn = $first + $incr - 1
    while ( $firstn <= $last )

      \rm $SCRATCH/${data_input}_${start}_r${nc}_n${firstn}.${extension} >& /dev/null
      \rm $SCRATCH/${data_input}_${start}_r${nc}_n${firstn}.res >& /dev/null

      @ firstn = $firstn + $incr
      @ lastn = $lastn + $incr
      if ( $lastn >= $last ) set lastn = $last

    end

    \rm $SCRATCH/${data_input}_${prev}_r${nc}_map1.${extension} >& /dev/null
    \rm $SCRATCH/${data_input}_${prev}_r${nc}_map2.${extension} >& /dev/null
    \rm $SCRATCH/${data_input}_${prev}_r${nc}_weights >& /dev/null
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
  \rm $SCRATCH/pbs_*.com* >& /dev/null
endif

echo "Cycle "$start": reconstruction done." >> frealign.log

grep reconstruct $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) then
    echo "monitor term" >> $SCRATCH/pid.log
  endif
endif

date
