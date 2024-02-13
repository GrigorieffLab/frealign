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

if ( ${1} != "" ) then
  set start = ${1}
else
  set start		= `grep start_process $SCRATCH/mparameters_run | awk '{print $2}'`
endif
set end			= `grep end_process $SCRATCH/mparameters_run | awk '{print $2}'`
set first		= `grep first_particle $SCRATCH/mparameters_run | awk '{print $2}'`
set last		= `grep last_particle $SCRATCH/mparameters_run | awk '{print $2}'`
set data_input		= `grep data_input $SCRATCH/mparameters_run | awk '{print $2}'`
set occ_helical		= `grep occ_helical $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $occ_helical == "" ) then
  set occ_helical	= F
endif
set bin_dir		= `grep frealign_bin_dir $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $bin_dir == "" ) then
  set bin_dir		= `which frealign_v9.exe`
  set bin_dir		= ${bin_dir:h}
endif
set cluster_type	= `cat $SCRATCH/cluster_type.log`
if ( $status || $cluster_type == "" ) then
  set cluster_type	= `grep cluster_type $SCRATCH/mparameters_run | awk '{print $2}'`
endif
set mem_per_cpu		= `grep mem_per_cpu $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $mem_per_cpu == "" ) then
  set mem_per_cpu	= 4096
endif
set stn			= `grep qsub_string_ref mparameters | awk -F\" '{print $2}'`
set no_delete		= `grep delete_scratch $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $no_delete == "F" ) then
  set no_delete = 1
else
  set no_delete = 0
endif

set raw_images		= `grep raw_images_ref $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $raw_images == "" ) set raw_images = `grep raw_images $SCRATCH/mparameters_run | awk '{print $2}'`
set raw_images = `echo ${raw_images:r}`
set extension = `ls $raw_images.* | head -1`
if ( $extension == "" ) then
  set raw_images = ${working_directory}/${raw_images}
  set extension = `ls $raw_images.* | head -1`
endif
set extension = `echo ${extension:e}`
set nx = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text NX | awk '{print $4}'`
set mem_big = `echo $nx | awk '{print int(10 * $1^3 * 4 * 66 /1024^3 + 1)/10}'`
if ( `echo $mem_big | awk '{print int(1024 * $1)}'` > $mem_per_cpu ) then
  set mem_big = `echo $nx | awk '{print int(10 * $1^3 * 4 * 3 /1024^3 + 1)/10}'`
endif
set mem_big = `echo $mem_big | awk '{if ($1 < 1) {print 1} else {print $1} }'`
set mem_small = `echo $nx | awk '{print int(10 * $1^3 * 4 * 3 /1024^3 + 1)/10}'`
set mem_small = `echo $mem_small | awk '{if ($1 < 1) {print 1} else {print $1} }'`

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then
 set mem_big = `echo $mem_big | awk '{print int($1 * 1024)}'`
 set mem_small = `echo $mem_small | awk '{print int($1 * 1024)}'`
endif

if ( $first == "" || $last == "" ) then
  set first = 1
  set last = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text NX | awk '{print $6}'`
endif

sleep 2
echo "monitor refine" >> $SCRATCH/pid.log

# Check what type of parameter file is available
set pix                 = `grep pix_size $SCRATCH/mparameters_run | awk '{print $2}'`
set dstep               = `grep dstep $SCRATCH/mparameters_run | awk '{print $2}'`
set sym                 = `grep Symmetry $SCRATCH/mparameters_run | awk '{print $2}'`
@ prev = $start - 1
set nc = `grep --binary-files=text -v C ${working_directory}/${data_input}_${prev}_r1.par | head -1 | wc | awk '{print $2}'`
if ( $nc == 4 ) then
  echo "Found startup parameter file with film#, defocus 1, defocus 2, astig. angle." >> frealign.log
  set mag = `echo $dstep $pix | awk '{print int(10000 * $1 / $2)}'`
  echo "Setting magnification to "${mag}"." >> frealign.log
  mv ${data_input}_${prev}_r1.par ${data_input}_${prev}.par
  echo "C Startup parameters" > ${data_input}_${prev}_r1.par
  echo "C           PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE" >> ${data_input}_${prev}_r1.par
  if ( $sym == "H" || $sym == "HP" ) then
    grep --binary-files=text -v C ${data_input}_${prev}.par | awk '{printf "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%7.2f%11d%11.4f%8.2f%8.2f\n",NR,360.0*rand(),90.0,360.0*rand(),0.0,0.0,'${mag}',$1,$2,$3,$4,100.0,0,0.5,100.0,0.0 }' >> ${data_input}_${prev}_r1.par
  else
    grep --binary-files=text -v C ${data_input}_${prev}.par | awk '{printf "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%7.2f%11d%11.4f%8.2f%8.2f\n",NR,360.0*rand(),360.0*rand(),360.0*rand(),0.0,0.0,'${mag}',$1,$2,$3,$4,100.0,0,0.5,100.0,0.0 }' >> ${data_input}_${prev}_r1.par
  endif
endif

set nc = `grep --binary-files=text -v C ${working_directory}/${data_input}_${prev}_r1.par | head -1 | wc | awk '{print $2}'`
if ( $nc < 12 ) then
  echo "ERROR: Something wrong with parameter file." >> frealign.log
  echo "If input file contained film#, defocus 1, defocus 2, astig. angle," >> frealign.log
  echo "please make sure that numbers are separated by spaces, not commas." >> frealign.log
  echo "monitor term" >> $SCRATCH/pid.log
endif

mainloop:

set restart		= `grep restart_after_crash $SCRATCH/mparameters_run | awk '{print $2}'`
set nclass		= `grep nclasses $SCRATCH/mparameters_run | awk '{print $2}'`
set pmask		= `grep parameter_mask $SCRATCH/mparameters_run | awk -F\" '{print $2}'`
set pmask = `echo $pmask "1 1 1 1 1" | awk '{print $1,$2,$3,$4,$5}'`
set pshift = `echo $pmask | awk '{print 0,0,0,$4,$5}'`
set pangle = `echo $pmask | awk '{print $1,$2,$3,0,0}'`
@ prev = $start - 1

# Split into N classes if desired
set nc = `ls ${data_input}_${prev}_r[0-9].par ${data_input}_${prev}_r[0-9][0-9].par ${data_input}_${prev}_r[0-9][0-9][0-9].par | wc -l`
if ($nc != $nclass) then
  ${bin_dir}/split_classes.com $start
endif

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

set nc = `ls ${data_input}_${prev}_r[0-9].par ${data_input}_${prev}_r[0-9][0-9].par ${data_input}_${prev}_r[0-9][0-9][0-9].par | wc -l`
set i = `ls ${data_input}_${prev}_r[0-9].${extension} ${data_input}_${prev}_r[0-9][0-9].${extension} ${data_input}_${prev}_r[0-9][0-9][0-9].${extension} | wc -l`
if ($nc != $i) then
  ${bin_dir}/mult_reconstruct.com $first $last $prev

  if ($nclass <= 1) goto skip_occ3

  echo $nclass > $SCRATCH/calc_occ.in
  echo 0.0 >> $SCRATCH/calc_occ.in
  set nc = 1
  while ($nc <= $nclass)
    echo ${data_input}_${prev}_r${nc}.par >> $SCRATCH/calc_occ.in
    @ nc++
  end
  set nc = 1
  while ($nc <= $nclass)
    echo ${data_input}_${prev}_r${nc}.par >> $SCRATCH/calc_occ.in
    @ nc++
  end

  if ( $sym == "H" || $sym == "HP" ) then
    if ( $occ_helical == "T" ) then
      time ${bin_dir}/calc_occ_helical.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
    else
      time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
    endif
  else
    time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
  endif

skip_occ3:

endif

cp mparameters $SCRATCH/mparameters_run
set nproc		= `grep nprocessor_ref $SCRATCH/mparameters_run | awk '{print $2}'`
set itmax		= `grep ITMAX $SCRATCH/mparameters_run | awk '{print $2}'`
set refineshiftinc	= `grep refineshiftinc $SCRATCH/mparameters_run | awk '{print $2}'`
set refineangleinc	= `grep refineangleinc $SCRATCH/mparameters_run | awk '{print $2}'`
set mode		= `grep MODE $SCRATCH/mparameters_run | awk '{print $2}'`
set sym			= `grep Symmetry $SCRATCH/mparameters_run | awk '{print $2}'`

set shiftiter		= `expr ${start} % ${refineshiftinc}`
set angleiter		= `expr ${start} % ${refineangleinc}`

if ( $nclass < 2 ) then
  set shiftiter = 0
  set angleiter = 0
endif

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

#if ( $mode != 1 ) then
#  if (${shiftiter} == 0) then
#    if (${angleiter} == 0) then
#      set incr = `expr ${incr} / ${itmax} \* 3`
#    endif
#  endif
#endif

\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
# apply external mask

set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
set fm = "M"
if ( $form == "SPIDER" ) set fm = "S"
if ( $form == "IMAGIC" ) set fm = "I"

set mask_file		= `grep mask_file $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $mask_file == "" ) then
  set nc = 1
  @ prev2 = $prev - 1
  while ( $nc <= $nclass )
    \rm $SCRATCH/${data_input}_${prev2}_r${nc}.${extension} >& /dev/null
    \rm $SCRATCH/${data_input}_${prev}_r${nc}.${extension} >& /dev/null
    ln -s ${working_directory}/${data_input}_${prev}_r${nc}.${extension} $SCRATCH/${data_input}_${prev}_r${nc}.${extension}
    if ( $form == "IMAGIC" ) ln -s ${working_directory}/${data_input}_${prev}_r${nc}.${extension} $SCRATCH/${data_input}_${prev}_r${nc}.${extension}
    @ nc++
  end
else
  set mask_edge			= `grep mask_edge $SCRATCH/mparameters_run | awk '{print $2}'`
  set mask_outside_weight	= `grep mask_outside_weight $SCRATCH/mparameters_run | awk '{print $2}'`
  set mask_filt_res		= `grep mask_filt_res $SCRATCH/mparameters_run | awk '{print $2}'`
  set mask_filt_edge		= `grep mask_filt_edge $SCRATCH/mparameters_run | awk '{print $2}'`
  set nc = 1
  @ prev2 = $prev - 1
  while ( $nc <= $nclass )
    grep kill $SCRATCH/pid.log >& /dev/null
    if ( ! $status ) exit

    \rm $SCRATCH/${data_input}_${prev2}_r${nc}.${extension} >& /dev/null

    echo $fm > $SCRATCH/apply_mask_r${nc}.in
    echo ${working_directory}/${data_input}_${prev}_r${nc}.${extension} >> $SCRATCH/apply_mask_r${nc}.in
    echo $pix >> $SCRATCH/apply_mask_r${nc}.in
    echo $mask_file >> $SCRATCH/apply_mask_r${nc}.in
    echo $mask_edge >> $SCRATCH/apply_mask_r${nc}.in
    echo $mask_outside_weight >> $SCRATCH/apply_mask_r${nc}.in
    if ( `echo $mask_filt_res | awk '{if ($1 == 0) print 0}'` == "0" ) then
      echo 0 >> $SCRATCH/apply_mask_r${nc}.in
    else
      echo 2 >> $SCRATCH/apply_mask_r${nc}.in
#      echo $mask_filt_res | awk '{if ($1 != 0) print $1}' >> $SCRATCH/apply_mask_r${nc}.in
      echo $mask_filt_res >> $SCRATCH/apply_mask_r${nc}.in
      echo $mask_filt_edge >> $SCRATCH/apply_mask_r${nc}.in
    endif
    echo $SCRATCH/${data_input}_${prev}_r${nc}.${extension} >> $SCRATCH/apply_mask_r${nc}.in

    echo "Cycle "$start": applying mask file to ref $nc "`date` >> frealign.log

    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
      set pid = `qsub $stn -l mem_free=${mem_small}G -cwd -r y -V -e /dev/null -o /dev/null -N m${nc} ${bin_dir}/apply_mask.com ${nc} | awk '{print $3}'`
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
      set pid = `sbatch $stn -D $working_directory -J m${nc} -o /dev/null ${bin_dir}/apply_mask.com ${nc} | awk '{print $4}'`
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then
#      set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_small}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/apply_mask.com ${nc} | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
      set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/apply_mask.com ${nc} | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
      echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
      echo ${bin_dir}/apply_mask.com ${nc} \& >> $SCRATCH/slurm.com
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
      echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${nc}.com
      echo cd $working_directory >> $SCRATCH/pbs_${nc}.com
      echo ${bin_dir}/apply_mask.com ${nc} >> $SCRATCH/pbs_${nc}.com
      set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${nc}.com`
      set pid = `echo $pid | awk -F. '{print $1}'`
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
      echo "executable = ${bin_dir}/apply_mask.com" > $SCRATCH/condor.com
      echo "universe = vanilla"                    >> $SCRATCH/condor.com
      echo "request_memory = ${mem_small} GB"      >> $SCRATCH/condor.com
      echo "arguments = $nc" >> $SCRATCH/condor.com
      echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
      echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
      echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
      echo "queue" >> $SCRATCH/condor.com
    else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
      if ( $nh > $nhosts ) set nh = 1
      set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
      @ nh++
      ssh -f $rh "cd $working_directory; nohup ${bin_dir}/apply_mask.com ${nc} > /dev/null &" >> $SCRATCH/pid_ssh.log
      set pid = $rh
    else
      ${bin_dir}/apply_mask.com ${nc} &
      set pid = $!
    endif
    
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "SSH" ) then
      echo $pid apply_mask_r${nc}.log >> $SCRATCH/pid_temp.log
    endif
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
      echo apply_mask_r${nc}.log $pid >> $SCRATCH/pid_temp.log
    endif

    @ nc++
  end

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
    chmod +x $SCRATCH/slurm.com
    set incs = `grep com $SCRATCH/slurm.com | wc -l`
    set nnode = `echo $incs | awk '{print int(($1+15)/16)}'`
    set pid = `sbatch  -D $working_directory -J merge -N $nnode -n $incs -p normal -o $SCRATCH/frealign.o%j -t 12:00:00 ${bin_dir}/launcher.slurm`
    set pid = `echo $pid | awk '{ print $NF }'`
    echo $pid apply_mask.log >> $SCRATCH/pid_temp.log
  endif
  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then   
    set pid = `condor_submit $SCRATCH/condor.com | tail -1`
    set pid = `echo $pid | awk '{print $6}' | awk -F. '{print $1}'`
    set i = 0
    foreach nc (`grep arguments $SCRATCH/condor.com | awk '{print $3}'`)
      echo ${pid}.${i} apply_mask_r${nc}.log >> $SCRATCH/pid_temp.log
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

    grep --binary-files=text "apply_mask.com finished" $SCRATCH/apply_mask_r${nc}.log >& /dev/null
    if ($status) goto checkdonem

    ls $SCRATCH/${data_input}_${prev}_r${nc}.${extension} >& /dev/null
    if ($status) goto checkdonem

    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE") then
      set pid = `grep "apply_mask_r${nc}.log" $SCRATCH/pid.log | awk '{print $1}'`
      cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
      if ($pid != "") grep -v "apply_mask_r${nc}.log" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
      \rm $SCRATCH/pid_temp.log >& /dev/null
    endif

    sleep 1

    echo "Cycle "$start": finished applying mask to ref $nc "`date` >> frealign.log
    cat $SCRATCH/${data_input}_${start}_r${nc}.res >> ${working_directory}/${data_input}_${start}_r${nc}.par
    mv $SCRATCH/${data_input}_${start}_r${nc}.${extension} ${working_directory}/${data_input}_${start}_r${nc}.${extension}

    set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
    if ( $form == "IMAGIC" ) mv $SCRATCH/${data_input}_${start}_r${nc}.img ${working_directory}/${data_input}_${start}_r${nc}.img

    @ nc++

  end

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit 

  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
    set pid = `grep "apply_mask" $SCRATCH/pid.log | awk '{print $1}'`
    cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
    if ($pid != "") grep -v "apply_mask" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
    \rm $SCRATCH/pid_temp.log >& /dev/null
  endif

endif

# queue submission loop

\rm $SCRATCH/pid_temp.log >& /dev/null
\rm $SCRATCH/pid_ssh.log >& /dev/null
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  echo \#\!\/bin\/csh -f > $SCRATCH/slurm.com
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  if ( $sym == "H" || $sym == "HP" ) then
    echo "executable = ${bin_dir}/mult_hrefine_n.com" > $SCRATCH/condor.com
  else
    echo "executable = ${bin_dir}/mult_refine_n.com"  > $SCRATCH/condor.com
  endif
  echo "universe = vanilla"                          >> $SCRATCH/condor.com
  echo "request_memory = ${mem_big} GB"              >> $SCRATCH/condor.com
endif

set firstn = $first
@ lastn = $first + $incr - 1
if ( $lastn >= $last ) set lastn = $last
set nc = 1
if ( $restart == "T" ) then
  if ( -e $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} || -e ${data_input}_${start}_r${nc}.par ) then
    echo "Attempting to recovering from crashed run..." >> frealign.log
  else
    echo "Restarting crashed job failed: no intermediate files found."
    exit
  endif
endif

subloop:
set firstn = $first
@ lastn = $first + $incr - 1
if ( $lastn >= $last ) set lastn = $last

while ( $lastn <= $last )
#  echo "Cycle "$start": refining particles" $firstn to $lastn on `date` >> frealign.log

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  set nc = 1
  while ( $nc <= $nclass )

    set fst = $firstn
    if ( $restart == "T" ) then
      if ( -e $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} ) then
loop_recovery:
        set fst = `grep --binary-files=text -v C $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} | tail -1 | awk '{print $1}'`
        if ( $fst != "" ) then
          @ fst++
          if ( -e $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_${lastn} ) then
            grep --binary-files=text -v C $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_${lastn} >> $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn}
            \rm $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_${lastn} >& /dev/null
            goto loop_recovery
          endif
        else
          set fst = $firstn
          \rm $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} >& /dev/null
        endif
      else if ( -e ${data_input}_${start}_r${nc}.par ) then
        set fst = `grep --binary-files=text -v C ${data_input}_${start}_r${nc}.par | tail -1 | awk '{print $1}'`
        if ( $fst > $firstn ) then
          echo "Cycle "$start": skipping particles" $firstn to $lastn", class "${nc}" "`date` >> frealign.log
          goto skipjob
        endif
      endif
    endif

      if ( $lastn >= $fst ) then
      echo "Cycle "$start": refining particles" $fst to $lastn", class "${nc}" "`date` >> frealign.log

    if (${shiftiter} == 0) then
      if (${angleiter} == 0) then
        if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask | awk '{print $3}'`
          else
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask | awk '{print $3}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask | awk '{print $4}'`
          else
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask | awk '{print $4}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
          if ( $sym == "H" || $sym == "HP" ) then
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          else
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask \& >> $SCRATCH/slurm.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask \& >> $SCRATCH/slurm.com
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
          echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${fst}_${nc}.com
          echo cd $working_directory >> $SCRATCH/pbs_${fst}_${nc}.com
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask >> $SCRATCH/pbs_${fst}_${nc}.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask >> $SCRATCH/pbs_${fst}_${nc}.com
          endif
          set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${fst}_${nc}.com`
          set pid = `echo $pid | awk -F. '{print $1}'`
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
          echo "arguments = $fst $lastn $start $nc $mode $pmask" >> $SCRATCH/condor.com
          echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
          echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
          echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
          echo "queue" >> $SCRATCH/condor.com
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
          if ( $nh > $nhosts ) set nh = 1
          set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
          @ nh++
          if ( $sym == "H" || $sym == "HP" ) then
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask > /dev/null &" >> $SCRATCH/pid_ssh.log
          else
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask > /dev/null &" >> $SCRATCH/pid_ssh.log
          endif
          set pid = $rh
        else
          if ( $sym == "H" || $sym == "HP" ) then
            ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc $mode $pmask &
          else
            ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc $mode $pmask &
          endif
          set pid = $!
        endif
      else if (${angleiter} != 0) then
        if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift | awk '{print $3}'`
          else
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift | awk '{print $3}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift | awk '{print $4}'`
          else
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift | awk '{print $4}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
          if ( $sym == "H" || $sym == "HP" ) then
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          else
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift \& >> $SCRATCH/slurm.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift \& >> $SCRATCH/slurm.com
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
          echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${fst}_${nc}.com
          echo cd $working_directory >> $SCRATCH/pbs_${fst}_${nc}.com
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift >> $SCRATCH/pbs_${fst}_${nc}.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift >> $SCRATCH/pbs_${fst}_${nc}.com
          endif
          set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${fst}_${nc}.com`
          set pid = `echo $pid | awk -F. '{print $1}'`
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
          echo "arguments = $fst $lastn $start $nc 1 $pshift" >> $SCRATCH/condor.com
          echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
          echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
          echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
          echo "queue" >> $SCRATCH/condor.com
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
          if ( $nh > $nhosts ) set nh = 1
          set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
          @ nh++
          if ( $sym == "H" || $sym == "HP" ) then
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift > /dev/null &" >> $SCRATCH/pid_ssh.log
          else
            ssh -f $rh "cd $working_directory; nobup ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift > /dev/null &" >> $SCRATCH/pid_ssh.log
          endif
          set pid = $rh
        else
          if ( $sym == "H" || $sym == "HP" ) then
            ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pshift &
          else
            ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pshift &
          endif
          set pid = $!
        endif
      endif
    else if (${shiftiter} != 0) then
      if (${angleiter} == 0) then
        if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle | awk '{print $3}'`
          else
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle | awk '{print $3}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle | awk '{print $4}'`
          else
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle | awk '{print $4}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
          if ( $sym == "H" || $sym == "HP" ) then
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          else
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle \& >> $SCRATCH/slurm.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle \& >> $SCRATCH/slurm.com
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
          echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${fst}_${nc}.com
          echo cd $working_directory >> $SCRATCH/pbs_${fst}_${nc}.com
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle >> $SCRATCH/pbs_${fst}_${nc}.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle >> $SCRATCH/pbs_${fst}_${nc}.com
          endif
          set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${fst}_${nc}.com`
          set pid = `echo $pid | awk -F. '{print $1}'`
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
          echo "arguments = $fst $lastn $start $nc 1 $pangle" >> $SCRATCH/condor.com
          echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
          echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
          echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
          echo "queue" >> $SCRATCH/condor.com
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
          if ( $nh > $nhosts ) set nh = 1
          set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
          @ nh++
          if ( $sym == "H" || $sym == "HP" ) then
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle > /dev/null &" >> $SCRATCH/pid_ssh.log
          else
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle > /dev/null &" >> $SCRATCH/pid_ssh.log
          endif
          set pid = $rh
        else
          if ( $sym == "H" || $sym == "HP" ) then
            ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 $pangle &
          else
            ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 $pangle &
          endif
          set pid = $!
        endif
      else if (${angleiter} != 0) then
        if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk '{print $3}'`
          else
            set pid = `qsub $stn -l mem_free=${mem_big}G -cwd -r y -V -e /dev/null -o /dev/null -N n${nc}_${fst} ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk '{print $3}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
          if ( $sym == "H" || $sym == "HP" ) then
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk '{print $4}'`
          else
            set pid = `sbatch $stn -D $working_directory -J n${nc}_${fst} -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk '{print $4}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then 
          if ( $sym == "H" || $sym == "HP" ) then
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          else
#            set pid = `bsub -W 12:00 -q short -R "select[model!=Opteron2216]" -R "rusage[mem=${mem_big}]" -cwd  "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
            set pid = `bsub $stn -cwd "$working_directory" -r -e /dev/null -o /dev/null ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 | awk -F\< '{print $2}' | awk -F\> '{print $1}'`
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 \& >> $SCRATCH/slurm.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 \& >> $SCRATCH/slurm.com
          endif
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
          echo \#\!\/bin\/csh -f > $SCRATCH/pbs_${fst}_${nc}.com
          echo cd $working_directory >> $SCRATCH/pbs_${fst}_${nc}.com
          if ( $sym == "H" || $sym == "HP" ) then
            echo ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 >> $SCRATCH/pbs_${fst}_${nc}.com
          else
            echo ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 >> $SCRATCH/pbs_${fst}_${nc}.com
          endif
          set pid = `qsub $stn -l walltime=100:00:00,nodes=1:ppn=1 -o $SCRATCH -e $SCRATCH $SCRATCH/pbs_${fst}_${nc}.com`
          set pid = `echo $pid | awk -F. '{print $1}'`
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
          echo "arguments = $fst $lastn $start $nc 1 0 0 0 0 0" >> $SCRATCH/condor.com
          echo "output = $SCRATCH/condor.out" >> $SCRATCH/condor.com
          echo "error = $SCRATCH/condor.err" >> $SCRATCH/condor.com
          echo "log = $SCRATCH/condor.log" >> $SCRATCH/condor.com
          echo "queue" >> $SCRATCH/condor.com
        else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
          if ( $nh > $nhosts ) set nh = 1
          set rh = `awk 'FNR=='$nh'' $SCRATCH/hosts`
          @ nh++
          if ( $sym == "H" || $sym == "HP" ) then
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 > /dev/null &" >> $SCRATCH/pid_ssh.log
          else
            ssh -f $rh "cd $working_directory; nohup ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 > /dev/null &" >> $SCRATCH/pid_ssh.log
          endif
          set pid = $rh
        else
          if ( $sym == "H" || $sym == "HP" ) then
            ${bin_dir}/mult_hrefine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 &
          else
            ${bin_dir}/mult_refine_n.com $fst $lastn $start $nc 1 0 0 0 0 0 &
          endif
          set pid = $!
        endif
      endif
    endif
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "CONDOR" && `echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "SSH" ) then
      echo $pid ${data_input}_mult_refine_n_r${nc}.log_${fst}_${lastn} >> $SCRATCH/pid_temp.log
    endif
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
      echo ${data_input}_mult_refine_n_r${nc}.log_${fst}_${lastn} $pid >> $SCRATCH/pid_temp.log
    endif

    endif

skipjob:

    @ nc++
  end

#  sleep 2

  if ( $lastn == $last ) then
    goto alignment_done
  endif
  @ firstn = $firstn + $incr
  @ lastn = $lastn + $incr
  if ( $firstn >= $last ) set firstn = $last
  if ( $lastn >= $last ) set lastn = $last
end

alignment_done:

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  chmod +x $SCRATCH/slurm.com
  set incs = `grep com $SCRATCH/slurm.com | wc -l`
  set nnode = `echo $incs | awk '{print int(($1+15)/16)}'`
  set pid = `sbatch  -D $working_directory -J refine -N $nnode -n $incs -p normal -o $SCRATCH/frealign.o%j -t 12:00:00 ${bin_dir}/launcher.slurm`
  set pid = `echo $pid | awk '{ print $NF }'`
  echo $pid ${data_input}_mult_refine_n.log_${first}_${last} >> $SCRATCH/pid_temp.log
endif
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  set pid = `condor_submit $SCRATCH/condor.com | tail -1`
  set pid = `echo $pid | awk '{print $6}' | awk -F. '{print $1}'`
  set i = 0
  foreach fst (`grep arguments $SCRATCH/condor.com | awk '{print $3}'`)
    set lastn = `grep arguments scratch/condor.com | cat -n | awk '{if ($1 - 1 == '$i') print $5}'`
    set nc = `grep arguments scratch/condor.com | cat -n | awk '{if ($1 - 1 == '$i') print $7}'`
    echo ${pid}.${i} ${data_input}_mult_refine_n_r${nc}.log_${fst}_${lastn} >> $SCRATCH/pid_temp.log
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

if ( $restart != "T" ) then
  set nc = 1
  while ( $nc <= $nclass )
    \rm ${data_input}_${start}_r${nc}.par >& /dev/null
    @ nc++
  end
endif

set nc = 1
while ( $nc <= $nclass )

  set firstn = $first
  if ( $restart == "T" && -e ${data_input}_${start}_r${nc}.par ) then
    set firstn = `grep --binary-files=text -v C ${data_input}_${start}_r${nc}.par | tail -1 | awk '{print $1}'`
    @ firstn++
  endif
  @ lastn = $first + $incr - 1
  if ( $lastn >= $last ) set lastn = $last

checkdone:

  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) exit

  sleep 5
  while ( $firstn <= $last )
    set fst = `grep --binary-files=text -v C $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} | tail -1 | awk '{print $1}'`
    if ( $lastn > $fst ) then
      @ fst++
      grep --binary-files=text "overall score" $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_$lastn >& /dev/null
      if ( ! $status) then
        grep --binary-files=text -v C $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_$lastn >> $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_$lastn
        grep --binary-files=text "overall score" $SCRATCH/${data_input}_${start}_r${nc}.par_${fst}_$lastn >> $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_$lastn
      endif
    else
      set fst = $firstn
    endif

    grep --binary-files=text "overall score" $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_$lastn >& /dev/null

    if ($status) goto checkdone

    echo "Cycle "$start": particles $firstn to $lastn, ref $nc, finished "`date` >> frealign.log
    if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` != "STAMPEDE" ) then
      set pid = `grep "${data_input}_mult_refine_n_r${nc}.log_${fst}_${lastn}" $SCRATCH/pid.log | awk '{print $1}'`
      cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
      if ($pid != "") grep -v "${data_input}_mult_refine_n_r${nc}.log_${fst}_${lastn}" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
      \rm $SCRATCH/pid_temp.log >& /dev/null
    endif
    if ($firstn == $first ) head -68 $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} | grep --binary-files=text C >> ${data_input}_${start}_r${nc}.par

    grep --binary-files=text -v C --binary-files=text $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} >> ${data_input}_${start}_r${nc}.par
    if ( ! $no_delete ) \rm $SCRATCH/${data_input}_${start}_r${nc}.par_${firstn}_${lastn} >& /dev/null
    if ( ! $no_delete ) \rm $SCRATCH/${data_input}_${start}_r${nc}.shft_${firstn}_${lastn} >& /dev/null

    @ firstn = $firstn + $incr
    @ lastn = $lastn + $incr
    if ( $lastn >= $last ) set lastn = $last
  end

collect_done:

  @ nc++

end

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  set pid = `grep "${data_input}_mult_refine_n.log_${first}_${last}" $SCRATCH/pid.log | awk '{print $1}'`
  cp $SCRATCH/pid.log $SCRATCH/pid_temp.log
  if ($pid != "") grep -v "${data_input}_mult_refine_n.log_${first}_${last}" $SCRATCH/pid_temp.log > $SCRATCH/pid.log
  \rm $SCRATCH/pid_temp.log >& /dev/null
endif

if ($nclass <= 1) goto skip_occ

echo $nclass > $SCRATCH/calc_occ.in
echo 1.0 >> $SCRATCH/calc_occ.in
set nc = 1
while ($nc <= $nclass)
  echo ${data_input}_${start}_r${nc}.par >> $SCRATCH/calc_occ.in
  @ nc++
end
set nc = 1
while ($nc <= $nclass)
  echo ${data_input}_${start}_r${nc}.par >> $SCRATCH/calc_occ.in
  @ nc++
end

if ( $sym == "H" || $sym == "HP" ) then
  if ( $occ_helical == "T" ) then
    time ${bin_dir}/calc_occ_helical.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
  else
    time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
  endif
else
  time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
endif

skip_occ:

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

reconstruct:
echo "Cycle "$start": calculating 3D structures...." >> frealign.log

${bin_dir}/mult_reconstruct.com $first $last $start

# sleep 20

if ($nclass <= 1) goto skip_occ2

echo $nclass > $SCRATCH/calc_occ.in
echo 0.0 >> $SCRATCH/calc_occ.in
set nc = 1
while ($nc <= $nclass)
  echo ${data_input}_${start}_r${nc}.par >> $SCRATCH/calc_occ.in
  @ nc++
end
set nc = 1
while ($nc <= $nclass)
  echo ${data_input}_${start}_r${nc}.par >> $SCRATCH/calc_occ.in
  @ nc++
end

if ( $sym == "H" || $sym == "HP" ) then
  if ( $occ_helical == "T" ) then
    time ${bin_dir}/calc_occ_helical.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
  else
    time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
  endif
else
  time ${bin_dir}/calc_occ.exe < $SCRATCH/calc_occ.in >& $SCRATCH/calc_occ.log
endif

skip_occ2:

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  \rm $SCRATCH/slurm.com >& /dev/null
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  \rm $SCRATCH/condor.* >& /dev/null
endif

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
  \rm $SCRATCH/pbs_*.com* >& /dev/null
endif

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) exit

cp mparameters $SCRATCH/mparameters_run
set end		= `grep end_process $SCRATCH/mparameters_run | awk '{print $2}'`
set mem_per_cpu		= `grep mem_per_cpu $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $status || $mem_per_cpu == "" ) then
  set mem_per_cpu	= 4096
endif
set stn			= `grep qsub_string_ref mparameters | awk -F\" '{print $2}'`
set no_delete		= `grep delete_scratch $SCRATCH/mparameters_run | awk '{print $2}'`
if ( $no_delete == "F" ) then
  set no_delete = 1
else
  set no_delete = 0
endif

if ($start < $end ) then
  @ start++
  goto mainloop
endif

echo "Normal termination of refinement "`date` >> frealign.log
echo "monitor term" >> $SCRATCH/pid.log

date
