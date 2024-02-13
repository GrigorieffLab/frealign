#!/bin/csh -f
#
#   Control script to monitor Frealign jobs
#
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
set cluster_type	= `cat $SCRATCH/cluster_type.log`
set cpid		= ${1}
set restart		= `grep restart_after_crash $SCRATCH/mparameters_run | awk '{print $2}'`

wait_for_pid:

sleep 5
touch $SCRATCH/monitor_frealign.log
sleep 5
touch $SCRATCH/monitor_frealign.log

set maxtry = 10
if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") set maxtry = 60

monitorloop:

if ( ! -e $SCRATCH/pid.log ) goto wait_for_pid

grep term $SCRATCH/pid.log >& /dev/null
if ( ! $status ) then
  echo "Normal termination of frealign run "`date` >> frealign.log
  set restart = "F"
  goto terminate
endif

grep kill $SCRATCH/pid.log >& /dev/null
if ( ! $status ) then
  echo "Killing frealign run..." >> frealign.log
  set restart = "F"
  goto terminate
endif

ps ${cpid} >& /dev/null
if ( $status ) then
  sleep 5
  touch $SCRATCH/monitor_frealign.log
  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) then
    echo "Killing frealign run..." >> frealign.log
    set restart = "F"
    goto terminate
  endif
  grep term $SCRATCH/pid.log >& /dev/null
  if ( $status ) then
    grep "monitor reconstruct" $SCRATCH/pid.log >> /dev/null
    if ( $status ) then
      echo "Frealign run script crashed "`date` >> frealign.log
      echo "Terminating..." >> frealign.log
# Remove the following line if job should be restarted automatically
      set restart = "F"
    else
      echo "Normal termination of frealign run "`date` >> frealign.log
      set restart = "F"
    endif
    goto terminate
  endif
endif

foreach pid (`grep -v monitor $SCRATCH/pid.log | awk '{print $1}'`)

  set try = 1

testloop:
  grep kill $SCRATCH/pid.log >& /dev/null
  if ( ! $status ) then
    echo "Killing frealign run..." >> frealign.log
    set restart = "F"
    goto terminate
  endif
  if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
    qstat -j $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
    squeue -j $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then
    bjobs $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
    squeue -j $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
    qstat $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
    condor_q | grep $pid >& /dev/null
    set err = $status
  else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
    set rh = `grep " $pid " $SCRATCH/pid.log | awk '{print $3}'`
    ssh $rh ps $pid >& /dev/null
    set err = $status
  else
    ps $pid >& /dev/null
    set err = $status
  endif
  if ($err) then
    grep term $SCRATCH/pid.log >& /dev/null
    if ( ! $status ) then
      echo "Normal termination of frealign run "`date` >> frealign.log
      set restart = "F"
      goto terminate
    endif
    sleep 5
    touch $SCRATCH/monitor_frealign.log
    grep $pid $SCRATCH/pid.log >& /dev/null
    if ( ! $status ) sleep 5
    touch $SCRATCH/monitor_frealign.log
    grep $pid $SCRATCH/pid.log >& /dev/null    
    if ( ! $status ) then
      set logfile = `grep $pid $SCRATCH/pid.log | awk '{print $2}'`
      grep "Normal termination" $SCRATCH/${logfile} >& /dev/null
      if ($status) then
        if ( $try < $maxtry ) then
          @ try++
          goto testloop
        endif
        echo "Job $pid crashed "`date` >> frealign.log
        sleep 5
        touch $SCRATCH/monitor_frealign.log
        echo "Logfile "$SCRATCH/${logfile} >> frealign.log
        echo "Final lines:" >> frealign.log
        echo "" >> frealign.log
        tail $SCRATCH/${logfile} >> frealign.log
        echo "Terminating..." >> frealign.log
        goto terminate
      endif
    endif
  endif

end

sleep 2
touch $SCRATCH/monitor_frealign.log

goto monitorloop

terminate:

if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SGE") then
  kill `pstree -p ${cpid} | sed 's/(/\n(/g' | grep '(' | sed 's/(\(.*\)).*/\1/' | tr "\n" " "`
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  qdel `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SLURM") then
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  scancel `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "LSF") then
  kill `pstree -p ${cpid} | sed 's/(/\n(/g' | grep '(' | sed 's/(\(.*\)).*/\1/' | tr "\n" " "`
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  bkill `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "STAMPEDE") then
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  scancel `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "PBS") then
  kill `pstree -p ${cpid} | sed 's/(/\n(/g' | grep '(' | sed 's/(\(.*\)).*/\1/' | tr "\n" " "`
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  qdel `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "CONDOR") then
  kill `pstree -p ${cpid} | sed 's/(/\n(/g' | grep '(' | sed 's/(\(.*\)).*/\1/' | tr "\n" " "`
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  set pid = `tail -1 $SCRATCH/pid.log | awk -F. '{print $1}'`
  condor_rm $pid
else if (`echo ${cluster_type} | tr '[a-z]' '[A-Z]'` == "SSH") then
  ps -eo pid,ppid | grep ${cpid} | awk '{print $1}' > $SCRATCH/p.log
  foreach p (`cat $SCRATCH/p.log`)
    kill -9 `ps -eo pid,ppid | grep ${p} | grep frealign | awk '{print $1}'`
  end
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  foreach p (`grep -v monitor $SCRATCH/pid.log | awk '{print $1}'`)
    kill -9 `ps -eo pid,ppid | grep ${p} | awk '{print $1}'`
  end
  foreach h (`cat hosts | awk '{print $1}'`)
    ssh $h "killall frealign_v9.exe"
  end
  kill -9 ${cpid}
else
##  kill `pstree -p ${cpid} | sed 's/(/\n(/g' | grep '(' | sed 's/(\(.*\)).*/\1/' | tr "\n" " "`
##  kill -9 `grep -v monitor $SCRATCH/pid.log | awk '{print $1}'` >& /dev/null
  ps -eo pid,ppid | grep ${cpid} | awk '{print $1}' > $SCRATCH/p.log
  foreach p (`cat $SCRATCH/p.log`)
    kill -9 `ps -eo pid,ppid | grep ${p} | grep frealign | awk '{print $1}'`
  end
  cat $SCRATCH/pid_temp.log >> $SCRATCH/pid.log
  foreach p (`grep -v monitor $SCRATCH/pid.log | awk '{print $1}'`)
    kill -9 `ps -eo pid,ppid | grep ${p} | awk '{print $1}'`
  end
endif

echo "monitor kill" >> $SCRATCH/pid.log

sleep 5
touch $SCRATCH/monitor_frealign.log
sleep 5
touch $SCRATCH/monitor_frealign.log
if ( $restart == "T" ) then
  grep "monitor refine" $SCRATCH/pid.log >> /dev/null
  if ( ! $status ) then
    echo "Restarting refinement..." >> frealign.log
    set cycle = `grep Cycle frealign.log | tail -1 | awk '{print $2}' | awk -F: '{print $1}'`
    mv frealign.log frealign_crash.log
    ${bin_dir}/frealign_run_refine restart cycle
  endif
endif
