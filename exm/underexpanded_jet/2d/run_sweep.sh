#!/bin/bash
cd /home/qiaoj/testcerisse/cerisse/exm/underexpanded_jet/2d
RANKS=16
STOP=5.0e-3
run_one(){ d=$1; ( cd $d && mpirun --oversubscribe -np $RANKS ./main2d.gnu.TPROF.MPI.ex inputs_ext stop_time=$STOP max_step=400000 cns.nstep_screen_output=500 > run.log 2>&1; echo "DONE $d t=$(grep -oE 'TIME = [0-9.e+-]+' run.log | tail -1)" >> /tmp/sweep_status.txt ); }
echo "SWEEP START $(date)" > /tmp/sweep_status.txt
run_one visc_ext_npr3  & run_one visc_ext_npr4   & wait
run_one visc_ext_npr5  & run_one visc_ext_npr7p5 & wait
run_one visc_ext_npr15 & run_one visc_ext_npr20  & wait
echo "SWEEP ALL DONE $(date)" >> /tmp/sweep_status.txt
