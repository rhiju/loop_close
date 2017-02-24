#!/bin/bash
#PBS -N _Users_rhiju_Dropbox_projects_RNA_loop_close_runs_swm_P5b_connect_SWM_1
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /Users/rhiju/Dropbox/projects/RNA/loop_close/runs/swm/P5b_connect

/Users/rhiju/src/rosetta//main/source/bin/stepwise -s P5b_connect_START1_1gid_RNAA.pdb -native P5b_connect_NATIVE_1gid_RNAA.pdb -extra_min_res A:150 A:153 -fasta P5b_connect.fasta -nstruct 20 -save_times -score:weights stepwise/rna/rna_res_level_energy4.wts -cycles 2000 -new_move_selector -submotif_frequency 0.2 -atr_rep_screen_for_docking false -superimpose_over_all -motif_mode -submotif_frequency 0.0 -out:file:silent SWM/1/swm_rebuild.out -cycles 200 > 1.log 2> 1.err 
