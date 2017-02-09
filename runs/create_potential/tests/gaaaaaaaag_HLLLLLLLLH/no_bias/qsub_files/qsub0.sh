#!/bin/bash
#PBS -N _biox3_scratch_users_rhiju_projects_RNA_nn_rule_chain_closure_loop_close_runs_create_potential_tests_gaaaaaaaag_HLLLLLLLLH_no_bias_.
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /biox3/scratch/users/rhiju/projects/RNA/nn_rule_chain_closure/loop_close/runs/create_potential/tests/gaaaaaaaag_HLLLLLLLLH/no_bias

/home/rhiju/src/rosetta//main/source/bin/rna_denovo -temperature 1.0 -lores_scorefxn just_chainbreak.wts -close_loops false -minimize_rna false -no_filters -nstruct 10 -farna:rounds 1 -cycles 100000 -output_jump_o3p_to_o5p -output_rotation_vector -output_score_frequency 100 -save_jump_histogram -sequence gaaaaaaaag -out:file:silent gaaaaaaaag.out -output_jump_res 1 10 -secstruct_legacy HLLLLLLLLH > /dev/null 2> /dev/null 
