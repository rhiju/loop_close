nohup /home/rhiju/src/rosetta//main/source/bin/rna_denovo -temperature 1.0 -lores_scorefxn bias.wts -close_loops false -minimize_rna false -no_filters -nstruct 10 -farna:rounds 1 -cycles 100000 -output_jump_o3p_to_o5p -output_rotation_vector -output_score_frequency 100 -target_xyz 0 0 0 -save_jump_histogram -sequence gag -out:file:silent gag.out -output_jump_res 1 3 -secstruct_legacy HLH > /dev/null 2> /dev/null &
