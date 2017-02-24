# loop_close
Scripts used to derive 6D potential for RNA loop_close

# Directory structure
+ `runs` has a bunch of different Rosetta runs to derive and then test 6D potentials
+ `matlab` has matlab scripts for visualizing 6D potentials and also running WHAM to compute them from Rosetta runs.
+ `notes` has scans of handwritten notes, as well as a running log of this project.
+ `check_rotations` has some silly scripts to check derivatives of 6D potentials (should move to MATLAB)
More documentation coming soon...

# Getting started
To see how to define a 6D potential from Rosetta, check out:
```
runs/create_potential/wham/README_SETUP
```
This sets up Rosetta `rna_denovo` runs to compute end-to-end statistics with different loop lengths.
There are different runs with biasing potentials that ensure sampling of loops with ends near different `(x,y,z)` translations.

To quickly explain, here's a randomly chosen example from `loop05/`:
```
rna_denovo @flags -sequence gaaaaag -output_jump_res 1 7 -secstruct_legacy HLLLLLH -lores_scorefxn bias.wts -target_xyz 0 0 0 -out:file:silent bias_0_0_0/default.out -output_histogram_file bias_0_0_0/HISTOGRAM.bin.gz
```
This simulates a sequence `gaaaaag` (you could use other ones to steer fragments to, e.g., pyrimidine loops). 

The `-secstruct_legacy` flag specifies that the first and last residues should only draw from fragments that were originally in Watson-Crick helices in the crystallographic database ('H'). 
The other nucleotides must be drawn from fragments that were _not_ in helices ('L'). You can supply 'X' if you don't care. 'N' means the nucleotide must not have formed any pairs (but might be stacked!) in the original database structure.

The `-output_jump_res 1 7` specifies that the jumps (i.e., rotation/translations) should be computed starting at the 5' g and ending at the 3' g.

The `-lores_scorefxn` means use `bias.wts`, which just includes a single term `rna_stub_coord_hack` which implements the biasing potential.
`target_x_y_z` specifies the central translation for the bias. Other flags specify output files, and finally,
the `flags` file is:
```
-temperature 1.0
-close_loops false
-minimize_rna false
-no_filters
-nstruct 10
-farna:rounds 1
-cycles 1000000
-output_jump_o3p_to_o5p
-output_rotation_vector
-output_score_frequency 10
-output_score_file none
-lores_scorefxn none.wts
-save_jump_histogram
-jump_histogram_boxsize 45.000000
-jump_histogram_binwidth 4.500000
```
Here, temperature is set to `1.0` since that's assumed throughout the analysis, no loop closure is on (not relevant to this partiular case),
there is no minimizing (just use fragment assembly monte carlo), there is a single round (this ensures that all score terms are on at full strength, rather than getting titrated up).
`-output_jump_o3p_to_o5p` specifies the atoms at which to compute the rigid_body_transform (centered at O3' at 5' end and O5' at the 3' end) and
the biasing potential. `output_rotation_vector` means use axis-angle representation and not Euler angles.

After doing these runs on a cluster, go into a directory like `loop05` and run `analyze_histograms` in MATLAB; make sure all subdirectories of `matlab/` are in your path.


