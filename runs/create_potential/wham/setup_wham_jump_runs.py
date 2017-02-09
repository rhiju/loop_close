#!/usr/bin/python

import argparse
import os

parser = argparse.ArgumentParser(description='Run rna_helix to build a helix')
parser.add_argument('-length', required=True, help='Length of helix', nargs='?', type=int )
parser.add_argument('-points',default=[-10,10],help='points in x (y,z) for biasing',nargs='+',type=int )
parser.add_argument('-boxsize',default=40.0,help='boxsize for +/-x,y,z',nargs='?',type=float )
parser.add_argument('-binwidth',default=4.0,help='binwidth for x,y,z',nargs='?',type=float )
parser.add_argument('-skip_no_bias',help='skip the no-bias run',action='store_true' )
parser.add_argument('-append',help='append to README instead of overwriting',action='store_true' )
args = parser.parse_args()

assert( args.length >= 2 )
sequence = 'g'
for i in range( args.length-2 ): sequence += 'a'
sequence += 'g'

secstruct = 'H'
for i in range( args.length-2 ): secstruct += 'L'
secstruct += 'H'

outdir = 'loop_%02d/' % args.length
if not os.path.exists( outdir ): os.makedirs( outdir )

out_flags = open( outdir + 'flags', 'w' )
out_flags.write('-temperature 1.0\n')
out_flags.write('-close_loops false\n')
out_flags.write('-minimize_rna false\n')
out_flags.write('-no_filters\n')
out_flags.write('-nstruct 10\n')
out_flags.write('-farna:rounds 1\n')
out_flags.write('-cycles 1000000\n')
out_flags.write('-output_jump_o3p_to_o5p\n')
out_flags.write('-output_rotation_vector\n')
out_flags.write('-output_score_frequency 10\n')
out_flags.write('-output_score_file none\n')
out_flags.write('-lores_scorefxn none.wts\n')
out_flags.write('-save_jump_histogram\n')
out_flags.write('-jump_histogram_boxsize %f\n' % args.boxsize)
out_flags.write('-jump_histogram_binwidth %f\n' % args.binwidth)
out_flags.close()

out_wts = open( outdir + 'none.wts', 'w' )
out_wts.close()

out_wts = open( outdir + 'bias.wts', 'w' )
out_wts.write( 'rna_stub_coord_hack 0.01\n' )
out_wts.close()


if args.append:
    out = open( outdir+'README', 'a')
else:
    out = open( outdir+'README', 'w')

# baseline command
cmd_baseline = 'rna_denovo @flags -sequence %s -output_jump_res 1 %d -secstruct_legacy %s' % ( sequence,args.length,secstruct )
count = 0

# no-bias command
if not args.skip_no_bias:
    cmd = cmd_baseline + ' -output_histogram_file HISTOGRAM.bin.gz'
    out.write( cmd+'\n' )
    count += 1
# now put in biasing potentials
bias_pts = args.points
for x in bias_pts:
    for y in bias_pts:
        for z in bias_pts:
            bias_dir = 'bias_%d_%d_%d' % (x,y,z)
            if not os.path.exists( outdir+bias_dir ): os.makedirs( outdir + bias_dir )
            cmd = cmd_baseline + ' -lores_scorefxn bias.wts'
            cmd += ' -target_xyz %d %d %d' % (x,y,z)
            cmd += ' -out:file:silent %s/default.out' % (bias_dir)
            cmd += ' -output_histogram_file %s/HISTOGRAM.bin.gz' % (bias_dir)
            out.write( cmd+'\n' )
            count += 1
out.close()

print 'All %d command lines are in: %sREADME' % (count,outdir)


