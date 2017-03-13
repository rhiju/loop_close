function C_eff = analyze_chainbreak( scorefile, chainbreak_weight,  delx)
% C_eff = analyze_chainbreak( scorefile, chainbreak_weight, delx)
%
% If you have run with a bunch of different weights, use 
%  analyze_chainbreak_WHAM
%
% INPUTS:
%  scorefile = text file with scores -- must be chainbreak * chainbreak_weight 
%  chainbreak_weight = weight during rna_denovo/FARFAR (must be > 0.0)
%  delx = histogram spacing
%
if ~exist( 'delx', 'var' ) delx = 0.25; end;
x = [0:delx:100];
h = get_rmsd( x, delx, scorefile, chainbreak_weight );
get_effective_molarity_chainbreak( x,  h, scorefile );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = get_rmsd( x, delx, scorefile, chainbreak_weight );

d = load( scorefile );
% this is pretty cryptic -- should directly output RMSDs instead of weighted scores!
scores = d(:,3);
rmsd = sqrt( scores/chainbreak_weight/3 ); Nvals = length( rmsd );
h_raw = hist( rmsd, x ) / Nvals / delx; % dp/dr.
% correct for chainbreak bonus
weights = exp( 1.0 * scores );
h = (h_raw.*exp( 3 * chainbreak_weight * x.^2 ) ) * Nvals/sum(weights);

