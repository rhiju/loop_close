function C_eff = analyze_chainbreak_WHAM( filetag, delx )
% C_eff = analyze_chainbreak_WHAM( filename, delx );
%
% INPUTS:
%  filetag = tag of files with scores  [e.g., "cc_gag_B2_chainbreak" if files are cc_gag_B2_chainbreak*SCORES.txt]
%             these must be text files with lists of chainbreak * chainbreak_weight            
%  delx = histogram spacing
%
if ~exist( 'delx', 'var' ) delx = 0.1; end;

dirs = dir( [filetag,'*SCORES.txt'] );
Nexpt = length(dirs);
rmsd = {}; chainbreak_weights = [];
for i = 1:Nexpt
    scorefile = dirs(i).name;
    col_start = strfind( scorefile, '_chainbreak' ) + length( '_chainbreak' );
    col_end = strfind( scorefile, '.SCORES.txt' );
    cols = strsplit( scorefile( col_start:col_end-1 ), '_' );
    chainbreak_weights(i) = str2double( cols{1} );
    %fprintf( 'Reading in scorefile %s with assumed chainbreak weights %f\n', scorefile, chainbreak_weights(i) );
    data{i} = load( scorefile );
    scores = data{i}(:,3);
    % this is pretty cryptic -- should directly output RMSDs instead of weighted scores!
    rmsd{i} = sqrt( scores/chainbreak_weights(i)/3 ); Nvals = length( rmsd );
end

delx = 0.1;
x = [0:delx:100];
% WHAM
for i = 1:Nexpt
    c(:,i) = hist( rmsd{i}, x ) / length( rmsd ) / delx; % dp/dr.
    if ( chainbreak_weights(i) > 0.005 ) c(:,i) = c(:,i)*0.005/chainbreak_weights(i); end;
    w(:,i) = exp( -3 * chainbreak_weights(i) * x.^2 );
end


% probability distribution
f = zeros( size( w, 1 ), 1 );
% normalization constants
z = ones( Nexpt, 1 );
% from run_wham

for iter = 1:100
    f = squeeze( sum( c, 2 ) )./( w * z );
    f = f/sum(f)/delx;
    z_prev = z;
    z = squeeze( sum( c, 1 ) )'./( w' * f );
    if ( norm( z_prev - z )/norm(z) < 1.0e-12 ) break; end;
end
fprintf( 'Ran WHAM on (x,y,z) for %d iterations\n', iter );

get_effective_molarity_chainbreak( x, f, filetag );
hold on; plot( x, (c./w) * diag(1./z) ); hold off

