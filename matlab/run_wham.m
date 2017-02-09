function T = run_wham( h, bias_strength, bias_xyz, json );
% T = run_wham( h, bias_strength, bias_xyz );
%
% See analyze_histograms.m for how to get set up to run this command.
%
% INPUTS
%  h = 7D tensor of histograms (Nv, Nv, Nv, Nr, Nr, Nr, Nexpt),
%       where Nv is no. bins in rotation vector,
%             Nr is no. bins in translation,
%             Nexpt is no. runs with different bias potentials
% bias_strength = weight on harmonic bias, (1, Nexpt)
% bias_xyz      = target vector for harmonic bias, (3, Nexpt)
% json          = .json struct from read_tensor of histogram file
%
% OUTPUT
%  T = 6D tensor with final probability distribution (and json as field)
%
% (C) Rhiju Das, Stanford University

assert( length( size( h ) ) == 7 );
h_size = size( h );
Nexpt = size( h, 7 );

% since bias is over xyz, faster to just
% project over rotations to get xyz:
hxyz = squeeze( sum( sum( sum( h, 4 ), 5 ), 6 ) );
c = reshape( hxyz, [prod(h_size(1:3)), Nexpt] );
xbins = [json.minval(1) : json.binwidth(1) : json.maxval(1) ];
ybins = [json.minval(2) : json.binwidth(2) : json.maxval(2) ];
zbins = [json.minval(3) : json.binwidth(3) : json.maxval(3) ];
[X,Y,Z] = ndgrid( xbins, ybins, zbins );
for i = 1:Nexpt
    wxyz(:,:,:,i) = exp( -bias_strength( i ) * ...
        ( (X - bias_xyz(1,i)).^2 + ...
          (Y - bias_xyz(2,i)).^2 + ...
          (Z - bias_xyz(3,i)).^2 ) );
end
w = reshape( wxyz, [prod(h_size(1:3)), Nexpt] );

% probability distribution
f = zeros( size( w, 1 ), 1 );
% normalization constants
z = ones( Nexpt, 1 );

for iter = 1:100
    f = squeeze( sum( c, 2 ) )./( w * z );
    f = f/sum(f);
    z_prev = z;
    z = squeeze( sum( c, 1 ) )'./( w' * f );
    if ( norm( z_prev - z )/norm(z) < 1.0e-12 ) break; end;
end
fprintf( 'Ran WHAM on (x,y,z) for %d iterations\n', iter );
fxyz = reshape( f, h_size(1:3) );

% where is max point?
fprintf( 'Figure out max points...\n', iter );
[w,idx] = max(fxyz(:));
[i,j,k]=ind2sub(size(fxyz),idx);
maxpt = [xbins(i),ybins(j),zbins(k)];

% now apply the normalization coefficients z to the full 6D tensor
fprintf( 'Calculating full 6D tensor\n', iter );
c = reshape( h, [prod(h_size(1:6)), Nexpt ] );
w_full = repmat( wxyz, [1 1 1 h_size(4:6) 1] );
w = reshape( w_full, [prod(h_size(1:6)), Nexpt ] );
f = squeeze( sum( c, 2 ) )./( w * z );
f = reshape( f, h_size(1:6) );

T.tensor = f/sum(f(:));
T.json   = json;

% plot the results
figure(1)
clf; plot3( maxpt(1), maxpt(2), maxpt(3), 'kx' );
plot_6d_hist_projxyz( T, 5 );
figure(2)
clf; plot_6d_hist_rotvector( T, maxpt);
figure(3)
visualize_6D_potential( T );
title( '-log(stats) (projection at z=0, vz=0)' )




