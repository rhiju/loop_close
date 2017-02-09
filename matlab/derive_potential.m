function E = derive_potential( T )
% E = derive_potential( T )
%
%
nrotbins = size( T.tensor, 4 );
assert( size( T.tensor, 5 ) == nrotbins );
assert( size( T.tensor, 6 ) == nrotbins );
assert( T.json.minval(4) == -180.0 );
assert( T.json.minval(5) == -180.0 );
assert( T.json.minval(6) == -180.0 );
assert( T.json.maxval(4) == +180.0 );
assert( T.json.maxval(5) == +180.0 );
assert( T.json.maxval(6) == +180.0 );

% following will have value 1/8pi^2 at (0,0,0) grid point. Its a
% probability density in 3D rotation space.
d_ref = get_rotvector_ref( nrotbins );

% Convert to probability *density* at 1 M -- expand into 3 translation
% dimensions.
fprintf( 'Taking log ratio...\n' );
T_size = size( T.tensor );
d_ref_6d = permute( repmat( d_ref, [1 1 1 T_size(1:3)] ), [ 4 5 6 1 2 3] );
%d_ref_6d = repmat( d_ref, [T_size(1:3) 1 1 1 ] );

F = -log( T.tensor/ prod(T.json.binwidth(1:3)) /(6.022e23) /1e-27 ./ d_ref_6d );
F(isnan(F)) = Inf; % Inf is a sign of missing data.

E = T;
E.tensor = F;
E.json.type = class( F ); % double

% let's pad by one extra bin in each of the rotation vector directions
% this will allow me to ensure continuous interpolation upon traversing
% through sphere of radis pi
E = pad_at_pi( E );

% let's do fills of V>pi points.
E = fill_2pi( E );

% get rid of inf -- use nearest non-inf neighbor.
E = fill_inf_by_knn( E );

visualize_6D_potential( E, 0, 0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = pad_at_pi( E );
fprintf( 'Padding rotation vector dimensions by 1 bin beyond pi\n' );
F = E.tensor;
F = padarray( F, [0 0 0 1 1 1], Inf );
E.tensor = F;

json = E.json;
json.minval( 4:6 ) = json.minval( 4:6 ) - json.binwidth( 4:6 );
json.maxval( 4:6 ) = json.maxval( 4:6 ) + json.binwidth( 4:6 );
json.n_bins = size( E.tensor );
E.json = json;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = fill_2pi( E )
% In rotation-vector space, there is a mapping:
%   v -> v + (2 pi) v/|v|
%
F = E.tensor;
fprintf( 'Running fill based on 2pi symmetry mates...\n' );
[x,y,z,Vx,Vy,Vz] = get_bin_centers( E.json );
[VX,VY,VZ] = ndgrid(Vx,Vy,Vz);
V = sqrt( VX.^2 + VY.^2 + VZ.^2 );
Vmap = [ VX(:), VY(:), VZ(:) ];
Vunit  = diag( 1./V(:) ) * Vmap; % unit vectors
Vunit( isnan( Vunit ) ) = 0.0;
VunitX = diag( V(:)>180.0 ) * Vunit; % unit vectors only for V>pi points
Vmap = Vmap - 360.0*VunitX;

% do the mapping
[X,Y,Z,VX,VY,VZ] = ndgrid( x,y,z,Vx,Vy,Vz );
F_size = size( F );
Nxyz = prod( F_size(1:3) );
Vmapx = repmat(Vmap(:,1),[1 Nxyz])';
Vmapy = repmat(Vmap(:,2),[1 Nxyz])';
Vmapz = repmat(Vmap(:,3),[1 Nxyz])';
Fmap = F(:);
bp = find( isinf( F ) );
Fmap(bp) = interpn( X,Y,Z,VX,VY,VZ, F,...
    X(bp),Y(bp),Z(bp),Vmapx(bp),Vmapy(bp),Vmapz(bp) );

% reshape back out to 6D
F = reshape(Fmap,F_size);

% try to fill?
%fprintf( 'Running imfill...\n' );
%F = imfill(F); %,isinf(F));

E.tensor = F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = fill_inf_by_knn( E )
fprintf( 'Filling points with no stats based on nearest neighbor...\n' );
[x,y,z,Vx,Vy,Vz] = get_bin_centers( E.json );
[X,Y,Z,VX,VY,VZ] = ndgrid( x,y,z,Vx,Vy,Vz );

F = E.tensor;
gp = find( ~isinf( F ) ); % good points

% set an (X,Y,Z) within which we actually have counts
R = sqrt( X.^2+Y.^2+Z.^2 );
Rmax = max( R( gp ) );
bp = find( isinf( F ) & (R < Rmax) ); % bad points
fprintf( 'Filling %d no-statistics points, within radius of %f Angstroms, based on nearest neighbor.\n',...
    length(bp), Rmax );
K = 10; % there are 12 possible direct neighbors -- let's sample from all.
idx = knnsearch( ...
    [X(gp),Y(gp),Z(gp),VX(gp),VY(gp),VZ(gp)],...
    [X(bp),Y(bp),Z(bp),VX(bp),VY(bp),VZ(bp)], 'K', K );
             
% Use max since this was a low statistics point -- energy should be
% 'worst-case' estimate.
F( bp ) = max(F( gp(idx) ),[],2);

Fmax = max( F( ~isinf(F) ) );
F( isinf(F) ) = Fmax;
fprintf( 'Filling rest of no-statistics points to be max energy %f.\n', Fmax );


E.tensor = F;

function [x,y,z,Vx,Vy,Vz] = get_bin_centers( json );
x  = [json.minval(1): json.binwidth(1) : json.maxval(1)];
y  = [json.minval(2): json.binwidth(2) : json.maxval(2)];
z  = [json.minval(3): json.binwidth(3) : json.maxval(3)];
Vx = [json.minval(4): json.binwidth(4) : json.maxval(4)];
Vy = [json.minval(5): json.binwidth(5) : json.maxval(5)];
Vz = [json.minval(6): json.binwidth(6) : json.maxval(6)];
