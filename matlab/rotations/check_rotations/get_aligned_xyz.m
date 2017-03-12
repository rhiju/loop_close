function [xyz, Ixyz] = get_aligned_xyz( xyz_file );

if ~exist( 'xyz_file' ) xyz_file = 'xyz_O3p_P_O5p.txt'; end
xyz = load( xyz_file );

% zero center
xyz = xyz - repmat( mean(xyz), [size(xyz,1),1] );
[V,D] = eig( xyz' * xyz ); % moments of inertia
xyz = xyz*V;
Ixyz = diag(D);