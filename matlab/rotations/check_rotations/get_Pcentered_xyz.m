function xyz = get_Pcentered_xyz( xyz_file );

if ~exist( 'xyz_file' ) xyz_file = 'xyz_O3p_P_O5p.txt'; end
xyz = load( xyz_file );

% zero center
xyz = xyz - repmat( xyz(2,:), [size(xyz,1),1] );

% want P --> O5' to go along z -- i.e. set up coordinate system used by
%  rna_denovo to define rotations and translations:
z = xyz(3,:) - xyz(2,:);
z = z/norm(z);
x = xyz(2,:) - xyz(1,:);
y = cross(z,x);
y = y/norm(y);
x = cross(y,z);
R = [x;y;z]
xyz = xyz * R';

