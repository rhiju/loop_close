function plot_xyz_hist( scorefilename, maxdev, binsize, axis_angle );
% plot_xyz_hist( scorefilename, maxdev );
%
% Inputs:
%  scorefilename = output of Rosetta with 
%                   iter, score, round, x, y, z, alpha, beta, gamma
%  maxdev        = max deviation (Angstroms) from origin. (40 A default)
%  binsize       = bin size for plotting (2 A default)
%  axis_angle    = 0: use Euler_angle
%                  1: convert Euler to axis angle (i.e. "rotation vector")
%                  2: scorefile already in axis angle (not Euler)
%
if ~exist( 'scorefilename','var');  scorefilename = 'gag.SCORES.txt'; end
if ~exist( 'maxdev','var') maxdev = 40; end;
if ~exist( 'binsize','var') binsize = 2; end;
if ~exist( 'axis_angle','var') axis_angle = 0; end;

d = load( scorefilename );
x = d(:,4);
y = d(:,5);
z = d(:,6);
if axis_angle < 2
 phi   = d(:,7); 
 theta = d(:,8); 
 psi   = d(:,9);
 rotation_vector = get_rotation_vector_from_euler( phi, theta, psi ) * pi / 180.0;
else
 rotation_vector = d(:,7:9);
 mag = sqrt( sum( rotation_vector.*rotation_vector, 2 ) );
 ev = rotation_vector ./ repmat( mag, [1,3] );
 ev(:,4) = mag; % must be in degrees for SpinCalc
 ea = SpinCalc( 'EVtoEA313',ev,1.0e-6,0 );
 phi   = ea(:,1);
 theta = ea(:,2);
 psi   = ea(:,3);
end

maxval = max( [abs(x),abs(y),abs(z)] );
if maxval > maxdev
    fprintf( 'max x,y,z %f is greater than maxdev %f, so increase maxdev\n', maxval,maxdev);
    return;
end

figure(1)
clf;
plot3( x, y, z , '.' );
gp = (theta <= 90.0 );
plot3( x(gp), y(gp), z(gp) , '.' ); hold on
gp = (theta >= 90.0 );
plot3( x(gp), y(gp), z(gp) , '.' ); hold on
plot3( 0,0,0,'ko' );
xlabel( 'x' ); ylabel( 'y' ); zlabel( 'z' );
set(gcf, 'PaperPositionMode','auto','color','white');
hold off
axis vis3d
axis equal
view(-150,12);

% histogram
bins = [-maxdev:binsize:maxdev]; n = numel( bins );
xr = interp1( bins, 1:n, x, 'nearest' );
yr = interp1( bins, 1:n, y, 'nearest' );
zr = interp1( bins, 1:n, z, 'nearest' );
W = accumarray( [xr yr zr], 1, [n n n ]);
[w,i] = max(W(:));
[i,j,k]=ind2sub(size(W),i);

% Convert to effective molarity
C = w/length(x)/(binsize^3)/(6.022e23)/1e-27;
fprintf( 'Maximum density %d (C=%f M) at %f,%f,%f\n',w,C,bins(i),bins(j),bins(k));

figure(2)
clf;
gp = find( xr==i & yr==j & zr==k ); 
plot_rotations( phi, theta, psi, rotation_vector, gp, axis_angle );
title( 'Rotations at maximally populated x,y,z' );


% go back to figure(1) and draw a cube showing this bin
figure(1)
verts = ([0 0 0;0 1 0;1 1 0;1 0 0;0 0 1;0 1 1;1 1 1;1 0 1]-0.5) * binsize + repmat(bins([i,j,k]),[8 1]);
face = [1 2 3 4;5 6 7 8;3 4 8 7;1 2 6 5;2 3 7 6;4 1 8 5];
h = patch('Faces',face,'Vertices',verts,'FaceColor','b','EdgeColor','k','FaceAlpha',0,'linewidth',1.5);

histogram_file = strrep( scorefilename, 'SCORES.txt','HISTOGRAM.bin.gz' );
if exist( histogram_file, 'file' ) 
    T = read_tensor( histogram_file );
    figure(1);
    plot_6d_hist_projxyz( T );
    figure(2);
    plot_6d_hist_rotvector( T, bins( [i,j,k] )  );
end


return;
figure(3)
gp = find( z < 0 ); 
plot_rotations( phi, theta, psi, rotation_vector, gp, axis_angle );
title( 'Rotations for z < 0' );

figure(4)
% special position ('direct stacking')
xb = interp1( bins, 1:n, 2, 'nearest' );
yb = interp1( bins, 1:n, 0, 'nearest' );
zb = interp1( bins, 1:n, 4, 'nearest' );
gp = find( xr == xb & yr == yb & zr == zb );
plot_rotations( phi, theta, psi, rotation_vector, gp, axis_angle );
title( 'euler angles at (2,0,4) [direct stack location]' );
C = length(gp)/length(x)/(binsize^3)/(6.022e23)/1e-27;
fprintf( 'Density %d (C=%f M) at %f,%f,%f\n',w,C,2,0,4 );

