function C_eff = plot_6d_hist_rotvector( T, xyzval, N_contours )
% T = 6D tensor output by RNA_FragmentMonteCarlo,
%      as read into MATLAB with read_tensor.
% xyzval = x,y,z position from which to draw rotation vector distributions
%
if ~exist( 'N_contours','var' ) N_contours = 2; end;
h = T.tensor;
h_size = size( h );
assert( length( h_size ) == 6 ); % better be 6D tensor
h = h/sum(h(:));
h_trans = squeeze( sum(reshape( h, [h_size(1:3) prod(h_size(4:6))] ), 2) );

% figure out where xyzval is:
binsizes = T.json.binwidth;
xbins = [T.json.minval(1) : T.json.binwidth(1) : T.json.maxval(1) ];
ybins = [T.json.minval(2) : T.json.binwidth(2) : T.json.maxval(2) ];
zbins = [T.json.minval(3) : T.json.binwidth(3) : T.json.maxval(3) ];
xb = interp1( xbins, 1:h_size(1), xyzval(1), 'nearest' );
yb = interp1( ybins, 1:h_size(2), xyzval(2), 'nearest' );
zb = interp1( zbins, 1:h_size(3), xyzval(3), 'nearest' );

% just rotational stats that correspond to that xyzval translation
hr = squeeze( h( xb, yb, zb, :, :, : ) );
C = sum(hr(:))/prod(binsizes(1:3))/(6.022e23)/1e-27;
fprintf( 'Effective molarity (averaged over rotations) at (%f,%f,%f) is: %f M\n', xyzval(1:3), C );

vxbins = [T.json.minval(4) : T.json.binwidth(4) : T.json.maxval(4) ];
vybins = [T.json.minval(5) : T.json.binwidth(5) : T.json.maxval(5) ];
vzbins = [T.json.minval(6) : T.json.binwidth(6) : T.json.maxval(6) ];
[VX,VY,VZ] = ndgrid( vxbins, vybins, vzbins );

% if we convert to radians, we can compare to uniform distribution that
% would form at 1 M. In rotation coordinates, would expect 
%
%   1/(8pi^2) * sinc( V/2 )
%
% Actually, should be more careful that binsize is large -- do integral of
% above over bin, zeroing at V > pi (!), rather than just evaluate at center.
%
Cr = hr/prod(binsizes(1:6))/(6.022e23)/1e-27/(pi/180)^3;
V = sqrt( VX.^2 + VY.^2 + VZ.^2 )*(pi/180.0);
% Cr_no_sinc = Cr/(1/(8*pi^2));
uniform_rot_density = (1/(8*pi^2)) * (sinc(V/2/pi)).^2;
Cr = Cr./uniform_rot_density;
C_eff = max(Cr(:));
if length( xyzval ) == 3 
    fprintf( 'Effective molarity (at most favored rotation): %f M\n', max(Cr(:)) );
else
    assert( length( xyzval ) == 6 )
    vxb = interp1( vxbins, 1:h_size(4), xyzval(4), 'nearest' );
    vyb = interp1( vybins, 1:h_size(5), xyzval(5), 'nearest' );
    vzb = interp1( vzbins, 1:h_size(6), xyzval(6), 'nearest' );
    fprintf( 'Effective molarity (at rotation %f, %f, %f): %f M\n', xyzval(4:6),Cr(vxb,vyb,vzb));  
    C_eff = Cr(vxb,vyb,vzb);
end

% draw it -- note the permute is to ensure x and y aren't MATLAB-flipped.
hold on
contours = [ 20.0, 0.2, 2e-2, 2e-3, 2e-4, 2e-5 ]; 
alpha    = [ 0.5, 0.3, 0.1, 0.05, 0.02, 0.02, 0.02]; 
colors = {'black','blue','cyan','green','yellow',[1 0.5 0],'red'};
for i = 1:min( N_contours, length( contours ) )
    p = patch( isosurface(VX,VY,VZ,permute(Cr,[2,1,3]),contours(i)) );
    p.FaceColor = colors{i}; p.EdgeColor = 'none'; p.FaceAlpha = alpha(i);
    hold on
end

axis( [min(vxbins) max(vxbins) min(vybins) max(vybins) min(vzbins) max(vzbins)] );
plot3( 0, 0, 0, 'ko' );
camlight; lighting phong
axis vis3d
xlabel( 'vx' ); ylabel( 'vy' ); zlabel( 'vz' );
