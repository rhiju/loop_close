function plot_6d_hist_projxyz( T, N_contours )
% T = 6D tensor output by RNA_FragmentMonteCarlo,
%      as read into MATLAB with read_tensor.
% N_contours = # contours to show for 2, 0.2, 0.02, ... [default 2]
%
if ~exist( 'N_contours','var' ) N_contours = 2; end;
assert( length( size( T.tensor ) ) == 6 ); % assume in 6D format, already permuted.
h_trans = squeeze( sum(sum(sum(T.tensor,4),5),6) );
h_trans = h_trans/ sum( h_trans(:) );
xbins = [T.json.minval(1) : T.json.binwidth(1) : T.json.maxval(1) ];
ybins = [T.json.minval(2) : T.json.binwidth(2) : T.json.maxval(2) ];
zbins = [T.json.minval(3) : T.json.binwidth(3) : T.json.maxval(3) ];
[X,Y,Z] = ndgrid( xbins, ybins, zbins );
C = h_trans/prod(T.json.binwidth(1:3))/(6.022e23)/1e-27;

% draw it -- note the permute is to ensure x and y aren't MATLAB-flipped.
contours = [ 2.0, 0.2, 2e-2, 2e-3, 2e-4, 2e-5, 2e-6 ]; 
alpha    = [ 0.5, 0.3, 0.1, 0.05, 0.02, 0.02, 0.02]; 
colors = {'black','blue','cyan','green','yellow',[1 0.5 0],'red'};
for i = 1:min( N_contours, length( contours ) )
    p = patch( isosurface(X,Y,Z,permute(C,[2,1,3]),contours(i)) );
    p.FaceColor = colors{i}; p.EdgeColor = 'none'; p.FaceAlpha = alpha(i);
    hold on
end
axis( [min(xbins) max(xbins) min(ybins) max(ybins) min(zbins) max(zbins)] );
plot3( 0, 0, 0, 'ko' );
camlight; lighting phong
axis vis3d
view(-150,12);
xlabel( 'x' ); ylabel( 'y' ); zlabel( 'z' );
set(gcf, 'PaperPositionMode','auto','color','white');
title( 'Effective molarity contours, projected to x,y,z (2 M, 0.2 M, ...)' )

% double-check 
% max(C(:))
