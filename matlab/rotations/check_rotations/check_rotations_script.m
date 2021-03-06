N = 100000;
M = get_uniform_random_matrices( N );
EA = SpinCalc( 'DCMtoEA313', M, 1.0e-6, 0 );
figure(1)
% convert to axis_angle and plot -- uniform? yes.
plot_rotations( EA );

% axis angles
figure(2);
% convert to axis_angle and plot -- smooth? yes
v = plot_rotations( EA, 1 );

figure(3)
clf;
binsize = 0.1;
vbins = [0:binsize:pi];
vmag = (pi/180) * sqrt( sum(v.*v,2) );
h = hist( vmag, vbins );
% density in axis-angle space:
h = h / sum( h ) / binsize ./ ( 4 * pi * vbins.^2);
plot( vbins, h, 'x' );
hold on
plot( vbins, (sin(vbins/2).^2)./(vbins/2).^2 * 1/(8*pi^2) );
legend( 'numerical', '1/(8\pi^2)sinc^2(v/2)');
xlabel( 'v (angle of rotation in radians)' );
ylabel( 'Density in rotation vector space' );
hold off



