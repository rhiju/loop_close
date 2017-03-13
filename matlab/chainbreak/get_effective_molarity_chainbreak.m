function C_eff = get_effective_molarity_chainbreak( x, h, tag )
% C_eff = get_effective_molarity_chainbreak( x, h, tag )
%
%  x = bin centers of histogram
%  h = histogram
%  tag = name of scorefile (will show on plot title)
%
%
xyzfile = [fileparts(which(mfilename)),'/xyz.txt']; % this better have O3'-P-P5' coordinates
xyz  = load( xyzfile );
N = size( xyz, 1 );
% re-center -- note that original code did not have this recentering (bug!).
xyz = xyz - repmat( mean( xyz, 1 ), [N 1]); 
Ixyz = eig( xyz' * xyz ); % moments of inertia

loglog( x, h, 'o' ); hold on; 

% compare to expected.
molar = (1e27/6.022e23);
pred_curve1 = (1.0/molar) * N^(3/2) * x.^5 * pi /8 / sqrt(Ixyz(1)+Ixyz(2)) / sqrt(Ixyz(1)+Ixyz(3)) / sqrt(Ixyz(2)+Ixyz(3));

%x_cutoff = 3;
%x_cutoff = 5;
% find where we have stats, and pick 4th bin.
gp = find( h > 0 ); x_cutoff = x( gp( 4 ) );
C_eff = sum( h( find( x < x_cutoff) ) )/ sum( pred_curve1( find( x < x_cutoff ) ) );
plot( x, C_eff * pred_curve1, 'k' ); hold off
legend( 'Measured', 'Predicted for O3''-P-O5'' at C_{eff}');
xlabel( 'RMSD (Angstroms)');
ylabel( 'Probability distribution (dP/dRMSD)');
out_string = sprintf( 'C(eff) %60s: %9.6f M ', tag, C_eff );
title( out_string, 'interpreter','none','fontsize',7 );
fprintf( [ out_string, '\n' ] );