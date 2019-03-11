function C_eff = make_6D_plots( T, maxpt, dirname )

if ~exist( 'dirname','var') dirname = './'; end;

% plot the results
figure(1)
clf; plot3( maxpt(1), maxpt(2), maxpt(3), 'kx' );
plot_6d_hist_projxyz( T, 7 );
expfig( [dirname,'/projxyz_contours.png'] );

figure(2)
clf; C_eff = plot_6d_hist_rotvector( T, maxpt, 4);
expfig( [dirname, '/rotvec_contours.png'] );

figure(3)
visualize_6D_potential( T );
title( [strrep(pwd(),'_','\_'),'\newline-log(stats) (projection at z=0, vz=0)'] )
expfig( [dirname, '/proj4d_logstats.pdf'] );




