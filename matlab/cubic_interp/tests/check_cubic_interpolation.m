setup_test_functions;

subplot(4,2,1);
imagesc( ctrs_fine, ctrs_fine, Ff', [0 6]);
title( 'target function' );
ylim = get(gca,'ylim'); xlim = get(gca,'xlim');

subplot(4,2,2);
imagesc( ctrs, ctrs, F', [0 6]);
title( 'sampled points' );
set(gca,'ylim',ylim,'xlim',xlim);

subplot(4,1,2);
clear F_1Di*
for i = 1:length( ctrs_fine )
    F_1Dip( i ) = nCubicInterpolate( F_1D, ctrs_fine(i), minval, binwidth, 'periodic' );
    F_1Dif( i ) = nCubicInterpolate( F_1D, ctrs_fine(i), minval, binwidth, 'flat' );
    F_1Dil( i ) = nCubicInterpolate( F_1D, ctrs_fine(i), minval, binwidth, 'linear' );
    F_1Dic( i ) = nCubicInterpolate( F_1D, ctrs_fine(i), minval, binwidth, 'cubic' );
end
plot( ctrs_fine, Ff_1D,'k' );
hold on
plot( ctrs, F_1D, 'o' );
plot( ctrs_fine, F_1Dip );
plot( ctrs_fine, F_1Dif );
plot( ctrs_fine, F_1Dil );
plot( ctrs_fine, F_1Dic );
title( 'transect through x = 0' );
legend( 'target','sampled','interp (periodic)','interp (extrap flat)',...
        'interp (extrap-linear)','interp (extrap-cubic)' );
    
subplot(4,1,3);
for i = 1:length( ctrs_fine )
    for j = 1:length( ctrs_fine )
        Fip( i, j ) = nCubicInterpolate( F, ctrs_fine([i,j]), [minval,minval], [binwidth,binwidth], 'periodic' );
        Fif( i, j ) = nCubicInterpolate( F, ctrs_fine([i,j]), [minval,minval], [binwidth,binwidth], 'flat' );
        Fil( i, j ) = nCubicInterpolate( F, ctrs_fine([i,j]), [minval,minval], [binwidth,binwidth], 'linear' );
        Fic( i, j )  = nCubicInterpolate( F, ctrs_fine([i,j]), [minval,minval], [binwidth,binwidth], 'cubic' );
    end
end

subplot(4,2,5);
imagesc( ctrs_fine, ctrs_fine, Fip', [0 6]);
title( 'interp (periodic)' );

subplot(4,2,6);
imagesc( ctrs_fine, ctrs_fine, Fif', [0 6]);
title( 'interp (extrap flat)' );

subplot(4,2,7);
imagesc( ctrs_fine, ctrs_fine, Fil', [0 6]);
title( 'interp (extrap linear)' );

subplot(4,2,8);
imagesc( ctrs_fine, ctrs_fine, Fic', [0 6]);
title( 'interp (extrap cubic)' );
