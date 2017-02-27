setup_cubic_test_functions;

% check 1D
x = 0.1211;
[val,deriv_analytic] = nCubicInterpolate( F_1D, x, minval, binwidth, 'periodic' );

delta = 1.0e-7;
val_delta = nCubicInterpolate( F_1D, x + delta, minval, binwidth, 'periodic' );
deriv_numeric = ( val_delta - val )/delta;

fprintf( 'Compare %f (analytic) to %f (numeric)\n',deriv_analytic,deriv_numeric);


% check 2D
x = [0.1211,0.5555];
[val,deriv_analytic] = nCubicInterpolate( F, x, [minval,minval], [binwidth,binwidth], 'periodic' );
Ndim = 2;
for n = 1 : Ndim
    x_delta = x;
    x_delta( n ) = x( n ) + delta;
    [val_delta,deriv_analytic] = nCubicInterpolate( F, x_delta, [minval,minval], [binwidth,binwidth], 'periodic' );
    deriv_numeric(n) =  (val_delta - val ) / delta;
end
fprintf( 'Compare %f %f (analytic) to %f %f (numeric)\n',deriv_analytic,deriv_numeric);

% check 3D
