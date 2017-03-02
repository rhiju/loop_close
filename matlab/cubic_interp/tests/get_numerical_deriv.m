function [val, deriv_analytic,deriv_numeric] = get_numerical_deriv( F, x, minval, binwidth, boundary );

interpolate = @nCubicInterpolate;
%interpolate = @nCubicInterpolate_WORKS;

[val,deriv_analytic] = interpolate( F, x, minval, binwidth, boundary);

Ndim = length( size( F ) );
if ( size(F,2) == 1 ) Ndim = 1; end;
delta = 1.0e-8;

for n = 1 : Ndim
    x_delta = x;  x_delta( n ) = x( n ) + delta;
    val_delta = interpolate( F, x_delta, minval, binwidth, boundary);
    deriv_numeric(n) =  (val_delta - val ) / delta;
end

fprintf( '\nValue at ' );
for n = 1:Ndim; fprintf(' %f', x(n)); end;
fprintf( ': %f\n', val);
fprintf( 'Compare ');
for n = 1:Ndim; fprintf(' %f', deriv_analytic(n)); end;
fprintf(' (analytic) to \n        ');
for n = 1:Ndim; fprintf(' %f', deriv_numeric(n)); end;
fprintf(' (numeric) with boundary %s\n', boundary);