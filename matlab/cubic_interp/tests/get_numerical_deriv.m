function [val, deriv_analytic,deriv_numeric] = get_numerical_deriv( F, x, minval, binwidth, boundary );

[val,deriv_analytic] = nCubicInterpolate( F, x, minval, binwidth, boundary);

Ndim = length( size( F ) );
if ( size(F,2) == 1 ) Ndim = 1; end;
delta = 1.0e-8;

for n = 1 : Ndim
    x_delta = x;  x_delta( n ) = x( n ) + delta;
    val_delta = nCubicInterpolate( F, x_delta, minval, binwidth, boundary);
    deriv_numeric(n) =  (val_delta - val ) / delta;
end

fprintf( 'Compare ');
for n = 1:Ndim; fprintf(' %f', deriv_analytic(n)); end;
fprintf(' (analytic) to ');
for n = 1:Ndim; fprintf(' %f', deriv_numeric(n)); end;
fprintf(' (numeric) with boundary %s\n', boundary);