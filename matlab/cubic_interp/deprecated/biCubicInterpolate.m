function val = biCubicInterpolate( F, x, minval, binwidth, periodic )
N = size( F, 1 );
idx_real = 1 + (x(1) - minval(1))/binwidth(1);
idx = floor( idx_real );
idx_array = idx + [-1, 0, 1, 2];
idx_array = mod( idx_array - 1, N ) + 1; % wrap
p = [ cubicInterpolate( F(idx_array(1),:)', x(2), minval(2), binwidth(2), periodic(2) ),
      cubicInterpolate( F(idx_array(2),:)', x(2), minval(2), binwidth(2), periodic(2) ),
      cubicInterpolate( F(idx_array(3),:)', x(2), minval(2), binwidth(2), periodic(2) ),
      cubicInterpolate( F(idx_array(4),:)', x(2), minval(2), binwidth(2), periodic(2) ) ];
d = idx_real - idx;
val = basicInterpolate( p, d );

