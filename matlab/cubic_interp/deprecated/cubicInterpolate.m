function val = cubicInterpolate( F, x, minval, binwidth, periodic )
idx_real = 1 + (x - minval)/binwidth;
idx = floor( idx_real );
N = size( F, 1 );
if periodic
    p = F( mod( [idx-1:idx+2]-1, N ) + 1 );
else % not periodic
    if ( idx < 2 )
        idx = 1;
        F0 = F(3) - 3*F(2) + 3*F(1);
        p = [ F0 F(1:3)' ];
    elseif ( idx >= N - 1 )
        idx = N - 1;
        FNplus1 = 3*F(N)- 3*F(N-1)+F(N-2);
        p = [ F(N-2:N)' FNplus1 ];
    else
        p = F( [idx-1: idx+2] );
    end
end

d = idx_real - idx;

val = basicInterpolate( p, d );
      