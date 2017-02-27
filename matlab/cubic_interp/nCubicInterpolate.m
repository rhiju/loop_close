function [val,deriv] = nCubicInterpolate( F, x, minval, binwidth, boundary, deriv)
% val = nCubicInterpolate( F, x, minval, binwidth, boundary)
%
F_size = size( F );
if ( F_size(1) == 1 );
    assert( F_size( 2 ) > 1 );
    assert( length( x ) == 1 );
    % MATLAB thing -- squeeze() in  prior recursions
    % disallows formation of a 'matrix' with single row.
    F = F'; 
end
do_basic = ( size( F, 2 ) == 1 );
if ~exist( 'deriv','var' ) deriv = []; end;
N = size( F, 1 );
if ischar( boundary ) 
    boundary_str = boundary;
    boundary = {};
    for i = 1:N; boundary{i} = boundary_str; end;
end

idx_real = 1 + (x(1) - minval(1))/binwidth(1);
idx = floor( idx_real );

[idx, Fslice ] = get_Fslice( idx, F, boundary{1} );
d = idx_real - idx;

if do_basic
    [val, dv_dx] = basicInterpolate( Fslice, d );
    deriv = [ dv_dx/binwidth(1) ];
else
    % here's the recursion: 
    [p(1),dp_dy(1)] = nCubicInterpolate( squeeze(Fslice(1,:)), x(2:end), minval(2:end), binwidth(2:end), boundary(2:end) );
    [p(2),dp_dy(2)] = nCubicInterpolate( squeeze(Fslice(2,:)), x(2:end), minval(2:end), binwidth(2:end), boundary(2:end) );
    [p(3),dp_dy(3)] = nCubicInterpolate( squeeze(Fslice(3,:)), x(2:end), minval(2:end), binwidth(2:end), boundary(2:end) );
    [p(4),dp_dy(4)] = nCubicInterpolate( squeeze(Fslice(4,:)), x(2:end), minval(2:end), binwidth(2:end), boundary(2:end) );
    [val, dv_dx, dv_dp] = basicInterpolate( p, d );
    dv_dy = dv_dp * dp_dy';
    deriv = [ dv_dx/binwidth(1), dv_dy ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx, Fslice ] = get_Fslice( idx, F, boundary );
N = size( F, 1 );
switch boundary
    case 'periodic'
        idx_array = idx + [-1, 0, 1, 2];
        idx_array = mod( idx_array - 1, N ) + 1; % wrap
        % 4 x N2 x N3 x ...
        Fslice = F( idx_array, : );
    case 'flat'
        idxs = [idx-1: idx+2];
        gp = find( idxs >= 1 & idxs <= N );
        Fslice( gp, : ) = F( idxs( gp ), : );    
        bp = find( idxs < 1 );
        Fslice( bp, : ) = repmat( F( 1, : ), [length(bp),1]);
        bp = find( idxs > N );
        Fslice( bp, : ) = repmat( F( N, : ), [length(bp),1] );
    case 'linear'
        idxs = [idx-1: idx+2];
        gp = find( idxs >= 1 & idxs <= N );
        Fslice( gp, : ) = F( idxs( gp ), : );    
        bp = find( idxs < 1 );
        for m = bp;
            Fslice( m,  : ) =  F( 1, : ) + (idxs(m) - 1) * ( F( 2, :) - F( 1, : ) );
        end
        bp = find( idxs > N );
        for m = bp
            Fslice( m, : ) = F( N, : ) + (idxs(m) - N) * ( F( N, :) - F( N-1, : ) );
        end
    case 'cubic' 
        % This is the specification in Keys, 1981:
        % not that this is a bit wasteful since we only really need
        % the part of the tensor whose second indices are near x(2).
        if ( idx < 2 )
            idx = 1;
            F0 = F(3,:) - 3*F(2,:) + 3*F(1,:);
            Fslice = [ F0; F(1:3,:)];
        elseif ( idx >= N - 1 )
            idx = N - 1;
            FNplus1 = 3*F(N,:)- 3*F(N-1,:)+F(N-2,:);
            Fslice = [ F(N-2:N,:); FNplus1 ];
        else
            Fslice = F( [idx-1: idx+2], : );
        end
    otherwise
        val = 0.0;
        fprintf( 'Unrecognized boundary condition: %s\n', boundary{ 1 } );
        return;
end
