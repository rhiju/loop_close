function [val,deriv] = nCubicInterpolate( F, x, minval, binwidth, boundary )
% val = nCubicInterpolate( F, x, minval, binwidth, boundary)
%
Ndim = length( x );
if ischar( boundary ); boundary = repmat({boundary},[1,Ndim]); end
if length( minval ) < Ndim; minval = repmat( minval, [1,Ndim] ); end;
if length( binwidth ) < Ndim; binwidth = repmat( binwidth, [1,Ndim] ); end;
    
idx_real = 1 + (x(1) - minval(1))/binwidth(1);
idx = floor( idx_real );

[idx, Fslice ] = get_Fslice( idx, F, boundary{1} );
d = idx_real - idx;

if ( Ndim == 1 )
    [val, dv_dx] = basicInterpolate( Fslice, d );
    deriv = [ dv_dx/binwidth(1) ];
else
    % here's the recursion: 
    F_size = size( F );
    F_size_onefewerdim = F_size(2:end);
    if Ndim == 2;  F_size_onefewerdim = [ F_size(2), 1 ]; end; % reshape requires 2 indices. Silly MATLAB.
    for i = 1:4
        % need to look at i-th slice of Fslice -- the reshape
        % lets us pretend that its now a tensor with 1 fewer dimension.
        Fslice_i = reshape(Fslice(i,:),F_size_onefewerdim);
        [p(i),dp_dy(:,i)] = nCubicInterpolate(Fslice_i, x(2:end), minval(2:end), binwidth(2:end), boundary(2:end) );
    end
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
