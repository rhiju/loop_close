function [val,deriv] = nCubicInterpolate( F, x, minval, binwidth, boundary )
% [val, deriv] = nCubicInterpolate( F, x, minval, binwidth, boundary)
%
Ndim = length( x );
if ischar( boundary ); boundary = repmat({boundary},[1,Ndim]); end
if length( minval ) < Ndim; minval = repmat( minval, [1,Ndim] ); end;
if length( binwidth ) < Ndim; binwidth = repmat( binwidth, [1,Ndim] ); end;

% make 4x4x4...4 slice of input tensor --
% need to fill out boundaries if we're beyond
% where convolution is defined
[ Fslice, d ] = get_Fslice( F, x, minval, binwidth, boundary);
[val, deriv ] = tensorInterpolate( Fslice, d );
deriv = deriv ./ binwidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val, deriv ] = tensorInterpolate( Fslice, d );
Ndim = length( d );
if ( Ndim == 1 )
    [val, dv_dx] = basicInterpolate( Fslice, d );
    deriv = [ dv_dx ];
else
    % here's the recursion: 
    F_size = size( Fslice );
    F_size_onefewerdim = F_size(2:end);
    if Ndim == 2;  F_size_onefewerdim = [ F_size(2), 1 ]; end; % reshape requires 2 indices. Silly MATLAB.
    for i = 1:4
        % need to look at i-th slice of Fslice -- the reshape
        % lets us pretend that its now a tensor with 1 fewer dimension.
        Fslice_i = reshape(Fslice(i,:),F_size_onefewerdim);
        [p(i),dp_dy(:,i)] = tensorInterpolate(Fslice_i, d(2:end) );
    end
    [val, dv_dx, dv_dp] = basicInterpolate( p, d(1) );
    dv_dy = dv_dp * dp_dy'; % chain rule
    deriv = [ dv_dx, dv_dy ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Fslice, d ] = get_Fslice( F, x, minval, binwidth, boundary );
% 4 x 4 x 4 ... 4 patch of values from F, with neighbors and
% nearest-neighbors of point x
Ndim = length( x );

slice_size = 4 * ones(1,Ndim);
if ( Ndim == 1 ) slice_size = [4 1]; end;
Fslice = reshape( zeros( 4^Ndim, 1), slice_size );

idx_real = 1 + (x - minval)./binwidth;
idx = floor( idx_real );
d = idx_real - idx;
for k = 1 : length (Fslice(:) ) % 4^Ndim
    [p{1:Ndim}] = ind2sub( size( Fslice ), k ); % this 'decodes' the indices
    idx_check = idx + ( cell2mat( p ) - 2 );  % convert p of 1,2,3,4 to -1,0,1,2
    Fslice( k ) = get_val( F, idx_check, boundary, 1 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = get_val( F, idx, boundary, which_dim );
%
% F   = tensor with Ndim dimensions
% idx = (x,y,..) with Ndim integer elements.
%
% Job of this function is to figure out value of F at
%  index idx. Usually this is just F( idx ). 
%
% But when extrapolating,
%  this function needs to figure out value as some kind of linear
%  combination of values where idx is in range. 
%
% Strategy: handle each dimension of the idx one by one through
%  a recursion. n tells us where we are in the recursion
%
if ( which_dim > length( idx ) ) % recursion is complete
    assert( all( idx >= 1 ) );
    F_size = size( F );
    assert( all( idx <= F_size( 1:length(idx) )  ) );
    for k = 1:length( idx ) idx_cell{k} = idx( k ); end;
    idx_overall = sub2ind( size( F ), idx_cell{:} );
    val = F( idx_overall );
else
    N = size( F, which_dim );
    switch boundary{ which_dim }
        case 'periodic'
            idx( which_dim ) = mod( idx( which_dim ) - 1, N ) + 1; % wrap
            val = get_val( F, idx, boundary, which_dim + 1 );
        case 'flat'
            if idx( which_dim ) <= 1
                idx( which_dim ) = 1;
            elseif idx( which_dim ) >= N
                idx( which_dim ) = N;
            end
            val = get_val( F, idx, boundary, which_dim + 1 );
        case {'linear','cubic'}
            if idx( which_dim ) <= 1
                idx_1 = idx; idx_1( which_dim ) = 1;
                idx_2 = idx; idx_2( which_dim ) = 2;
                val = (2 - idx(which_dim)) * get_val( F, idx_1, boundary, which_dim+1 ) + ...
                    (idx(which_dim) - 1) * get_val( F, idx_2, boundary, which_dim+1 );
            elseif idx( which_dim ) > N
                idx_N       = idx; idx_N( which_dim )  = N;
                idx_Nminus1 = idx; idx_Nminus1( which_dim )= N-1;
                val = (idx(which_dim) - N + 1) * get_val( F, idx_N, boundary, which_dim+1 ) + ...
                    (N - idx(which_dim) )    * get_val( F, idx_Nminus1, boundary, which_dim+1 );
            else
                val = get_val( F, idx, boundary, which_dim+1 );
            end
            % This is the specification in Keys, 1981
            %         if ( idx < 2 )
            %             idx = 1;
            %             F0 = F(3,:) - 3*F(2,:) + 3*F(1,:);
            %             Fslice = [ F0; F(1:3,:)];
            %         elseif ( idx >= N - 1 )
            %             idx = N - 1;
            %             FNplus1 = 3*F(N,:)- 3*F(N-1,:)+F(N-2,:);
            %             Fslice = [ F(N-2:N,:); FNplus1 ];
            %         else
            %             Fslice = F( [idx-1: idx+2], : );
            %         end
        otherwise
            val = 0.0;
            fprintf( 'Unrecognized boundary condition: %s\n', boundary{ which_dim } );
            return;
    end
end
