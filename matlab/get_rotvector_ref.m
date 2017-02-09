function d_ref = get_rotvector_ref( nbins, compute_mode )
% Figures out how much probability should be at
%  each point in an N x N x N grid
%  that tiles -pi to pi cube of 'rotation-vector' space.
% 
%
% Integrates (numerically):
%
%  (1/8pi^2) * sin( v/2 )^2/(v/2)^2
%
% over each voxel.
%
% INPUTS:
%   nbins = no. bins, centered at -pi, ..., 0,... +pi (default is 11)
%   compute_mode = 0 (use from disk if available, otherwise use integral3)
%                  1 (use integral3) -- slow
%                  2 (use voxel blocksum) -- fast
%   
%
% (C) R. Das, Stanford University, 2017

if ~exist( 'nbins', 'var' ) nbins = 11; end;
if ~exist( 'compute_mode' ) compute_mode = 0; end;

d_ref = [];

% name of cached d_ref
[pathstr,~,~] = fileparts(which( mfilename ));
workspace_file= [pathstr, '/rotvector_ref_',num2str(nbins),'bins.mat' ];
if ( compute_mode == 0 )
    if exist( workspace_file, 'file' )
        fprintf( 'Reading d_ref from %s\n', workspace_file );
        ws = load( workspace_file, 'd_ref' );
        d_ref = ws.d_ref;
        return;
    end
    compute_mode = 1;
end

if ( compute_mode == 1)
    d_ref = run_integral3( nbins );
    fprintf( 'Saving d_ref to %s\n', workspace_file );
    save( workspace_file, 'd_ref' );
    return;
else
    assert( compute_mode == 2 )
    d_ref = run_block_sum( nbins );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_ref = run_integral3( nbins )
% uh, I actually didn't finish this.
v = @(x,y,z) sqrt(x.*x + y.*y + z.*z);
% integral3 is not compatible with discontinuity at v = pi
% stepf = @(x,y,z) v(x,y,z)<pi
stepf = @(x,y,z) 1-sigmf( v(x,y,z), [100,pi] );
fun   = @(x,y,z) (1/8/pi^2) * sinc( v(x,y,z)/2/pi ).^2 .* stepf(x,y,z);

binwidth = 2*pi/(nbins-1);
xbounds = [-pi-(binwidth/2): binwidth : pi+(binwidth/2)];
assert( length( xbounds ) == nbins+1 );
for i = 1:nbins
    for j = 1:nbins
        for k = 1:nbins
            fprintf( 'Computing at voxel (%f,%f,%f)\n',...
                     (xbounds(i)+xbounds(i+1))/2,...
                     (xbounds(j)+xbounds(j+1))/2,...
                     (xbounds(k)+xbounds(k+1))/2 );
                     
            d_ref(i,j,k) = ...
                integral3( fun, xbounds(i), xbounds(i+1),...
                xbounds(j), xbounds(j+1),...
                xbounds(k),xbounds(k+1) );
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_ref = run_block_sum( nbins )

binwidth = 2*pi/(nbins-1);
nsub = 20;
binwidth_fine = binwidth/nsub;

Vx = [-pi-(binwidth/2)+(binwidth_fine/2) : binwidth_fine : pi+(binwidth/2)-(binwidth_fine/2)];
Vy = Vx;
Vz = Vx;
assert( length( Vx ) == nbins * nsub );

[VX,VY,VZ] = meshgrid( Vx, Vy, Vz );
V = sqrt( VX.^2 + VY.^2 + VZ.^2 );

% in MATLAB, sinc(x) = sin(x*pi)/(x*pi)
d = (1/8/pi^2) * sinc( V/2/pi ).^2. .* (V < pi) * (binwidth_fine^3); % probability *density* integrated over voxel

% check sum
fprintf( 'Check that integral %f is close to 1.0\n', sum(d(:)) );

% sum over nsub x nsub x nsub blocks
d_ref = squeeze( sum( sum( sum( reshape(d,[nsub nbins nsub nbins nsub nbins]), 1 ), 3 ), 5 ) );


