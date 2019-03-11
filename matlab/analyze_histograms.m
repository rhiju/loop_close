function [T,C_eff] = analyze_histograms( dirname, test_point, look_for_bias_runs )
% T = analyze_histograms( dirname )
% 
%  dirname = directory with HISTOGRAM.bin.gz file (default = current
%                     directory)
%  test_point = (x,y,z,vx,vy,vz) at which to make 3D projection plots and
%                  estimate C_eff (effective concentration). 
%               Give 3 vals (x,y,z) if you only know that. (C_eff will go
%                to max over vx,vy,vz)
%               Give empty set or nothing, and C_eff will be calculatd at
%                 max.
%
% (C) Rhiju Das, Stanford University 2017.
if ~exist( 'dirname', 'var' ) dirname = './'; end;
if ~exist( 'test_point', 'var' ) test_point = []; end;
if ~exist( 'look_for_bias_runs','var') look_for_bias_runs = 0; end;
    
% no bias case
clear dirnames; i = 1;
dirnames{i} = dirname;
bias_strength = 0.0; bias_xyz = [0,0,0]';

% look for bias histograms
if look_for_bias_runs
    dirs = dir( [dirname,'/bias_*']);
    for j = 1:length( dirs )
        i = i+1;
        dirnames{i} = [dirname,'/',dirs(j).name];
        bias_strength(i) = 0.01;
        cols = strsplit( dirs(j).name, '_' );
        bias_xyz(:,i) = str2double(cols(2:4));
        if length( cols ) > 4 & strcmp( cols{5}, 'weight' )
            bias_strength(i) = str2double( cols{6} );
        end
    end
end

for i = 1:length( dirnames )
    filename = [dirnames{i},'/HISTOGRAM.bin.gz'];
    fprintf( 'Reading in: %s\n', filename );
    T = read_tensor( filename );
    h(:,:,:,:,:,:,i) = T.tensor;
    if isfield( T.json, 'xyz_bias_weight' ); bias_strength(i) = T.json.xyz_bias_weight; end
    if isfield( T.json, 'target_xyz' ); bias_xyz(:,i) = T.json.target_xyz; end
    json_info = T.json;
end

C_eff = 0.0;
if length( dirnames ) == 1 & isempty( test_point )
    test_point = [0,0,0,0,0,0];
end
[T, C_eff] = run_wham( h, bias_strength, bias_xyz, json_info, test_point );
if isempty( test_point )
    figure(4);
    E = derive_potential( T );
    write_tensor( E, 'potential.bin' );
    title( 'Derived 6D potential (projection at z=0, vz=0)' )
    expfig( 'proj4d_potential.pdf' );
end



