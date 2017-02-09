% (C) Rhiju Das, Stanford University
% no bias case
clear dirnames; i = 1;
bias_strength(i) = 0.0;
bias_xyz(:,i) = [0 0 0 ];
dirnames{i} = './';

% look for bias histograms
dirs = dir( 'bias_*');
for j = 1:length( dirs )
    i = i+1;
    dirnames{i} = dirs(j).name;
    bias_strength(i) = 0.01;
    cols = strsplit( dirnames{i}, '_' );
    bias_xyz(:,i) = str2double(cols(2:4));
end

for i = 1:length( dirnames )
    filename = [dirnames{i},'/HISTOGRAM.bin.gz'];
    fprintf( 'Reading in: %s\n', filename );
    T = read_tensor( filename );
    h(:,:,:,:,:,:,i) = T.tensor;
    json_info = T.json;
end

T = run_wham( h, bias_strength, bias_xyz, json_info );

figure(4);
E = derive_potential( T );
title( 'Derived 6D potential (projection at z=0, vz=0)' )

