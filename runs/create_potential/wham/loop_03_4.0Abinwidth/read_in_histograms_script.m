% no bias case
clear dirnames
i = 1;
bias_strength(i) = 0.0;
bias_xyz(:,i) = [0 0 0 ];
dirnames{i} = './';

for x = [-5 5];
    for y = [-5 5];
        for z = [-5 5];
            i = i+1;
            bias_strength(i) = 0.01;
            bias_xyz(:,i) = [x y z];
            dirnames{i} = sprintf('bias_%d_%d_%d',x,y,z);
        end
    end
end

for x = [-10 10];
    for y = [-10 10];
        for z = [-10 10];
            i = i+1;
            bias_strength(i) = 0.01;
            bias_xyz(:,i) = [x y z];
            dirnames{i} = sprintf('bias_%d_%d_%d',x,y,z);
        end
    end
end

Nexpt = length( dirnames );
Nvxyz = 11;
Nxyz = 21;
if ~exist( 'h', 'var' )
    for i = 1:Nexpt
        filename = [dirnames{i},'/HISTOGRAM.bin.gz'];
        fprintf( 'Reading in: %s\n', filename );
        T = read_tensor( filename );
        h(:,:,:,:,:,:,i) = T.tensor;
        json_info = T.json;
    end
end

f = run_wham( h, bias_strength, bias_xyz, json_info );

