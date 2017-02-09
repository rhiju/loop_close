function T = read_tensor( filename )
% read_tensor( filename )
%
%  Read tensor from file with name like "my_tensor.bin" or
%  "my_tensor.bin.gz"
%
%  and look for accompanying .json file with n_bins information.
%
% (C) Rhiju Das, Stanford University

T.tensor = [];
if (~strcmp(filename(end-3:end),'.bin') & ...
    ~strcmp(filename(end-6:end),'.bin.gz') )
    fprintf( 'Filename should have .bin or .bin.gz as suffix\n');
    return;
end
% need to gunzip binary files
if strcmp( filename(end-2:end), '.gz' )
    filename = filename(1:end-3); % will gunzip in next block
end
if ~exist( filename, 'file' ); 
    filename_gz =  [filename, '.gz' ];
    if exist( filename_gz, 'file' )
        fprintf( 'Gunzipping: %s\n', filename_gz );
        gunzip( filename_gz ); 
    end
end

json_file = strrep( strrep( filename, '.bin', '.json' ), '.gz', '' );
if ~exist( json_file, 'file' ); 
    json_gz =  [json_file, '.gz' ];
    if ~exist( json_gz, 'file' )
        fprintf( 'Could not find JSON file with name like %s or %s.\n',json_file,json_gz );
        return;
    end
    gunzip( [json_file, '.gz'] ); 
end
json = loadjson( json_file );
assert( isfield( json, 'n_bins' ) );
assert( isfield( json, 'type' ) );

fid = fopen( filename, 'r' );
[X, n_data] = fread( fid, json.type );
if ( prod( json.n_bins ) ~= n_data ) 
    fprintf( 'Mismatch between n_data %d and expected value based on product of n_bins %d\n', n_data, prod( json.n_bins) );
    return;
end
fclose( fid );

% the ordering of MathNTensor output in Rosetta requires reshaping and
% permuting to get into 'reasonable' order.
Ndim = length( json.n_bins );
T.tensor = reshape( X, json.n_bins(Ndim:-1:1) );
T.tensor = permute( T.tensor, [Ndim:-1:1] );
T.json = json;
