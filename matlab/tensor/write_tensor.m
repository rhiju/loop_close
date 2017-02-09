function write_tensor( T, filename, do_gzip )
% write_tensor( T, filename )
%
%  Create tensor with name like "my_tensor.bin.gz"
%  and create accompanying .json file with n_bins information.
%
% INPUTS
%   T = struct with .tensor and .json fields
%   filename = output filename
%   do_gzip  = gzip the .bin file [default: 1]
%
% (C) Rhiju Das, Stanford University 2017
if ~exist( 'do_gzip', 'var' ) do_gzip = 1; end;

assert( strcmp( T.json.type, class( T.tensor ) ) );

if (~strcmp(filename(end-3:end),'.bin') & ...
    ~strcmp(filename(end-6:end),'.bin.gz') )
    fprintf( 'Filename should have .bin or .bin.gz as suffix\n');
    return;
end
if strcmp( filename(end-2:end), '.gz' )
    filename = filename(1:end-3);
end
json_file = strrep( strrep( filename, '.bin', '.json' ), '.gz', '' );

% the ordering of MathNTensor output in Rosetta requires reshaping and
% permuting to get into 'reasonable' order.
Ndim = length( size( T.tensor ) );
F = permute( T.tensor, [Ndim:-1:1] );

fid = fopen( filename, 'w');
count = fwrite( fid, F(:), class(F));
fprintf( 'Put %d vals into: %s \n', count, filename );
fclose( fid );
if do_gzip
    gzip( filename );
    delete( filename )
    fprintf( 'Gzipped into %s.gz\n', filename );
end

json = savejson( '', T.json, json_file );
fprintf( 'Created companion JSON file: %s\n', json_file );
