loop_dirnames = dir( 'loop_*' );
pwdname = pwd;
for q = 3:length( loop_dirnames)
    dirname = loop_dirnames(q).name;
    fprintf( 'Going into directory: %s\n', dirname );
    cd( dirname );
    %analyze_histograms;
    T = read_tensor( 'potential.bin.gz' );
    write_tensor( T, 'potential.txt.gz' );
    cd( pwdname );
end
