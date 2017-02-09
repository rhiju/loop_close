function expfig(  filename )
% expfig( filename )
if exist( 'export_fig', 'file' ) 
    fprintf( 'Creating: %s\n',filename );
    export_fig( filename ); 
end;
