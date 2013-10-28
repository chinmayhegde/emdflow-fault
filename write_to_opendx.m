function write_to_opendx(density,fname)

[nx,ny,nz] = size(density);

fid = fopen(fname,'w');

fprintf(fid,'object 1 class gridpositions counts %d %d %d\n',nx,ny,nz);
fprintf(fid,'origin %f %f %f\n',0.0,0.0,0.0);
fprintf(fid,'delta %d %d %d\n',1.0,0.0,0.0);
fprintf(fid,'delta %d %d %d\n',0.0,1.0,0.0);
fprintf(fid,'delta %d %d %d\n',0.0,0.0,1.0);
fprintf(fid,'object 2 class gridconnections counts %d %d %d\n',nx,ny,nz);
fprintf(fid,'object 3 class array type double rank 0 times %d\n',nx*ny*nz);
densityvec = density(:); densityvec = flipud(densityvec);
mm = 0;
while mm < numel(density)            
            fprintf(fid,'%f %f %f\n',densityvec(mm+1),densityvec(mm+2),densityvec(mm+3));
            mm = mm+3;
end
fprintf(fid,'attribute "dep" string "positions"\n');
fprintf(fid,'object "regular positions regular connections" class field\n');
fprintf(fid,'component "positions" value 1\n');
fprintf(fid,'component "connections" value 2\n');
fprintf(fid,'component "data" value 3');

fclose(fid);