function ai_write_json(file,value)
fid=fopen(file,'w');assert(fid>0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
end
