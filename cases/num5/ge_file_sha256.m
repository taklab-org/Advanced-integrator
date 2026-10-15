function hex=ge_file_sha256(file)
fid=fopen(file,'rb');assert(fid>=0,'GE:HashRead','Cannot read %s',file);closer=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while true
 bytes=fread(fid,8*1024*1024,'*uint8');if isempty(bytes),break;end
 md.update(typecast(bytes,'int8'));
end
hex=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
