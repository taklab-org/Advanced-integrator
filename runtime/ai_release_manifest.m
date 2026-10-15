function records=ai_release_manifest(root)
% Bind all distributed MATLAB source, including runtime and problem inputs.
f=dir(fullfile(root,'**','*.m'));records=struct('path',{},'sha256',{});
for k=1:numel(f)
 p=fullfile(f(k).folder,f(k).name);rel=p(numel(root)+2:end);
 if any(strcmp(strsplit(rel,filesep),'outputs')),continue;end
 fid=fopen(p,'rb');assert(fid>0);b=fread(fid,Inf,'*uint8');fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(b,'int8'));
 h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 records(end+1)=struct('path',strrep(rel,filesep,'/'),'sha256',h); %#ok<AGROW>
end
[~,ix]=sort({records.path});records=records(ix);
end
