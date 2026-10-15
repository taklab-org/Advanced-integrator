function check_package()
%CHECK_PACKAGE Source-only release audit and inventory; MATLAB only.
root=fileparts(fileparts(mfilename('fullpath')));
files=dir(fullfile(root,'**','*'));records=struct('path',{},'sha256',{});
for k=1:numel(files)
 if files(k).isdir,continue;end
 p=fullfile(files(k).folder,files(k).name);rel=p(numel(root)+2:end);
 if any(ismember(strsplit(rel,filesep),{'.git','outputs'})),continue;end
 [~,~,ext]=fileparts(p);
 assert(~ismember(lower(ext),{'.mat','.log','.dll','.so','.dylib','.zip'}) ...
  &&~startsWith(lower(ext),'.mex'),'Binary/calculation data in release: %s',rel);
 if strcmpi(ext,'.m')
  assert(isempty(regexp(fileread(p),'(?:[A-Z]:\\|/Users/|/home/)\w','once')), ...
   'Host-specific path in %s',rel);
 end
 rel=strrep(rel,filesep,'/');
 if strcmp(rel,'docs/release-files.json'),continue;end
 fid=fopen(p,'rb');assert(fid>0);b=fread(fid,Inf,'*uint8');fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(b,'int8'));
 hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 records(end+1)=struct('path',rel,'sha256',hash); %#ok<AGROW>
end
s=fileread(fullfile(root,'cases','num7','proof_GE_2DOK_core.m'));
assert(contains(s,'tau = 8;')&&contains(s,'K_F_list = 128;') ...
 &&contains(s,'m_a=[12,12];h_auto=0.1'));
assert(~contains(s,'8.625')&&~contains(s,'K_F_list = 138'));
assert(isfile(fullfile(root,'cases','num5','initial_num5.m')) ...
 &&isfile(fullfile(root,'cases','num7','initial_num7.m')));
[~,ix]=sort({records.path});records=records(ix);
fid=fopen(fullfile(root,'docs','release-files.json'),'wb');assert(fid>0);
closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(records,'PrettyPrint',true));
fprintf('PASS: %d source/document files; no binary computation data or host paths.\n',numel(records));
end
