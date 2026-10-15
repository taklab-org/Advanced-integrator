function manifest=ge_source_manifest(root)
f=dir(fullfile(root,'**','*.m'));names={};hashes={};
for k=1:numel(f)
 path=fullfile(f(k).folder,f(k).name);name=path(numel(root)+2:end);
 parts=strsplit(name,filesep);
 if any(ismember(parts,{'runs','archive','saddle2saddle_connexion'})),continue;end
 names{end+1}=name;hashes{end+1}=ge_file_sha256(path);
end
[names,ix]=sort(names);manifest=struct('names',{names},'sha256',{hashes(ix)});
end
