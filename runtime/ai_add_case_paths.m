function ai_add_case_paths(c)
addpath(c);
[~,name]=fileparts(c);
if strcmp(name,'pbvp')
 dirs={'defect','defect/Chebyshev','High_Dim_Cosine_Conv', ...
  'CodeVariational','verify_solution','saddle2saddle_connexion'};
else
 dirs={'defect','defect/Chebyshev','High_Dim_Conv','CodeVariational_v2','verify_solution'};
end
for name=dirs
 p=fullfile(c,name{1});if isfolder(p),addpath(p);end
end
end
