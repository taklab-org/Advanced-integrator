function test_variational_v2()
% Requires initialized INTLAB and the selected case paths.
% Large-window doubles below are regression diagnostics, not certificates.
assert(contains(which('verify_variationaleq'),'verify_solution'));
root=fileparts(fileparts(mfilename('fullpath')));
for caseName={'num5','num7','pbvp'}
 c=fullfile(root,'cases',caseName{1},'CodeVariational_v2');
 for file={'compute_bounds_cheb.m','compute_Y_cheb.m','rig_variational_cheb.m'}
  assert(strcmp(fileread(fullfile(c,file{1})),fileread(fullfile(root,'cases','num5','CodeVariational_v2',file{1}))));
 end
end
nop=[0,0,0,0;0,0,0,0;0,0,0,-1];
a=zeros(2,2,4);a(1,1,1)=0.2;a(2,1,2)=0.03;
c=zeros(size(a));c(1,1,1)=0.7;c(2,2,2)=-0.1;
actual=N_op_lin_int(nop,intval(a),intval(c));
expected=-3*convapbqcos(intval(a),2,intval(c),1);
assert(all(inf(actual)<=sup(expected) & sup(actual)>=inf(expected),'all'));
for name={'N_op_lin.m','N_op_lin_int.m'}
 s=fileread(which(name{1}));assert(~contains(s,'Nop_coeff(i+1)'));
end
for amplitude=[0.2,0.5]
 ba=zeros(3,2,2);ba(1,1,1)=amplitude;ba(2,2,1)=0.05;ba(3,1,2)=0.02;
 R=v2_window_diagnostic(ba,0.0625,[2,2],[0,1,0.39],[1,1.1],3,16,'low_order',true);
 assert(all(R.r>=R.err),'AI:Regression','Low-order radius underestimates independent window error.');
 assert(R.Z0>=R.tZ0 && R.Zinf>=R.tZinf && R.Zle>=R.tZle && R.ZA>=R.tZA && R.Z1>=R.tZ1);
end
fprintf('VARIATIONAL_V2_REGRESSION_PASSED\n');
end
