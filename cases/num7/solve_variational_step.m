function [bUn, bUn_endpt, r1, verification] = solve_variational_step(ba, h, m, params, L)
% Same trajectory polynomial and strict variational tests; enlarge only
% its exact-zero Chebyshev representation for the transient Phi proof.
original_rows=size(ba,1);minimum_rows=24;
if original_rows<minimum_rows
 original=ba;sz=size(ba);sz(1)=minimum_rows;
 if isa(ba,'intval'),ba=intval(zeros(sz));else,ba=zeros(sz,'like',original);end
 ba(1:original_rows,:,:)=original;
 assert(isequal(ba(1:original_rows,:,:),original));
 appended=intval(ba(original_rows+1:end,:,:));assert(all(inf(appended)==0 & sup(appended)==0,'all'));
end
[bUn,r1,verification]=verify_variationaleq(ba,h,m,params,L);
bUn_endpt=getting_Un_endpt(bUn,r1);
assert(all(isfinite(sup(abs(bUn_endpt))),'all'),'GE:VariationalEndpoint','Nonfinite endpoint enclosure.');
verification.representation_original_rows=original_rows;
verification.representation_rows=size(ba,1);
end
