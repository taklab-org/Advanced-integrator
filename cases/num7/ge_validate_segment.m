function ge_validate_segment(s)
assert(isfinite(s.schur_margin) && s.schur_margin>0,'GE:SegmentSchur','Invalid Schur margin');
names={'ba','banorm','F1_fin','F1_inf','F2_fin','F2_inf','int_semi','int_semi_q', ...
 'psi_vec','c_out_pm','VinfJ','VinfqJ','g22','end_error','a_end'};
for k=1:numel(names)
 x=intval(s.(names{k}));assert(all(isfinite(sup(abs(x))),'all'),'GE:SegmentFinite','Invalid segment field %s',names{k});
end
for n={'F1_fin','F1_inf','F2_fin','F2_inf','banorm'}
 assert(all(inf(intval(s.(n{1})))>=0,'all'),'GE:SegmentBound','Negative bound');
end
end
