function v = getting_integrated_exponential(lambda,h)
% Enclose integral_0^h exp(lambda*s) ds, including lambda=0.
lambda=intval(lambda);hi=intval(h);assert(inf(hi)>=0);
if isinf(sup(hi))
    assert(all(sup(lambda(:))<0));v=-1./lambda;return;
end
v=intval(zeros(size(lambda)));nonzero=inf(lambda)>0 | sup(lambda)<0;
v(nonzero)=(exp(lambda(nonzero)*hi)-1)./lambda(nonzero);
v(~nonzero)=hi.*exp(lambda(~nonzero)*infsup(0,sup(hi)));
v=abs(v);
end
