function test = test_cubic_radius(Y,Vq,anorm)
% Common cubic invariance/contraction test; Y is a verified defect bound.
Y=intval(Y);
test=struct('passed',false,'Y',Y,'Pr',intval(Inf),'Zv',intval(Inf),'radius',NaN);
y=sup(Y);v=sup(intval(Vq));a=sup(intval(anorm));candidate=y;
for it=1:60
    C=v*(3*a+candidate);disc=1-4*C*y;
    if ~isfinite(disc) || disc<=0,return;end
    next=2*y/(1+sqrt(disc)); % Avoid cancellation in 1-sqrt(disc).
    if ~isfinite(next),return;end
    if abs(next-candidate)<=1e-12*max(next,1e-30)
        candidate=next;break;
    end
    candidate=next;
end
test.radius=candidate*(1+1e-4)+1e-18;ri=intval(test.radius);
test.Pr=Y+intval(Vq)*(3*intval(anorm)+ri)*ri^2-ri;
test.Zv=intval(Vq)*3*(2*intval(anorm)+ri)*ri;
test.passed=sup(test.Pr)<0 && sup(test.Zv)<1;
end
