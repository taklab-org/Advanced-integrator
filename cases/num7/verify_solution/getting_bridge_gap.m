function z = getting_bridge_gap(error,frame,nu)
% A signed endpoint/initial mismatch is transformed BEFORE taking abs.
m=frame.m;n=prod(m);e=intval(error);ef=e(1:m(1),1:m(2));
z=intval(zeros(n+1,1));z(1:n)=abs(frame.Vi*frame.D*ef(:));
e(1:m(1),1:m(2))=0;z(end)=wnorm(e,intval(nu));
end
