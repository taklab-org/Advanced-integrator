function cm = getting_cm(ba_sup,m,nu)
m1 = m(1); m2 = m(2);
%
psi = intval(zeros(m));
%
for k1=1:m1
  for k2=1:m2
    psi(k1,k2) = psi_k(ba_sup,k1-1,k2-1,m,nu);
  end
end
%
cm = 3*wnorm(psi,nu);
