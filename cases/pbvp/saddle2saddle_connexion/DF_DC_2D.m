function DF = DF_DC_2D(a,params,N,L)
% Computes an approximation of DF using central differences
% a does not contain the zero mode
% input size(a) = (N1*N2-1, 1): vector
% output size(F) = (N1*N2-1, 1): vector form
% params = [a0,sigma,epsilon];

m = N(1)*N(2)-1; %total size
DF = zeros(m);
h = 1e-6; % finite difference parameter
I = eye(m);

for j=1:m
    a_plus = a + h*I(:,j);
    a_minus = a - h*I(:,j);
    DF(:,j) = (F_DC_2D(a_plus,params,N,L)-F_DC_2D(a_minus,params,N,L))/(2*h);
end

end
