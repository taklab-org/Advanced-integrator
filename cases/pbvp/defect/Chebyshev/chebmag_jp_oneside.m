function [f_mag] = chebmag_jp_oneside(a)
% returns the extrema f_min and f_max of a f(t) given as a Chebyshev series
% a: One-sided Chebyshev coefficients a = (a_0,...a_N)
% NOTE :: available for interval inputs, e.g., chebmag(intval(a),I))
N_FFT = 2^14; %power of two and size of the mesh

N = length(a)-1;
% Build the two-sided sequence
a = [flip(a(2:end));a(1);a(2:end)];
n = (-N:N)';

delta = infsup(0,2*pi/N_FFT);
alpha = a.*exp(1i*n*delta);
% We pad alpha with to apply the FFT
alpha = [zeros(N_FFT/2-N,1);alpha;zeros(N_FFT/2-N-1,1)];
f_mag = max(sup(abs(real(verifyfft(ifftshift(alpha))))));
%f_min = min(inf(f_on_mesh));
%f_max = max(sup(f_on_mesh));

end
