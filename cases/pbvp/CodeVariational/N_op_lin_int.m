function [Nout] = N_op_lin_int(Nop_coeff,a_bar,c)
Nout = intval(zeros(size(c)));
% N = size(c)-1;
% dim=length(N);
% for j=1:dim
%     s{j}=1:N(j)+1;
% end
for i = 1:(size(Nop_coeff,2)-1)
%     for  j = 0:(size(Nop_coeff,1)-1)
        if Nop_coeff(i+1) ~= 0
                apc = convapbqcos(a_bar,i-1,c,1);
                Nout = Nout + i*Nop_coeff(i+1)*apc;
        end
%     end
end
