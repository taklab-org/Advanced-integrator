function [Nout] = N_op_lin(Nop_coeff,a_bar,c)
assert_supported_nop_coeff(Nop_coeff);
Nout = zeros(size(c));
for i = 1:(size(Nop_coeff,2)-1)
    for j = 0:(size(Nop_coeff,1)-1)
        if Nop_coeff(j+1,i+1) ~= 0
            apc = convapbqcos(a_bar,i-1,c,1);
            Nout = Nout + i*Nop_coeff(j+1,i+1)*apc;
        end
    end
end
end
