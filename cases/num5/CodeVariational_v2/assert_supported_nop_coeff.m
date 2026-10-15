function assert_supported_nop_coeff(Nop_coeff)
% The current variational code supports a single active differential-order row.
if size(Nop_coeff, 1) <= 1
    return
end

lower_rows = Nop_coeff(1:end-1, 2:end);
if isa(lower_rows, 'intval')
    lower_rows = mid(lower_rows);
end

if any(lower_rows(:) ~= 0)
    error(['Nop_coeff with nonzero entries below the highest differential-order row ' ...
        'is not supported rigorously by this implementation.'])
end
end
