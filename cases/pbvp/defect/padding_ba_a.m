function [ba_padded,a_padded,max_rows] = padding_ba_a(ba,a)
max_rows = max(size(ba,1), size(a,1));
ba_padded = intval(zeros(max_rows, size(ba,2)));
a_padded  = zeros(max_rows, size(a,2));
ba_padded(1:size(ba,1), :) = ba;
a_padded(1:size(a,1), :) = a;
