function [Ta] = T_op_abs_int(a)
N = size(a)-1;
Ta = intval(zeros(N+1));
otherdims = repmat({':'},1,length(N)-1);
if N(end) >= 2
    Ta(otherdims{:}, 2:N(end)) = abs(a(otherdims{:}, 3:N(end)+1)) + abs(a(otherdims{:}, 1:N(end)-1));
end
Ta(otherdims{:}, N(end)+1) = abs(a(otherdims{:}, N(end)));
end
