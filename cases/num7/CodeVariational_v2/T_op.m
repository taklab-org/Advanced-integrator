function [Ta] = T_op(a)
N = size(a)-1;
Ta = zeros(N+1);
otherdims = repmat({':'},1,length(N)-1);
for i = 1:N(end)-1
    Ta(otherdims{:}, i+1) =  a(otherdims{:}, i+2) - a(otherdims{:}, i);
end
Ta(otherdims{:}, N(end)+1) =  -a(otherdims{:}, N(end));
end
