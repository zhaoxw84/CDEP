function [newlabel] = LabelNorm(label)
%�������㷨����õ��������ı�ǩ��Ϣ���±���Ϊ1��K����ʽ

[N,C] = size(label);
if N < C
    label = label';
    N = C;
end
value = unique(label);
K = length(value);
newlabel = label;

for i=1:K
    index =  label == value(i);
    newlabel(index) = i;
end
end

