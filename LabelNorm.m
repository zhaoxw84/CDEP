function [newlabel] = LabelNorm(label)
%将各个算法计算得到的类结果的标签信息重新编码为1：K的形式

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

