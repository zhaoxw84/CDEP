function Q=modularity_metric(label,adj)

nedges=numedges(adj); % total number of edges
Kvalue = unique(label);
K = length(Kvalue);
modules = cell(K,1);
for i=1:K
    modules{i} = find(label == Kvalue(i));
end

Q = 0;
for m=1:length(modules)
    mm=sum(sum(adj(modules{m},modules{m})))/2;
    e_mm=mm/nedges;
    a_m=sum(sum(adj(modules{m},:)))/(2*nedges);
    Q = Q + (e_mm - a_m^2);
  
end

% Kvalue = unique(label);
% K = length(Kvalue);
% modules = cell(K,1);
% for i=1:K
%     modules{i} = find(label == Kvalue(i));
% end
% 
% nedges=numedges(adj); % compute the total number of edges
% 
% Q = 0;
% for m=1:length(modules)  
% 
%   e_mm=numedges(adj(modules{m},modules{m}))/nedges;
%   a_m=sum(sum(adj(modules{m},:)))/nedges - e_mm; % counting e_mm only once
%   Q = Q + (e_mm - a_m^2);
% end