function [Label,Times] = CDEP(adj)
%Description: convert the vertex folding of degree 1 and 2 in the original network graph into a new network graph, and then automatically select the community center and propagate the label
t0 = clock;
N = size(adj,1); 
%%%%%%%%%%Count the degree of each vertex
degree = zeros(N,2);
degree(:,1) = 1:N;  
degree(:,2) = sum(adj(:,:));
Olddegree= degree;

Label =zeros(N,1);
outlier = find(degree(:,2) == 0); 
Label(outlier) = -1;                

if ~isempty(outlier)
      fprintf('This data has outliers!!!!!');
end
KNN = cell(N,1);     
Quality = ones(N,1); 
Nodelist = num2cell(1:N);

%%%%=============The first stage: compress the nodes with degree 1 and 2
[degreeoneid,~] = find(degree(:,2)==1);
[degreetwoid,~] = find(degree(:,2)==2);
flag = 1;  
while(flag)
    flag =0;
%%%%The node with degree 1 is compressed
    while ~isempty(degreeoneid) 
          currentnode = degreeoneid(1);
          %neighindex = KNN{currentnode}; 
          neighindex = find(adj(currentnode,:)); 
          
          Nodelist{neighindex}  = union(Nodelist{neighindex},Nodelist{currentnode}); 
         % KNN{neighindex}  = setdiff(KNN{neighindex},[currentnode]); 
          Quality(neighindex) = Quality(neighindex) + Quality(currentnode);
          degree(neighindex,2) = degree(neighindex,2) -1;
          if (degree(neighindex,2) == 2)
              degreetwoid = union(degreetwoid,neighindex);
          elseif (degree(neighindex,2) == 1)
             degreeoneid = union(degreeoneid,neighindex);
             degreetwoid = setdiff(degreetwoid,neighindex);
          elseif (degree(neighindex,2) == 0)
             degreeoneid = setdiff(degreeoneid,neighindex);
          end

          Nodelist{currentnode}  = [];
         % KNN{currentnode}  = []; 
          Quality(currentnode) = 0;
          degree(currentnode,2) = 0;
          degreeoneid = setdiff(degreeoneid,currentnode);

          adj(currentnode,neighindex) = 0;
          adj(neighindex,currentnode) = 0;
    end 

%%%The node with degree 2 is compressed
    while ~isempty(degreetwoid)  
          currentnode = degreetwoid(1);
         % neighindex = KNN{currentnode}; 
          neighindex = find(adj(currentnode,:)); 
       %   if (sum(adj(neighindex(1),:) & adj(neighindex(2),:)) >1) && (adj(neighindex(1),neighindex(2))==1)  
          if (adj(neighindex(1),neighindex(2))==1)  
          [~,maxindex] = max(degree(neighindex,2));
          for i = 1:2 
            %  KNN{neighindex(i)}  = setdiff(KNN{neighindex(i)},currentnode);
              if i == maxindex
               Nodelist{neighindex(i)}  = union(Nodelist{neighindex(i)},Nodelist{currentnode}); 
               Quality(neighindex(i)) = Quality(neighindex(i)) + Quality(currentnode);
              end
          end
          if(adj(neighindex(1),neighindex(2)) ~= 0)
              degree(neighindex(1),2) = degree(neighindex(1),2) -1;
              if (degree(neighindex(1),2) == 2)
                degreetwoid = union(degreetwoid,neighindex(1));
              elseif (degree(neighindex(1),2) == 1)
                degreeoneid = union(degreeoneid,neighindex(1));
                degreetwoid = setdiff(degreetwoid,neighindex(1));
              end

               degree(neighindex(2),2) = degree(neighindex(2),2) -1;
              if (degree(neighindex(2),2) ==2)
                degreetwoid = union(degreetwoid,neighindex(2));
              elseif (degree(neighindex(2),2) == 1)
                degreeoneid = union(degreeoneid,neighindex(2));
                degreetwoid = setdiff(degreetwoid,neighindex(2));
              end
          end

          %The weight change of the related edge after compressing the point with degree 2
          adj(neighindex(1),neighindex(2)) = adj(neighindex(1),neighindex(2)) + 0.5 * adj(neighindex(1),currentnode) * adj(neighindex(2),currentnode); %对折叠后的边权重进行修改
          adj(neighindex(2),neighindex(1)) = adj(neighindex(1),neighindex(2));
          adj(currentnode,neighindex(1)) = 0;  adj(neighindex(1),currentnode) = 0;
          adj(neighindex(2),currentnode) = 0;  adj(currentnode,neighindex(2)) = 0;
          
          
          degreetwoid = setdiff(degreetwoid,currentnode);
          Nodelist{currentnode}  = [];
          Quality(currentnode) = 0;
          degree(currentnode,2) = 0;   
          else
              degreetwoid = setdiff(degreetwoid,currentnode);
         end
    end 

    if (~isempty(degreetwoid) || ~isempty(degreeoneid))
        flag = 1;
    end
end
foldweight = adj; 



%=============The second stage: determine the initial community center according to the index of quality and degree

 node_no = []; 
 NN = 0;  
for i=1:N
    if ~isempty(Nodelist{i})
        NN = NN + 1;
        node_no = [node_no,i];
        Plotquality(NN) = Quality(i); 
        Plotdegree(NN) = degree(i,2);
        if degree(i,2) == 0 && Quality(i) ~= 1
           Plotdegree(NN) =  mean(Olddegree(Nodelist{i},2)); 
        end
    end
 end
Plotquality = Plotquality/max(Plotquality);
Plotdegree = Plotdegree/max(Plotdegree);
fprintf('The number of nodes after Compressed: %i \n', NN);
edges = (size(find(foldweight(node_no,node_no) ~= 0),1))/2;
fprintf('The number of edges after Compressed: %i \n', edges);

for i = 1:NN
   KNN{node_no(i)} =  node_no(find(adj(node_no,node_no(i))));
         % q=KNN{i}  
end

%node_no(i) 
%Plotquality(i) 

% % %Drawing decision diagrams
% plot(Plotquality,Plotdegree,'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
% xlabel ('Quality')
% ylabel ('Degree')
% hold on
% for i=1:NN
%     c = num2str(node_no(i));
%     c=[' ',c];
%     text(Plotquality(i),Plotdegree(i),c)
% end
% node_no(i)Represents the node number in the original graph corresponding to the i-th node after folding
%Get minimum threshold
%rect = getrect(1);
% qualitymin=rect(1);
% degreemin=rect(2);
% clf;


%Automatic determination of community number based on difference
garma = Plotquality.*Plotdegree; 
[ordergrama,~] = sort(garma,'descend');
SN =NN;
if NN>100
SN = floor(sqrt(NN));  %Only the first sqrt (NN) points are considered
end

chafen1=zeros(1,SN-1);
for i=1:SN-1
    chafen1(i) = abs(ordergrama(i)-ordergrama(i+1));
end
chafen2=zeros(1,SN-2);
for i=1:SN-2
    chafen2(i) = abs(chafen1(i)-chafen1(i+1));
end

plot([1:SN-2],chafen2);
%plot([1:NN-1],chafen,'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
% lower = ceil(NN*0.05);
% uper = floor(NN*0.95);
[~,maxindex] = max(chafen2(2:SN-2));  
%maxindex =6;
value = ordergrama(maxindex+1);
[~,foldcenterindex] = find(garma >= value);

centers = node_no(foldcenterindex) 
K = length(centers);
fprintf('The number of community centers: %i \n', K);

TT2 = etime(clock,t0)

%=============The third stage: calculate the new similarity between vertices according to the information of nodes and edges owned by common neighbors, and then propagate labels based on this

%%%%%%%%Label propagation
NodeImportant = sum(foldweight(:,:));  
Commonobject =cell(K,K);
CurrentUnlabelNeighbors = cell(K,1);  
Cluster = cell(K,1);
stop = K;  

for i=1:K
     Label(centers(i)) = i;
     Cluster{i} = centers(i);
end
TempCluster =  Cluster; 

while (stop) 
    stop = K;
    CurrentCluster =  TempCluster; 
    TempCluster = cell(K,1);
    
    for i=1:K
        indexid = CurrentCluster{i};  
        numindex = length(indexid); 
        neighindex=[];
        for p=1:numindex
            tempKNN = KNN{indexid(p)}; 
            tempindex = Label(tempKNN) == 0;
            neighindex = union(neighindex,tempKNN(tempindex));
        end
        CurrentUnlabelNeighbors{i} = neighindex;
        
        if (isempty(CurrentUnlabelNeighbors{i}))
            stop = stop - 1;
        end
    end
    
    if (stop == 0)
        break;
    end
    
    %Count the classes each common neighbor belongs to
    commonobject = [];
    Comobject_cluster=[];
       for i= 1:K-1 
        for j= i+1:K
           tempcommonobject = intersect( CurrentUnlabelNeighbors{i}, CurrentUnlabelNeighbors{j});
            if ~isempty(tempcommonobject)
                commonobject = union(commonobject,tempcommonobject); 
                [num,numcol] = size(tempcommonobject);
                if numcol > 1
                    tempcommonobject = tempcommonobject';
                    num = numcol;
                end
                Comobject_cluster = cat(1,Comobject_cluster,[tempcommonobject,repmat(i,num,1)]);
                Comobject_cluster = cat(1,Comobject_cluster,[tempcommonobject,repmat(j,num,1)]);
            end
        end
       end
    %Dealing with common neighbors
    numcommonobject = length(commonobject);
    for i=1:numcommonobject
        maxweight=0;
        maxindex = 0;
        commonindex = Comobject_cluster(:,1) == commonobject(i);
        commoncluster = unique(Comobject_cluster(commonindex,2)); 
        commonclusternum = length(commoncluster);
        for j = 1:commonclusternum   %   commoncluster(j)
            currentcenter = CurrentCluster{commoncluster(j)};   
            
            connectobjectindex = foldweight(currentcenter,commonobject(i)) ~= 0 ;
            connectobjectid = currentcenter(connectobjectindex);
            cn = length(connectobjectid); 
            tempsim = 0;
            for p=1:cn
               tempsim = tempsim + foldweight(commonobject(i),connectobjectid(p));
               commonneigh = intersect(KNN{commonobject(i)},KNN{connectobjectid(p)});
               if ~isempty(commonneigh)
                   tempsim = tempsim + sum(1./NodeImportant(commonneigh)); 
               end
            end
            tempweight = tempsim ;
            if tempweight > maxweight
                maxweight = tempweight;
                maxindex = commoncluster(j);
            end
        end
        Label(commonobject(i)) = maxindex;
        Cluster{maxindex}  = union(Cluster{maxindex},commonobject(i));
        TempCluster{maxindex}  = union(TempCluster{maxindex},commonobject(i));
    end
   
    %Handle objects with neighbors of each class that have no public neighbors
    for i=1:K
        object = setdiff(CurrentUnlabelNeighbors{i},commonobject);
        Label(object) = i;
        Cluster{i} = union(Cluster{i},object);
        TempCluster{i} = union(TempCluster{i},object);
    end
end

%Label compressed objects
% for i=1:NN
%     if Nodelist{node_no(i)} >1
%       Label(Nodelist{node_no(i)}) = Label(node_no(i));
%     end
% end
for i=1:N
    Label(Nodelist{i}) = Label(i);
end

Unlabelobject = find(Label == 0);
if (~isempty(Unlabelobject))
    Label(Unlabelobject) = K+1;
end

Times = etime(clock,t0)