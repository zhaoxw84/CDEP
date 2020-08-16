function CDEPTest(dataname,adj,adjclass)

if nargin<3
        adjclass = 0;
end
   [VV,times]= fold(adj);
   [label] = LabelNorm(VV);
   NMI =0;
   if adjclass ~=0
      NMI = nmi(adjclass',label');
   end
   Q=modularity_metric(label,adj);
   TK = length(unique(label));
   file_name = ['Fold_',dataname, '_Result.mat'];
    save(file_name, 'NMI', 'Q','label','TK','times');   
end