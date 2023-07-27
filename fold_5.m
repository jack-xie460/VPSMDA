clc;clear; 
load('.\mirna_disease495with383\You_dataset.mat');
interaction = miRNA_disease_Y;
MFS = miRNA_Function_S;
DSS = disease_Sem_S;

DSS=DSS-diag(diag(DSS)-1);
MFS=MFS-diag(diag(MFS)-1);
logweight1 = 6;
logweight2 = -6;
logweight3 = 0;
K1 = 340;
K2 = 22;
[nm,nd] = size(interaction);
DSSP = ones(nd,nd);
MFSP = ones(nm,nm);
DSSP(DSS == 0) = 0;
MFSP(MFS == 0) = 0;
iteration = 1;
auc = zeros(1,iteration);
nfolds =5;

for iter = 1:iteration
iter
y=zeros(nm,nd);
crossval_idx1 = crossvalind('Kfold',interaction(:),nfolds);
 for fold = 1:nfolds
     y_train = interaction;
     test_idx1  = find(crossval_idx1==fold);
     y_train(test_idx1) = 0;
     [kd,km] = logsimilarity(y_train',nd,nm,logweight1,logweight2,logweight3);
     [sd,sm] = integratedsimilarity(MFS,MFSP,DSS,DSSP,kd,km);
     w = neighborhood_Com(sd,K1);
     v = neighborhood_Com(sm,K2);
     sd= sd.*w;
     sm= sm.*v;
     [F_1] = LapRLS_mb(sm,sd,y_train, 2^(-5),1);
     y(test_idx1)= F_1(test_idx1);
 end 
auc(iter) = roc_1(y(:),interaction(:),'red');
end
auc_ave = mean(auc);
auc_std = std(auc);
