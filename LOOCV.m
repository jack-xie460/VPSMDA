clc;clear; 
load('.\mirna_disease495with383\You_dataset.mat');
interaction = miRNA_disease_Y;
MFS = miRNA_Function_S;
DSS = disease_Sem_S;

DSS=DSS-diag(diag(DSS)-1);
logweight1 = -6;
logweight2 = 6;
logweight3 = 0;
[nm,nd] = size(interaction);
interaction = interaction';
DSSP = ones(nd,nd);
MFSP = ones(nm,nm);
DSSP(DSS == 0) = 0;
MFSP(MFS == 0) = 0;

[kd0,km0] = logsimilarity(interaction,nd,nm,logweight1,logweight2,logweight3);
[sd0,sm0] = integratedsimilarity(MFS,MFSP,DSS,DSSP,kd0,km0);
w0 = neighborhood_Com(sd0,340);
v0 = neighborhood_Com(sm0,22);
sd0= sd0.*w0;
sm0= sm0.*v0;
[y] = LapRLS_mb(sm0,sd0,interaction', 2^(-5),1);
y_1= y';
index1 = find(interaction == 1);

 for i = 1:length(index1)
     
     fprintf('i = %f\n',i); 
     y_train = interaction;
     y_train(index1(i)) = 0;
     [kd,km] = logsimilarity(y_train,nd,nm,logweight1,logweight2,logweight3);
     [sd,sm] = integratedsimilarity(MFS,MFSP,DSS,DSSP,kd,km);
     w = neighborhood_Com(sd,230);
     v = neighborhood_Com(sm,30);
     sd= sd.*w;
     sm= sm.*v;
     [F_1] = LapRLS_mb(sm,sd,y_train', 2^(-5),1);
     F_1 = F_1';
     y_1(index1(i))= F_1(index1(i));
 end       
auc = roc_1(y_1(:),interaction(:),'red');