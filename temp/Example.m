

addpath('C:/Users/stepan.grinek/Desktop\GSCA')


filenameData= sprintf('C:/Users/stepan.grinek/Desktop/GSCA/Data_recoded_Z.txt');
      
Z0= load(filenameData);
 


filenameW = sprintf('C:/Users/stepan.grinek/Desktop/GSCA/W_recoded_data.txt');
filenameB = sprintf('C:/Users/stepan.grinek/Desktop/GSCA/B.txt');


W0= load(filenameW);

B0= load(filenameB);
      
Ngen=41;
Nphen=8;
N=1707;

lambda_w=0;
lambda_b=0;
[Westim, Bestim, vecWestim, vecBestim, FIT, FIT_M, FIT_S, AFIT, GFI, SRMR] = GSCAestim(Z0, W0, B0, N, lambda_w, lambda_b); 


%%Tests
nperm=10;

gene=1;
path=1;
pval_b1=GSCA_test_gene(gene,path,nperm,Z0,W0,B0,N,lambda_w,lambda_b) ;  

gene=2;
path=2;
pval_b2=GSCA_test_gene(gene,path,nperm,Z0,W0,B0,N,lambda_w,lambda_b) ;  



gene=[6 7 8 9 10];
path=6;
pval_b6=GSCA_test_gene(gene,path,nperm,Z0,W0,B0,N,lambda_w,lambda_b) ;  




##########################################################

%%Regularized version

%% K-fold cross validation
 LW = [0, 0.1, 0.5, 1, 5, 10, 50, 100];
 LB = [0, 0.1, 0.5, 1, 5, 10, 50, 100];
 K = 5;
 [lambda_w,lambda_b, MATPE] = gsca_ridge_cv(Z0, W0, B0, K, LW, LB);

[Westim,Bestim,vecWestim,vecBestim,FIT,FIT_M,FIT_S,AFIT,GFI,SRMR] = GSCAestim(Z0, W0,B0,N,lambda_w,lambda_b); 
