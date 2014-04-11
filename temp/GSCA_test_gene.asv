function[pval]=GSCA_test_gene(gene,path,nperm,Z0,W0,B0,N,lambda_w,lambda_b)   

%gene=vector of n SNPs
%path=index the path coefficient to be tested

        
N = size(Z0,1);
n = size(gene,2);

bindex = find(B0 == 1);


MatBperm = zeros(nperm,length(bindex));

vec_FIT = zeros(nperm,1);
vec_FIT_m = zeros(nperm,1);
vec_FIT_s = zeros(nperm,1);
vec_AFIT = zeros(nperm,1);
vec_GFI = zeros(nperm,1);
vec_SRMR = zeros(nperm,1);

%Estimate the original parameters
[W,B,vecW,vecB,FIT,FIT_M,FIT_S,AFIT,GFI,SRMR] = GSCAestim(Z0, W0,B0,N,lambda_w,lambda_b);

%Permutations
for perm = 1:nperm

    Z0perm=Z0;
    indices=randsample(1:N,N,false);
      
    for i=1:n 
      Z0perm(:,gene(i))=Z0(indices,gene(i));
    end
 
    [W_perm,B_perm,vecW_perm,vecB_perm,FIT_perm,FIT_M_perm,FIT_S_perm,AFIT_perm,GFI_perm,SRMR_perm] = GSCAestim(Z0perm,W0,B0,N,lambda_w,lambda_b);        
                           
     
  % may be optimize storage into matrix form, amy be store only ones in
  % question (tested ones)
    MatBperm(perm,:) = B_perm(bindex)';

  
 end   
% Will need to let user specify latent variables rather than b's, pathway
% coefficients
%p-value calculation 
pval=mean(abs(vecB(path,1))<abs(MatBperm(:,path)));




  