function[Westim,Bestim,vecWestim,vecBestim,FIT,FIT_M,FIT_S,AFIT,GFI,SRMR]=GSCAestim(Z0,W0,B0,N,lambda_w,lambda_b)

% N = number of observations
% Ngen = number of observed Genotypes
% Nphen = number of observed Phenotypes
% J = number of observed variables (Genotypes+Phenotypes)=Ngen+Nphen
% P = number of latent variables  (Genes+Clinical pathways)
% T = total number of variables (= J + P)

% Z = N x J matrix of observed variables: Genotypes first then Phenotypes
% W0 = J x P matrix with 0's and 1's. W0[j,p]=1 indicates an arrow from 
%     the observed variable j to the latent variable p.

% B0 = P x P matrix with 0's and 1's. B0[p1,p2]=1 indicates an arrow from 
%     the latent variable p1 to the latent variable p2.
% f0 = 100000000;
% imp = 100000000;

% V = [eye(J), W]
% A = [C,B]
% No Kronecker products used in estimating all parameters (W,C,& B)

% Ridge penalties are added (no regularization)
% lambda_w = ridge parameter for W   default 0
% lamdbda_b = ridge parameter for B  default 0
% if lambda_w = lambda_b = 0 (no regularization)

%Westim= matrix with weight coefficients estimates
%B.estim= matrix with path coefficients estimates
%vecW.estim=vector with the weight coefficients estimates
%vecB.estim=vector with the path coefficients estimates
%TODO trasnform R.


   [J,P] = size(W0); 
                
    T = J + P;
       
    W = W0;
    B = B0;
    C= zeros(P,J);
    
    for p = 1:P
       windex_p = find(W0(:,p));  %find the non null elements
       if length(windex_p) > 1    
          W0(windex_p,p) = 99;
       end
       
    end

    B0(B0 == 1) = 99;

    windex = find(W0 == 99);
  
    bindex = find(B0 == 99);


    W(windex) = rand(length(windex),1);         
    B(bindex) = rand(length(bindex),1);

       
    V = [eye(J),W];
    % data standardization
    Z = zscore(Z0)*sqrt(N)/sqrt(N-1); 
    
    ZZ = Z;
    if rank(Z'*Z) == J
       Z = chol(Z'*Z);
    end
    sizez = size(Z,1);
    Gamma = Z*W;
    Psi = Z*V;
    
    f0 = 100000000;
    imp = 100000000;
 
    it = 0;             
 
     
    while imp > 0.000001 
          it = it+1;   
          
        

          % Step 1: Update B
          tr_b = 0;
          for p = 1:P
              ee = zeros(1,P);
              ee(p) = 1;
              LL = eye(P);
              LL(p,p) = 0;
              b0 = B0(:,p);
              bindex_p = find(b0 == 99);
              YY = Gamma - Gamma*B*LL;
              B(bindex_p,p) = (Gamma(:,bindex_p)'*Gamma(:,bindex_p) + lambda_b*eye(length(bindex_p)))\(Gamma(:,bindex_p)'*YY*ee');
              tr_b = tr_b + B(bindex_p,p)'*B(bindex_p,p);
          end              
          A = [C,B];
 
          % Step 2: Update W
          tr_w = 0;
          for p = 1:P
              t = J + p;
              windex_p = find(W0(:,p) == 99);
              m = zeros(1,T);
              m(t) = 1;
              a = A(p,:);
              beta = m - a;
              H1 = eye(P);
              H2 = eye(T);
              H1(p,p) = 0;
              H2(t,t) = 0;
              Delta = W*H1*A - V*H2;      
              Zp = Z(:,windex_p);
              theta = ((beta*beta')*(Zp'*Zp) + lambda_w*eye(length(windex_p)))\Zp'*(Z*Delta)*beta';
              zw = Zp*theta;
              theta = sqrt(N)*theta/norm(zw); 
              % handle sign changes
                      
              W(windex_p,p) = theta;
              V(windex_p,t) = theta;
              tr_w = tr_w + theta'*theta;
          end
          Gamma = Z*W;
          Psi = Z*V;
          dif = Psi-Gamma*A;                    
          f = trace(dif'*dif) + lambda_w*tr_w + lambda_b*tr_b;                    
          imp = f0-f;
          f0 = f;
    end
     
        
    
     NPAR = length(windex)+ length(bindex);
     Fit = 1- f/trace(Psi'*Psi);
     Afit = 1 - ((1-Fit)*(N*J)/(N*J - NPAR));        
     dif_m = Z - Gamma*C;
     dif_s = Gamma - Gamma*B;
     Fit_m = 1- trace(dif_m'*dif_m)/(J*N);
     Fit_s = 1 - trace(dif_s'*dif_s)/(P*N);        
     [Gfi, Srmr, COR_RES] = modelfit_mg(ZZ, W, A, J, P, 1, [1,N]); 
     
    
%output values                      
        NITER = it;
        %NPAR;
        FIT = Fit;
        FIT_M = Fit_m;
        FIT_S = Fit_s;
        AFIT = Afit;
        GFI = Gfi;
        SRMR = Srmr;
        
        Westim = W;
        Bestim = B;

        vecWestim = W(windex);
        vecBestim = B(bindex);