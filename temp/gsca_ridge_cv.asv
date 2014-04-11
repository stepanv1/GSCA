% K-fold cross validation for Ridge GSCA
% Heungsun Hwang
% July 11, 2013\
%TODO: check for optimization
function [lambda_w,lambda_b, MATPE] = gsca_ridge_cv(z0, W0, B0, num_folder, LW, LB)
N = size(z0,1);
[J,P] = size(W0);                   
T = J + P;
C0 = W0';
windex = find(W0);
cindex = find(C0);
bindex = find(B0);
W = W0;
C = C0;
B = B0;
%% CROSS VALIDATION STARTS HERE
MATPE = [];                     % PREDICTION ERROR
for b = 1:length(LW)
    lambda_w = LW(b);
    for j = 1:length(LB)
        lambda_b = LB(j);
        PE = 0;
        for nf = 1:num_folder
            % Generate the data indices for the testing data and the training data
            testindex = floor((nf-1)*N/num_folder)+1 : floor(nf*N/num_folder);
            trainindex = setdiff(1:N, testindex);
      
            Z_test = z0(testindex, :);          % testing sample
            Z_train = z0(trainindex, :);        % training sample
            
            % data standardization
            N_train = size(Z_train,1);
            Z_train = zscore(Z_train)*sqrt(N_train)/sqrt(N_train-1);    
            N_test = size(Z_test,1);
            Z_test = zscore(Z_test)*sqrt(N_test)/sqrt(N_test-1); 
            
            if rank(Z_train'*Z_train) == J
               Z_train = chol(Z_train'*Z_train);
            end
            sizez = size(Z_train,1);
            
            W(windex) = rand(length(windex),1);         
            C(cindex) = rand(length(cindex),1); 
            B(bindex) = rand(length(bindex),1);
            V = [eye(J),W];
            Gamma = Z_train*W;
            vecGamma = reshape(Gamma,sizez*P,1); 
            % ALS Algorithm
            it = 0;             
            f0 = 100000000;
            imp = 100000000;
            while it <= 100 && imp > 0.00001 
                  it = it+1;   
                  % Step 1: Update A (C & B)
                  for p = 1:P
                      c0 = C0(p,:);
                      cindex_p = find(c0);
%                       C(p,cindex_p) = Gamma(:,p)'*Z_train(:,cindex_p)/(Gamma(:,p)'*Gamma(:,p));
                      C(p,cindex_p) = Gamma(:,p)'*Z_train(:,cindex_p)/N_train;
                  end
                  Phi = kron(eye(P),Gamma);
                  Phi = Phi(:,bindex);
                  B(bindex) = (Phi'*Phi + lambda_b*eye(length(bindex)))\Phi'*vecGamma; 
                  A = [C,B];
                  % Step 2: Update W
                  tr_w = 0;
                  for p = 1:P
                      t = J + p;
                      windex_p = find(W0(:,p));
                      m = zeros(1,T);
                      m(t) = 1;
                      a = A(p,:);
                      beta = m - a;
                      H1 = eye(P);
                      H2 = eye(T);
                      H1(p,p) = 0;
                      H2(t,t) = 0;
                      Delta = W*H1*A - V*H2;      
                      Zp = Z_train(:,windex_p);
                      theta = ((beta*beta')*(Zp'*Zp) + lambda_w*eye(length(windex_p)))\Zp'*(Z_train*Delta)*beta';
                      zw = Zp*theta;
                      theta = sqrt(N_train)*theta/norm(zw);        
                      W(windex_p,p) = theta;
                      V(windex_p,t) = theta;
                      tr_w = tr_w + theta'*theta;
                  end
                  Gamma = Z_train*W;
                  Psi = Z_train*V;
                  dif = Psi - Gamma*A;                    
                  f = trace(dif'*dif) + lambda_w*tr_w + lambda_b*(B(bindex)'*B(bindex));                    
                  imp = f0-f;
                  f0 = f;
                  vecGamma = reshape(Gamma,sizez*P,1);
           end
           Gamma = Z_test*W;
           Psi = Z_test*V;
           dif = Psi - Gamma*A;         
           PE = PE + trace(dif'*dif);    
       end                                         % K-fold CV per lambda
       lambda = [lambda_w,lambda_b];
       MATPE = [MATPE; [lambda, PE]];
    end                                            % LB
end                                                % LW
[~, minindex] = min(MATPE(:,3));
lambda_w = MATPE(minindex,1);
lambda_b = MATPE(minindex,2);



