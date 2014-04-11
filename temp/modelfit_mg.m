function[GFI, SRMR COR_RES] = modelfit_mg(Z0, W, T, nvar, nlv, ng, case_index) 
% MODEL FIT MEASURES
gfi_1 = 0;
gfi_2 = 0;  
srmr_1 = 0;
kk = 0;
ss = 0;
COR_RES = [];  % correlation residual matrix
for g = 1:ng
    k = kk + 1;
    kk = kk + nvar;
    s = ss + 1;
    ss  = ss + nlv;
    
    zz = Z0(case_index(g,1):case_index(g,2),k:kk);
    w = W(k:kk, s:ss);
    v = [eye(nvar), w];
    t = T(s:ss,:);
    omega = v - w*t;
    ee = zz*omega;
        
    samcov = cov(zz);                            % sample covariances for each group
    samcorr = corrcoef(zz);
    precov = (omega*omega')\omega*diag(var(ee))*omega'/(omega*omega');  % predicted covariances for each group
%     precov = (omega*omega')\omega*cov(ee)*omega'/(omega*omega');  % predicted covariances for each group
    COV_RES = samcov - precov;
    prerij = precov;
    for i = 1:nvar
        for j = 1:nvar
            prerij(i,j) = precov(i,j)/sqrt(precov(i,i)*precov(j,j));
        end
    end
    srmr = 0;
    for i = 1:nvar
        for j = 1:nvar
            if j > i
               corr_residual = (samcorr(i,j) - prerij(i,j)).^2;
               srmr = srmr + corr_residual;
            end
        end
    end  
    srmr_1 = srmr_1 + srmr;
    gfi_1 = gfi_1 + trace(COV_RES.^2);
    gfi_2 = gfi_2 + trace(samcov.^2);  
    COR_RES = [COR_RES;samcorr - prerij];
end
nvar_tot = ng*nvar;
srmr_2 = nvar_tot*(nvar_tot+1)/2;
SRMR = sqrt(srmr_1/srmr_2);         % Standardized root mean square residual
GFI = 1 - (gfi_1/gfi_2);            % GFI-ULS

                                                   
     
