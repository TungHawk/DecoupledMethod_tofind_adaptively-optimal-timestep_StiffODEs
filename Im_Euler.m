function [Xnew,niter,Jpsi] = Im_Euler(X)
tol = 1e-8;
miter = 60;
niter = 1;
Xb = X;
while 1
    [psi,Jpsi] = func_psi(X,Xb);
    d = -Jpsi\psi; %compute Newton direction
    Xnew = X + d; %compute the updating state
    X = Xnew;
    
    fprintf('\n    iter                 |psi|  ');
    fprintf('\n  xxxxxxxx            xxxxxxxxxxx');
    fprintf('\n  %7.1e    %12.5e',niter,norm(psi));
    
    if norm(psi) < tol
        fprintf('\nApproximate solution found.\n');
        niter
        return
    end
    if  niter >=  miter
        fprintf('\nMaximum iteration number reached.\n');
        return
    end
  
    niter = niter+1;
end