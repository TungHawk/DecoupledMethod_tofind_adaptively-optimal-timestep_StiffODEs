function [psi,Jpsi] = func_psi(X,Xb)

global dt

[f,df1] = func_f(X);
psi = X - Xb - dt*f;
Jpsi = 1 - dt*df1;
   
end