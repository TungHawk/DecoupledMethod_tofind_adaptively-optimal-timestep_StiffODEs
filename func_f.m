function [f,df1,df2,df3] = func_f(X)

global fun a

switch fun
    case 1
        f = -a*X;
        df1 = -a;
        df2 = 0;
        df3 = 0;
    case 2
        f = -a*X*(1-X);
        df1 = -a*(1-2*X);
        df2 = 2*a;
        df3 = 0;
    case 3
        f = -a*sin(X);
        df1 = -a*cos(X);
        df2 = a*sin(X);
        df3 = a*cos(X);
    case 4
        f = -a/2*X*(1-X^2);
        df1 = -a/2 + 3*a/2*X^2;
        df2 = 3*a*X;
        df3 = 3*a;
    
end