function  [exact] = ExactSolution(X0,t)
global fun a
switch fun
    case 1
        exact = X0*exp(-a*t);
    case 2
        exact = X0./(X0+(1-X0)*exp(a*t));
    
    case 3
        if X0 >= 0
            c = 2/(1-cos(X0)) - 1;
            exact = acos(1 - 2./(1+c*exp(2*a*t)));     
        else
            c = 2/(1-cos(X0)) - 1;
            exact = -acos(1 - 2./(1+c*exp(2*a*t)));  
        end
    case 4
        c = X0^2/(1-X0^2);
        if X0 >= 0
            exact = sqrt(c*exp(-a*t)./(1 + c*exp(-a*t)));
        else
            exact = -sqrt(c*exp(-a*t)./(1 + c*exp(-a*t)));
        end
end