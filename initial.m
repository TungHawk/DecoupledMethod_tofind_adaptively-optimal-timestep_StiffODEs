function [X0,a,dt,rho,mu0,tol,max_iter] = initial(cases)

global fun

switch fun
    case 1 %linear function f
        switch cases
            case 1 %a>0 and 1+adt>1
                X0 = 0.3;
                dt= 0.1;
                a = 50;
%                 c = 1+a*dt;
                rho = 0.0001;
                mu0 = 1;
                tol = 1e-6;
                max_iter = 100;
            case 2 %a<0 and 0<1+adt<1
                X0 = 0.3;
                dt = 0.4;
                a = -2;
                %                 c = 1+a*dt;
                rho = 0.354;
                mu0 = 0;
                tol = 1e-6;
                max_iter = 1000;
            case 3 %a<0 and -1<1+adt<0
                X0 = 0.8;
                dt = 0.1;
                a = -15;
                % c = 1+a*dt;
                rho = 0.1826;
                mu0 = 0;
                tol = 1e-6;
                max_iter = 100;
            case 4 %a<0 and 1+adt<-1
                X0 = 0.8;
                dt= 0.1;
                a =-30;
%                 c = 1+a*dt;
                rho = 0.004149;
                mu0 = -0.015;
                tol = 1e-6;
                max_iter = 100;
        end
    case 2 %quadratic function f
        switch cases
            case 1 %X0 =1%
                X0 = 1;
                a = 3;
                dt = 0.1;
                rho = 0.0001;
                mu0 = 1;
                tol = 1e-6;
                max_iter = 100;
            case 2 % 0<a*dt<1 
                X0 = 0.8;
                dt= 0.2;
                a = 2;
                rho = 0.0001;
                mu0 = 1;
                tol = 1e-6;
                max_iter = 100;
            case 3 %a*dt>1 and psi'(X0)=1+a*dt*(1-2X0)>0 -> converges to the correct root
                X0 = 0.6;
                dt = 0.4;
                a = 3;
                rho = 0.0001;
                mu0 = 1;
                tol = 1e-6;
                max_iter = 100;
            case 4 %a*dt>1 and psi'(X0)=1+a*dt*(1-2X0)<0 -> converges to the wrong root
                X0 = 0.9;
                dt = 0.7;
                a = 3;
                rho = 0.01;
                mu0 = -0.01;
                tol = 1e-6;
                max_iter = 100;
        end
    case 3 %F(X) = -a*sin(X)
        switch cases
            case 1 %a>0 and a*dt<=1
                X0 = 2;
                dt= 0.2;
                a = 2;
                rho = 0;
                mu0 = 0;
                tol = 1e-4;
                max_iter = 100;
            case 2 %a*dt > 1 and psi'(X_0) > 0 -> has a chance converge to the correct root
                X0 = 1.6;
                dt= 0.2;
                a = 50;
                rho = 0.0298;
                mu0 = 0.102;
                tol = 1e-6;
                max_iter = 100;
            case 3 %a*dt > 1 and psi'(X0)<0
                X0 = 1.7;
                dt= 0.2;
                a = 50;
                rho = 0.01;
                mu0 = -0.1;
                tol = 1e-6;
                max_iter = 100;

        end
    case 4 %f = -a/2*X*(1-X^2)
        switch cases
            case 1 %-1<= a*dt < 0 
                X0 = 0.8;
                dt= 0.2;
                a = -3;
                rho = 0;
                mu0 = 0;
                tol = 1e-6;
                max_iter = 100;
            case 2 %a*dt < -1
                X0 = 0.3;
                dt= 0.2;
                a = -50;
                rho = 0.02;
                mu0 = 0;
                tol = 1e-6;
                max_iter = 100;
            case 3 %a*dt > 0
                X0 = 0.3;
                dt= 0.2;
                a = 50;
                rho = 0.02;
                mu0 = 0;
                tol = 1e-6;
                max_iter = 100;
        end
end