function [sol,multiplier,niter] = augmented_lagrangian(tol, max_iter,U)

global a dt eps
global rho mu0
global fun

% Parameters
if nargin < 1
    tol = 1e-6;
end
if nargin < 2
    max_iter = 1000;
end

% Initialization
Xb = U(1);
X = U(1);
v = U(2);
w = U(3);

mu = mu0;


sol = [X;v;w]; %assign the solution of the unconstrained minization problem
%-----------------------------------%
switch fun
    case 1
        niter = 0;
        for k = 1:max_iter
            inner_iter = 0; % Initialize the count for inner iterations
            % Solve unconstrained subproblem for x
            for kk = 1:max_iter
                v_old = v;
                v = max(0,(a*dt*(mu + rho*(1+a*dt-eps-w)))/(1+(a*dt)^2*rho));
                w = max(0, mu/rho + 1 - eps -(1-v)*a*dt);
                mu_old = mu;
                mu = mu + rho*(1+(1-v)*a*dt - eps -w);
                niter = niter + 1;

                if norm(v-v_old) < eps
                   fprintf('\n Convergence de la boucle interne pour v en   %7.1e    iterations ',kk);
                    break;
                end
                if norm(mu-mu_old) < eps
                   fprintf('\n Convergence de la boucle interne  pour mu en   %7.1e    iterations ',kk);
                    break;
                end
                inner_iter = inner_iter + 1; % Increment the count for inner iterations
            end
            niter = niter + inner_iter; % Increment niter by the number of inner iterations

            %%%%%%% test 
            [v_exact, w_exact] = exact_v_w(a,dt,eps);
            
            error_v = norm(v_exact-v);
            error_w = norm(w_exact-w);
            
              
            X_old = X; 
            X = Xb/(1 + (1-v)*a*dt);
            
%             sol = [X;v;w]; %Assign the updating
            
            % Update Lagrange multiplier
%             mu_old = mu;
           
%             multiplier = mu;
            
            
            % Check convergence
            fprintf('\n    iter                 |X-X_old|  ');
            fprintf('\n  xxxxxxxx            xxxxxxxxxxx');
            fprintf('\n  %7.1e    %12.5e',k,norm(X - X_old));
            if norm(X- X_old) < tol
                fprintf('\n Convergence de la boucle externe  en   %7.1e    iterations ',k);
                break;
            end
        end
        
        niter = niter+k;
        sol = [X;v;w];
        multiplier = mu; 
        
    otherwise
        niter = 0;
        for L = 1:max_iter %outer loop
            [f,df1,~,~] = func_f(X);
            for k = 1:max_iter
                v_old = v;
                v = max(0,-df1*dt*(mu + rho*(1-dt*df1-eps-w))/(1 + (df1*dt)^2*rho));
                w = max(0, mu/rho + 1 - eps + (1-v)*df1*dt);
                mu_old = mu;
                mu = mu + rho*(1 - (1-v)*df1*dt - eps -w);
                if norm(v-v_old) < eps
                   fprintf('\n Convergence de la boucle interne pour v en   %7.1e    iterations ',k);
                    break;
                end
                if norm(mu-mu_old) < eps
                   fprintf('\n Convergence de la boucle interne  pour mu en   %7.1e    iterations ',k);
                    break;
                end
                niter = niter + 1;
            end
            X_old = X;
            X = X - (X-Xb-(1-v)*dt*f)/(1-(1-v)*dt*df1);
            % Check convergence
            fprintf('\n    iter                 |X-X_old|  ');
            fprintf('\n  xxxxxxxx            xxxxxxxxxxx');
            fprintf('\n  %7.1e    %12.5e',k,norm(X - X_old));
            if norm(X- X_old) < tol
                fprintf('\n Convergence de la boucle externe  en   %7.1e    iterations ',k);
                break;
            end
        end
        niter = L + niter;
        sol = [X;v;w];
        multiplier = mu;
        
            
end