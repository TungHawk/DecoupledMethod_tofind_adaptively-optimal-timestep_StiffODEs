function g = func_g(U,Xb)

global dt eps rho_v rho_theta
global model

switch model
    case 1
        X = U(1);
        v = U(2);
        w = U(3);
        lamb = U(4);
        nu = U(5);
        theta = U(6);

        [f,df1,df2,~] = func_f(X);
        g = [X - Xb - (1-v)*dt*f; 
            1-(1-v)*dt*df1-eps-w;
            lamb*(1-(1-v)*dt*df1) + theta*(1-v)*dt*df2;
            1 + lamb*dt*f - theta*dt*df1 - nu;
            v - max(0,v-rho_v*nu);
            theta - max(0,theta-rho_theta*w)];
        
%         Jg = zeros(6,6);
%         Jg(1,1) = 1 - A*dt*df1;
%         Jg(1,2) = -B*dt*f;
%         Jg(2,1) = -A*dt*df2;
%         Jg(2,2) = -B*dt*df1;
%         Jg(2,3) = -1;
%         Jg(3,1) = theta*A*dt*df3;
%         Jg(3,2) = theta*B*dt*df2;
%         Jg(3,3) = lamb;
%         Jg(3,4) = w + eps;
%         Jg(3,6) = A*dt*df2;
%         Jg(4,1) = -lamb*B*dt*df1 + theta*B*dt*df2;
%         Jg(4,4) = -B*dt*f;
%         Jg(4,5) = -1;
%         Jg(4,6) = B*dt*df1;
%         if v < nu
%             Jg(5,2) = 1;
%         else
%             Jg(5,5) = 1;
%         end
%         if w < theta
%             Jg(6,3) = 1;
%         else
%             Jg(6,6) = 1;
%         end
    case 2
        X = U(1);
        v = U(2);
        w = U(3);
        nu = U(4);
        theta = U(5);

%         Xb = Ub(1);

        [f,df1,df2,df3] = func_f(X);
        fprintf('\n      f      df1          df2         df3');
        fprintf('\n         xxxxxxxx  xxxxxxxx    xxxxxxxx    xxxxxxxx');
        fprintf('\n  %12.5e     %12.5e      %12.5e      %12.5e',f,df1,df2,df3);
        fprintf('\n ');
        g = [X - Xb - (1-v)*dt*f; 
            1 - (1-v)*dt*df1-eps-w;
            1 - theta*dt*df1 - nu;
            min(nu,v);
            min(theta,w)];
        
        Jg = zeros(5,5);
        Jg(1,1) = 1 - (1-v)*dt*df1;
        Jg(1,2) = dt*f;
        Jg(2,2) = dt*df1;
        Jg(2,3) = -1;
        Jg(3,4) = -1;
        Jg(3,5) = -dt*df1;
        
        if v < nu
            Jg(4,2) = 1;
        else
            Jg(4,4) = 1;
        end
        if w < theta
            Jg(5,3) = 1;
        else
            Jg(5,5) = 1;
        end
    case 3
        X = U(1);
        v = U(2);
        w = U(3);
        eta = 0.4;
        A = (1-eta)*v + eta*(1-v);
        B = 1-2*eta;
        %        Xb = Ub(1);
        [f,df1,df2,~] = func_f(X);
        g = [X - Xb - A*dt*f;
            1 - A*dt*df1-eps-w;
            min(v,w)];
        
        Jg = zeros(3,3);
        Jg(1,1) = 1 - A*dt*df1;
        Jg(1,2) = -B*dt*f;
        Jg(2,1) = -A*dt*df2;
        Jg(2,2) = -B*dt*df1;
        Jg(2,3) = -1;
        if v< w
            Jg(3,2) = 1;
        else
            Jg(3,3) = 1;
        end
end
