clear all
close all
clc

global fun
global a dt eps rho mu0

eps = 0.1;
fun = 3; 
% fun = 1 - linear function f
% fun = 2 - a quadratic function f
% fun = 3 - a trigonometric function f
% fun = 4 - a polynomial degree 3

% Initialization
% [X0,a,dt,rho,mu0,tol,max_iter] = initial(2);
% input variable is the case for each value of fun
% fun = 1 - input is from 1 to 4
% fun = 2 - input is from 1 to 4
% fun = 3 - input is from 1 to 3
% fun = 4 - input is from 1 to 2

X0 = 1.7;
dt= 0.2;
a = 2;
rho = 0.01;
mu0 = 0.1;
tol = 1e-6;
max_iter = 100;





%--------------------Discretizing the time---------------------%
ti=0; tf=10;
t = ti:dt:tf;

%-----------------Compute the exact solution------------------------------%

X_true = ExactSolution(X0,t);

%------------------------Implicit Euler scheme----------------------------%
X_impl = zeros(1,length(t));
X_impl(1) = X0; 
iter_impl = [];

for i = 2:length(t)
   [X_impl(i),niter] = Im_Euler(X_impl(i-1));
   iter_impl = [iter_impl;niter];
end
ReError_impl = 100*(abs(X_true-X_impl))./(1+abs(X_true)); %relative error


%------------------------Decoupling method method-------------------------%
U_NewIterative = zeros(3,length(t));

U_NewIterative(:,1) = [X0;1;1-eps]; %initial value
iter_NewIterative = [];


for i = 2:length(t)
   [U_NewIterative(:,i),multipliers,niter] =  augmented_lagrangian(tol, max_iter,U_NewIterative(:,i-1));
   fprintf('\n');
   disp('Optimal solution:');
   disp(U_NewIterative(:,i));
   disp('Lagrange multiplier:');
   disp(multipliers);
   
   %Number of iterations at each time step
   iter_NewIterative = [iter_NewIterative;niter]; %sum of inner and outer iterations
end

X_NewIterative = U_NewIterative(1,:); % numerical solution of the ODE
v_NewIterative = U_NewIterative(2,:); % reduction factor of the time step
ReError_NewIterative = 100*abs(X_true-X_NewIterative)./(1+abs(X_true)); %relative error

%----------------------Plot the numerical results-------------------------%
folder = 'D:\PUF\Thesis\Figures\Decoupled_ver2';

figure(1)
hold on
plot(t,X_true,'-g','LineWidth',2)
plot(t,X_NewIterative,'-.xb')
plot(t,X_impl,'-om')
xlabel('Time - t')
ylabel('X')
axis([ti tf min(min(X_impl),min(X_true)) max(max(X_impl),max(X_true))]);
legend('Exact solution','Decoupling method','Euler implicit')
title({'Comparison between exact solution',...
    'and numerical solution'})
% saveas(gcf, fullfile(folder, '28.png'));


figure(2)
hold on
plot(t(2:end),iter_NewIterative,'-xb','LineWidth',1)
plot(t(2:end),iter_impl,'-om','LineWidth',1)
xlabel('Time - t')
ylabel('Number of iterations')
legend('Decoupling method','Euler implicit')
title('Number of iteration at each time step')
% saveas(gcf, fullfile(folder, '29.png'));

figure(3)
hold on
plot(t,v_NewIterative,'-xb','LineWidth',1)
xlabel('Time - t')
ylabel('The value of v')
legend('Decoupling method')
title('Evolution of the time-step reduction factor')
% saveas(gcf, fullfile(folder, '30.png'));


% figure(4) %figure of the graph of the implicit equation
% X = -10:0.2:10;
% hold on
% plot(X,func_psi(X,X0),'-r','LineWidth',1)
% plot(X,zeros(size(X)),'-k')
% plot([0 0], ylim, '-k')  % Plotting x=0
% xlabel('$X$', 'Interpreter', 'latex')
% ylabel('$\psi$', 'Interpreter', 'latex')

