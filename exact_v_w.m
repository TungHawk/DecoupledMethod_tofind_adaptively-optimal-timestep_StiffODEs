function [v,w] = exact_v_w(a,dt,eps)

if a*dt >=(eps-1)
    v = 0;
    w = 1+a*dt-eps;
else 
    v = 1+(1-eps)/(a*dt);
    w = 0;
end 
