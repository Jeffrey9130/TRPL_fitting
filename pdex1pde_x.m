function [c,f,s] = pdex1pde_x(x,t,u,DuDx,D)

% This function is used to define the in-plane diffuions equation, 
c = 1;
f = D*DuDx;
s = 0;
end