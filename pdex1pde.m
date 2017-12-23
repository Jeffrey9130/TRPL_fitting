function [c,f,s] = pdex1pde(x,t,u,DuDx,D,tau)

%global D
%global tau

c = 1;
f = D*DuDx;
s = -u/tau;
end