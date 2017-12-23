function out = xdif_pde(t,D)
% This function is used to generate the in-plane diffuison decay factor
% In-plane diffusion is simualted by pdesolver

% Dummy parameters for test
%t = 5e-9;
%D = 0.1;

% ==========================Simulation Method=============================

% Approximation1: Inital carrier distribution is independent in each
% direction.
% n0 = n0x*n0y*n0z;

% Approximation2: Carrier diffuion is indepedent in each direction
% G_total = gx*gy*gz

% Therefore, carrier dynamics in each direction could be calculated
% seperatly and then times together

% Carrier in-plane diffusion is a free diffusion process and is homogenous
% in all the in-plane direction. In a Cartesian system, we consider two
% identical x and y diffusions

% According to doi: 10.1063/1.4944597, Eq.(15), the overall diffusion Green
% function is given as
% G_total = gx(x,x';t)*gy(y,y';t)*gz(z,z';t)*g_tau(t)

% Carrier concentration is calculated as
% n(t) = int{dr'*G_total*n(x',y',z';t0)}
%      \appro int{dr'*gx(x,x';t)*gy(y,y';t)*gz(z,z';t)*g_tau(t)*n0x(x';t0)*n0y(y';t0)*n0z(z';t0)}
%      =int{dx'*gx(x,x';t)*n0x(x';t0)}*int{dy'*gy(y,y';t)*n0y(y';t0)}*int{dz'*gz(z,z';t)*n0z(z';t0)}*g_tau(t)
%      =nx(t)*ny(t)*nz(t)*g_tau(t)
%      =nx(t)*ny(t)*nz_tau(t)
% nz_tau(t) is calculated by the main program

% Here, we try to calculate the diffuision function 
% nx(t)*ny(t) = int{dx'*gx(x,x';t)*n0x(x';t0)}*int{dy'*gy(y,y';t)*n0y(y';t0)}
%             = {int{dx'*gx(x,x';t)*n0x(x';t0)}}^2

% In the x-drection
% n0x \appro  gaussmf(x,[300e-7/2.355 0]), 1-D gassuian distribution
% gx(x,t) = exo(-x^2/4/D/t)/(2*sqart(pi*D*t)), According to doi: 10.1063/1.4944597, Eq.(15)

% nx(t) = conv(gx(x,t), n0x), conv is the convolution

% The above distribution is the total diffusion result, 
% An additional factor is the finite detection volume
% dectV_total = dectVx*dectVy*detecVz
% dectVx = dectVy = gaussmf(x,[300e-7/2.355 0]), two 1-D gassuians
% dectVz = exp(-alpha_pl*z)
% Signal = int{dr*n^2*dectV_total}
%        = int{dr*[nx(t)]^2*[ny(t)]^2*[nz_tau(t)]^2*dectVx*dectVy*detecVz}
%        = int{dx*[nx(t)]^2*deteVx}*int{dy*[ny(t)]^2*deteVy}*int{dz*[nz_tau(t)]^2*deteVz}
%        = Sx*Sy*Sz
% Sz is calculated by the mian program

% Therefore, the total decay factor to the carrier concentration is
% (conv(gx(x,t), n0x)*deteV)^2
% This decay factor will be simulated at each time point
%==========================================================================

% In-plane dx
dx = sqrt(max(diff(t))*D*2);
if dx < 1e-7
    dx = 1e-7;
end
% Generate x axis
x = -900e-7:dx:900e-7;
if length(x) < 5
    x = lispace(-900e-7,900e-7,5);
end
% Simulte 1-D diffuion, no surface recombination at the boundary (bad approximation)
n_x_t =  pdepe(0, ...
               @(x,t,u,DuDx)pdex1pde_x(x,t,u,DuDx,D), ...
               @(x)pdex1ic_x(x), ...
               @(xl,ul,xr,ur,t)pdex1bc_x(xl,ul,xr,ur,t), ...
               x,t);
% Define 1-D detection volume
deteVx = gaussmf(x,[300e-7/2.355 0]);
% Simulate the x component of the signal
temp = bsxfun(@times,deteVx,n_x_t.^2);
Sx = sum(temp,2);
% Square for both direction
out = Sx.^2;
% Normalize the decay factor
out = out./max(out);
end

