function out = xdif(t,D)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%t = 5e-9;
%D = 0.1;

%t = [eps;t];

dx = 1e-8;

x = -600e-7:dx:600e-7;

n0 = gaussmf(x,[300e-7/2.355 0]); %inital radial distribution

%figure()
%plot(x,n0)

[X,T] = meshgrid(x,t);

g = exp(-(X).^2/4/D./T)./(2*sqrt(pi*D.*T));

%result = conv(g,n0)*dx; % excited field distribution
result = conv2(1,n0,g,'same')*dx;

%x2 = linspace(-(length(result)-1)/2*dx,(length(result)-1)/2*dx,length(result));

%figure()
%plot(x2,result)

dectV = gaussmf(x,[300e-7/2.355 0]);

%figure()
%plot(x2,result.*dectV)

temp = bsxfun(@times,dectV,result)*dx;

out = sum(temp,2).^2;

out = [sum((n0.*dectV)*dx).^2;out];

end

