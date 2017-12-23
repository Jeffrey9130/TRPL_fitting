function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t,S1,S2)

%global S
%global S2

pl = S1*ul;
ql = -1;

%no surface recombination
pr = S2*ur;
qr = -1;


end