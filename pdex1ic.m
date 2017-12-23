function u0 = pdex1ic(x,alpha_abs,N00,small_alpha)

%global alpha_abs
%global N00

%temp = 1;

u0 = N00*exp(-(alpha_abs*small_alpha)*x);

end