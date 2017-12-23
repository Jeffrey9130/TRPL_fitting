% This script is used to simulate the core codes
clear
clc

m = 0;

% Initalize the err
err = 0;


for i = 1:1
    
    D = 0.01; 

    tau = 0.001;
      
    Amp = 1;
    
    offset = 0;
       
    S1 = 100;
    
    S2 = 0;
    
    ex_intensity = 10;     % unit A.U.
    L = 600e-7;            % unit cm
                       

    z = linspace(0,L,200)';
    t = linspace(0,50,1000)*1e-9;
    
    small_alpha = 1;      
    
    N00 =1; % unit cm^-3  

    alpha_abs = 150000;
    
% Solve the pde equation for the z direction
    % The result 2D matrix is the carrier density at each z coordination
    % and t point. It's homogenous in plane.
    n_z_t =  pdepe(m, ...
                   @(z,t,u,DuDx)pdex1pde(z,t,u,DuDx,D,tau), ...
                   @(z)pdex1ic(z,alpha_abs,N00,small_alpha), ...
                   @(xl,ul,xr,ur,t)pdex1bc(xl,ul,xr,ur,t,S1,S2), ...
                   z, t);
               
    % Define detection volume in z direction
    deteVz = exp(-10000.*z)';
    temp = bsxfun(@times,deteVz,n_z_t.^2);
    Sz = sum(temp,2);
    
    % Introduce in-plane correction if necessary
    if 1 == 1
        % Calculte the in-plane diffusion decay factor
        % Due to the in-plane diffusion, carrier density in the detection
        % volume will be smaller than that of the otherwise case
        twod_f = xdif_pde(t,D);
        % Apply the decay factor to the carrier distribution
        dummy_pl = Sz.*twod_f;
    else
        dummy_pl = Sz;
    end

    % Amend the PL curve by the amplitude and offset
    dummy_pl = dummy_pl/max(dummy_pl)*Amp + offset;   
    % Calculate the error by the LSM
    %err = err + sum(( exp_pl(pro_index: end_index)-dummy_pl(pro_index: end_index)).^2);
    
    %if plot_tag == 1
        figure(1)
        hold on
        %plot(t(1: end_index)*1e9, exp_pl(1: end_index),'o')
        %hold on
        plot(t*1e9,dummy_pl)
        xlabel('Time (ns)')
        ylabel('Intensity')        
        %hold off
        %drawnow;
    %end

end
