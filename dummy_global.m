function err = dummy_global(p,p_pointer,central_factor,postdata,z,pro_index,in_plane_tag,plot_tag,save_tag)

% This programe is used to generate simulated PL decay curve and compare
% with the exprimental results.

% The diffusion of carriers is seperated into the z directon and in-plane
% directions, which are simulated by 1D pdesolver individually

m = 0;

% Read the number of data set, for the global fitting
data_vol = length(postdata);

% Intialize error output
err = 0;

% For each of the data set, generate the simulation curve and calculate err
for i = 1:data_vol
    % Front-side surface recombination velocity
    S1 = p(p_pointer(i).S1)*central_factor.S1;
    % Back-side surface recombination velocity
    S2 = p(p_pointer(i).S2)*central_factor.S2;
    % Amplitude of the curve
    Amp = p(p_pointer(i).Amp)*central_factor.Amp;
    % y-offset of the surve, no offset in default
    offset = p(p_pointer(i).offset)*central_factor.offset;
    % Diffusion coefficient
    D = p(p_pointer(i).D)*central_factor.D;
    % Bulk lifetime
    tau = p(p_pointer(i).tau)*central_factor.tau;
    % Initial distribution factor. Carrier distribution is given as exp(-z*alpha_abs*small_alpha)
    % Used to adjust the inital carrier distribution when t = 0
    small_alpha = p(p_pointer(i).small_alpha)*central_factor.small_alpha;
    
    % Read the exprimental data
    tempdata = postdata(i);
    % Initalize the pesudo carrier concentration.
    N00 = 1;
    % Read the time axis       
    t = tempdata.t;
    % Read the absoption coefficeint at the excitation wavelength
    alpha_abs = tempdata.alpha_abs;
    
    % Solve the pde equation for the z direction
    % The result 2D matrix is the carrier density at each z coordination
    % and t point. It's homogenous in plane.
    n_z_t =  pdepe(m, ...
                   @(z,t,u,DuDx)pdex1pde(z,t,u,DuDx,D,tau), ...
                   @(z)pdex1ic(z,alpha_abs,N00,small_alpha), ...
                   @(xl,ul,xr,ur,t)pdex1bc(xl,ul,xr,ur,t,S1,S2), ...
                   z,tempdata.t);
               
    % Define detection volume in z direction
    deteVz = exp(-tempdata.alpha_pl.*z)';
    temp = bsxfun(@times,deteVz,n_z_t.^2);
    Sz = sum(temp,2);
    
    % Introduce in-plane correction if necessary
    if in_plane_tag == 1
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
    dummy_pl = (dummy_pl/max(dummy_pl)*Amp + offset);   
    % Calculate the error by the LSM
    err = err + sum((tempdata.exp_pl(pro_index:end)-dummy_pl(pro_index:end)).^2);
      
    % Plot the curve if necessary
    if plot_tag == 1
        figure(10)
        subplot(data_vol,1,i)
        % experimental PL
        plot(t(1:end)*1e9,tempdata.exp_pl(1:end),'o')
        hold on
        % Simulated PL
        plot(t(1:end)*1e9,dummy_pl(1:end))
        hold off
        title(postdata(i).name)
        xlabel('Time (ns)')
        ylabel('Normalized PL')
        drawnow;
        % 2D carrier distribution
%         [Z,T] = meshgrid(z,t);
%         figure(i+10)
%         mesh(Z,T,n_z_t)
%         drawnow;
    end
    
    % Save the data if necessary
    if save_tag == 1
        curve_exp = [t(1:end)*1e9,tempdata.exp_pl(1:end)];
        xlswrite(['output/',postdata(i).serise_name],curve_exp,[postdata(i).name,'_exp_data']);
        curve_sim = [t(1:end)*1e9,dummy_pl(1:end)];
        xlswrite(['output/',postdata(i).serise_name],curve_sim,[postdata(i).name,'_sim_data']);
        %csvwrite('output/curve_sim.csv',curve_exp);
        save('output/p.mat','p');
    end

end
end

