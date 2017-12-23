clc
clear
postdata = data_processing(100e-9, 0, 0.5);

%% Control parameters
% Consider in-plane diffusion correction (full 3D diffusion)
in_plane_tag = 0;
% Plot the temporary results during the fitting
plot_tag = 0;
% Save the result to the "output/" folder
save_tag = 1;
% Refine the fitting result based on the last fitting
refine_tag = 1;
if ~exist('output/p.mat','file')
    error('No previous fitting result exist, set the refine_tag = 0')
end
% Start point index of the PL signal, only those after this point will be
% used for the error estimation
pro_index = 1;
% Thickness of the sample
L = 300e-7;            % unit cm
% Number of global parameters
num_global = 3;
% Centralization factor for each fitting parameter
c_S1 = 1e3;
c_S2 = 0;    % set 0 to disable this parameter
c_Amp = 1;
c_offset = 1;
c_D = 1e-2;
c_tau = 1e-6;
c_small_alpha = 0.5; % set the acutal small_alpha value here !
% Fitting boundary and initial value of each parameter, set the lower
% boundary equal to the upper bounary to freeze the parameter
l_S1 = 0.05;            p0_S1 = 0.3;            u_S1 = 0.5;
l_S2 = 0;               p0_S2 = 0;              u_S2 = 0;
l_D = 1e-3;             p0_D = 1;               u_D = 10;
                                                u_Amp = 2;  % Amp use the exp_pl(1) as ld and p0
% no ld, p0 and ud for offset, read from exp_pl
l_tau = 1e-3;           p0_tau = 2.5;           u_tau = 100;
l_small_alpha = 1;      p0_small_alpha = 1;     u_small_alpha = 1;
                     

%% Initialize the parameters
% Stability condition for 1D diffusion dt <= 1/2*dx^2/D
% Means dx >= sqart(2*D*dt)
Dtemp = 1e-4;
dz =  sqrt(max(diff(postdata(1).t))*Dtemp*2);
% Generate the thickness axis
z = (0:dz:L)';
if length(z) < 5
    z = linspace(0,L,5)';
end
% Read number of data set
data_vol = length(postdata);

% Define the fitting parameters, all the parameters below will be used in
% the fitting, the frist num_global paramters will be used as the global
% fitting paramters
% 'tau': bulk liftime
% 'D': diffusion coeffient
% 'S1': surface recombination velocity of the excitated surface
% 'S2': surface recombination velocity of the backside surface
% 'Amp': amplitude of the curve
% 'offset': y-offset of the curve
% 'small_alpha': smooth parameter of the initialized carrier distribution.
% See pdex1ic.m for more details.
pname= {'tau','Amp','offset','D','S1','S2','small_alpha'};
central_factor = struct('S1',c_S1,...
                        'S2',c_S2,...
                        'Amp',c_Amp,...
                        'offset',c_offset,...
                        'D',c_D,...
                        'tau',c_tau,...
                        'small_alpha',c_small_alpha);

% All the global and priviate fitting parameters will be save in the same 
% vector p, use the index saved in the p_pointer to indicate the meaning of
% these parameters
pointer_index = 1; 
% global index
for i = 1:data_vol
    for j = 1:num_global
        p_pointer(i).(pname{j}) = pointer_index;
        pointer_index = pointer_index + 1;
    end
end
% private index
for i = num_global+1 : length(pname)
    for j = 1:data_vol
        p_pointer(j).(pname{i}) = pointer_index;
    end 
    pointer_index = pointer_index+1;
end
pointer_index = pointer_index-1;

% Shuffle the random seed
rng shuffle
% Initlize the parameters by random number
if refine_tag == 0
    for i = 1 : pointer_index
        p0(i) = rand*10;
    end
else
    load('output/p.mat');
    p0 = p;
end

% Define the public lower boundary
LD = 0*p0+eps;
% Define the public higher boundary
UD = LD+Inf;

for i = 1:data_vol
    
    sample_names{i} = postdata(i).name;
    
    UD(p_pointer(i).S1) = u_S1;
    LD(p_pointer(i).S1) = l_S1;
    
    UD(p_pointer(i).Amp) = u_Amp-postdata(i).offset;
    LD(p_pointer(i).Amp) = postdata(i).exp_pl(1)-postdata(i).offset;
    
    UD(p_pointer(i).offset) = postdata(i).offset;
    LD(p_pointer(i).offset) = postdata(i).offset;
    
    UD(p_pointer(i).D) = u_D;
    LD(p_pointer(i).D) = l_D;
    
    UD(p_pointer(i).tau) = u_tau;
    LD(p_pointer(i).tau) = l_tau;
    
    UD(p_pointer(i).small_alpha) = u_small_alpha;
    LD(p_pointer(i).small_alpha) = l_small_alpha;
     
    if refine_tag == 0
     p0(p_pointer(i).S1) = p0_S1;
     p0(p_pointer(i).Amp) = postdata(i).exp_pl(1)-postdata(i).offset;
     p0(p_pointer(i).offset) = postdata(i).offset;
     p0(p_pointer(i).D) = p0_D;
     p0(p_pointer(i).tau) = p0_tau;
     p0(p_pointer(i).small_alpha) = p0_small_alpha;
    end

end

% Fitting parameters
options = optimoptions('fmincon', ...                      
                       'Algorithm','sqp',...
                       'MaxFunEvals',5000, ...                       
                       'Display','iter-detailed', ...
                       'TolFun',1e-6, ...
                       'TolX',1e-6, ...
                       'UseParallel',0);

% Test using the inital guesses                
test = dummy_global(p0,p_pointer,central_factor,postdata,z,pro_index,in_plane_tag,1,1);
disp('Check the dummy curves, press enter to continue...');
pause
disp('Run fitting...')
% main 
p = fmincon(@(p)dummy_global(p,p_pointer,central_factor,postdata,z,pro_index,in_plane_tag,plot_tag,0),p0,[],[],[],[],LD,UD,[],options);
disp('Show results...')
dummy_global(p,p_pointer,central_factor,postdata,z,pro_index,in_plane_tag,1,1);

for i = 1:data_vol
    result(i) = struct('S1',          p(p_pointer(i).S1)*central_factor.S1, ...
                       'S2',          p(p_pointer(i).S2)*central_factor.S2, ...
                       'D',           p(p_pointer(i).D)*central_factor.D, ...
                       'tau',         p(p_pointer(i).tau)*central_factor.tau, ...
                       'Amp',         p(p_pointer(i).Amp)*central_factor.Amp, ...
                       'offset',      p(p_pointer(i).offset)*central_factor.offset, ...
                       'small_alpha', p(p_pointer(i).small_alpha)*central_factor.small_alpha);
end

result_tab = struct2table(result);
result_tab.Properties.RowNames = sample_names
disp('Save results...')
writetable(result_tab,'output/fitting_result.csv');