function [solver, stimulation, cable] = specify_model_SC_HH(modprmtr)

PW_vec = [ kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),1e1];
% PW in ms, 1 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane); 25 points (4 decades)

if strncmp(modprmtr.model_name,'UF_',3)                  % Uniform Field
    if strncmp(modprmtr.model_name(4:end),'Axon',4)      % Axon
        R_vec = [0.25,0.5,1,2,4,8]*1e-4;        % Compartment radius, in cm;
        mod_str = 'Axon';
        TP_dim = 1;
    elseif strncmp(modprmtr.model_name(4:end),'Soma',4)  % Soma
        R_vec = [     0.5,1,2,4,8,16]*1e-4; 	% Compartment radius, in cm;
        mod_str = 'Soma';
        TP_dim = 2;
    end
end

[RR,PPWW] =  ndgrid(R_vec,PW_vec);
R = RR(modprmtr.id);
PW = PPWW(modprmtr.id);

log_str = sprintf('Radius\t%2.1f um.',R/1e-4);
title_str = ['$$ R=', num2str(R/1e-4,'%2.2f'),'\: \rm{\mu m}$$; '];
    
dt = 2e-3;              % Time step 2 us, in ms

solver.txt.log_txt = {  sprintf('Model type: \t%s', mod_str),...
                        sprintf('Parameter ID:\t%d',modprmtr.id), ...
                        sprintf('Pulse duration in ms:\t%2.4f',PW),...
                        log_str};
solver.txt.fig_title = [mod_str,'; ',...
                        'Parameter $$ID=',num2str(modprmtr.id),'$$; ',...
                        '$$ PW = ',num2str(PW,'%3.3f'),'\: \rm{ms}$$; ',...
                        title_str];
%% Cellular parameters
% Specific to C & PS model (Roth & Basser 1990; Schnabel & Johannes, 2001;
% Neu, 2016a,b)
T = 23.5;       % degree Celcius; Roth & Basser 1990: 18.5 C

% Conductivities of extra- and intra-cellular spaces
sigma_i = 28.2 ;            % Intracellular conductivity, in mS/cm;
% sigma_e = 10;               % Extracellular conductivity, in mS/cm;

c_m = 1;                    % Membrane capacitance, uF/cm^2
% Cellular time constant for cylindrical cell is 13.55 ns per um in radius;
% tau_c = R * c_m * (sigma_i^-1 + sigma_e^-1)
% 3 um radius, tau_c =  40.6 ns; 
% 0.1 us = 2.5*tau_c    ->  reaching 91.8% of IP
% 0.2 us = 5 tau_c      ->  reaching 99.3% of IP

% HH parameters: original Hodgkin-Huxley
V_rest = -70;   % mV; Roth & Basser 1990: -65 mV; Cartee 2000, Rattay & Aberham 1993: -70 mV
E_Na = V_rest + 115;        g_Na = 120;     % mV & mS/cm^2;
E_K  = V_rest - 12;         g_K  = 36;   
E_L  = V_rest + 10.6;       g_L  = 0.3;    
% g_bar = 0.6773 mS/cm^2 at rest -> r_m = 1.476 kOhm*cm^2 
% NOTE early simulations had V_rest = -70.156 due to setting E_L = -60, 
% while it should be E_L =-59.4 

%% Specify cable
N_theta = 15;               % Discretization of azimuthal angle
d_theta =  pi / N_theta;    % Interval for integration, in radian
theta = linspace( d_theta/2, pi - d_theta/2, N_theta );         % Integration points between 0 and pi, row vector

N_comp = 1;
z = 0;
dz = 2 * R;

R_i = dz / (sigma_i*pi*R^2);        % Axial resistance between nodes, in kOhm
Area = 2 * pi * R * dz;             % Element area; cm^2
C_m = c_m * Area;                   % Node membrane capacitance, in uF
ones_A = ones(N_comp, 1);           % Column vector; empty array
if TP_dim == 1
    TP_weight = ones(N_comp, N_theta) * d_theta / pi;
else
    TP_weight = sin(theta) * d_theta / 2;
end
cable = struct( 'N_comp',N_comp,...         % Number of compartments
                'z',z,...                   % Center coordinates of compartment, cm (along local longitudinal axis)
                'dz',dz,...                 % Compartment length, cm
                'R', R * ones_A,...       % Compartment radius, cm
                'R_i',R_i,...               % Compartment axial resistance, kOhm
                'R_i_left',R_i,...          % Axial resistance to left neighbor, kOhm
                'R_i_right',R_i,...         % Axial resistance to right neighbor, kOhm
                'Area',Area,...             % Compartment Area, cm^2
                'C_m',C_m,...               % Compartment capacitance, uF
                'TP_dim', TP_dim,...        % Dimension: 1 for cylindrical; 2 for sphercial
                'TP_weight',TP_weight...    % Integration weights, row vector for each compartment
                );

% Axial resistance; sealed ends are reflected by d_phi_e_left/right = 0 at terminals
cable.R_i_left(1) = inf; 
cable.R_i_left(2:end)    = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(1:end-1) = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(end) = inf; 

% Biophysics of cable
cable.V_rest = V_rest;

% Ion channels
cable.E_Na = E_Na;
cable.E_K  = E_K;
cable.E_L  = E_L;

cable.g_Na = g_Na;
cable.g_K  = g_K;
cable.g_L  = g_L;

Ex = 1;                      %  1 mV/cm

z_ind_AP = 1;                                               % Location for AP to reach for detection
cable.z_ind_no_act = [];

d_phi_e = [];
stimulation.d_phi_e_left =  [ 0 ; d_phi_e ];    % Finite difference in potenial, in mV
stimulation.d_phi_e_right = [ d_phi_e ; 0 ];    % Finite difference in potenial, in mV

stimulation.ER_TP = kron(   (1+1./cable.TP_dim) .* abs( Ex ) .*...
                             cable.R , cos(theta)  );

%% Time and stimulation waveform

n_bf_start = 5;
t_start = - n_bf_start * dt;                                % Pulse on-set delay, in ms
t_end =   ceil( ( PW  + 5)/ dt ) * dt;                      % Simulation end time, PW + ~ 5 ms for AP initiation + propagation time
t_vec = ( t_start : dt : t_end);                            % Time vector
if abs(round(PW/dt) - (PW/dt)) > 1e-3                       % If pulse width is not a multiple of dt
    t_vec = sort([t_vec,PW]);                               % Include PW in time vector
end    

stimulation.pulse_shape = zeros(size(t_vec));
stimulation.pulse_shape( (t_vec > 1e-6) & t_vec <= (PW + 1e-6) ) = 1;
stimulation.PW = PW;

%% Solver related parameters
solver.n_theta = N_theta;
solver.t_vec = t_vec;
solver.Temp = T;
solver.V_init = V_rest;
solver.h_func = str2func('simulate_SC_HH');
                                               
E_rh = 2.5 / R * (1 + randn(1) * 0.05 );	% mV/cm
t_ch = 0.3;         % ms
amp_init = E_rh / ( 1 - 2^ ( -PW/t_ch ) );  % mV/cm

% solver.plot_t_intv = [0.25:0.25:solver.t_vec(end)];      % 250 us interval for plotting; 4 points per ms

solver.thresh_find.unit_str = 'V/m';
solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m

factor = sqrt(2);

                         
solver.thresh_find.amp_init = amp_init;
solver.thresh_find.amp_th_acc = 0.1e-2;         % 0.1%, accuracy of threshold finding
solver.thresh_find.factor = factor;
solver.thresh_find.range = 10^4; 
solver.thresh_find.phi_AP = 0;                  % mV; threshold definition, phi_m to cross 0 mV
solver.thresh_find.z_ind_AP = z_ind_AP;

end