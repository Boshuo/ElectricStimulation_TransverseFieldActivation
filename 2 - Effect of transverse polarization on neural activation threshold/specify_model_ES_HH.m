function [solver, stimulation, cable] = specify_model_ES_HH(modprmtr)

PW_vec = [ kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),1e1];
% PW in ms, 1 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane); 25 points (4 decades)

if strncmp(modprmtr.model_name,'UF_',3)                  % Uniform Field
    if strncmp(modprmtr.model_name(4:end),'Trm',3)       % Axon TeRMinal
        mod_str = 'Axon terminal';
        L_straight = 1;         % Axon length, 20 mm in cm
        modprmtr.mod = 1;
        % Parameters: 600 = 24 * 25
        alpha_vec = [0:5:75,78:3:87,88,89,89.5,90]*pi/180;        
        % various interval, 24
    elseif strncmp(modprmtr.model_name(4:end),'Bnd',3)  % Axonal BeND
        mod_str = 'Bend axon';
        L_straight = 2;         % Axon length, 20 mm in cm
        L_flank = 1;
        R_bend = 0.05;          % Bend radius, 0.5 mm in cm;  curve section length pi/2 mm ~= 1 lambda (0.8 mm)
        if strcmp(modprmtr.model_name(7:end),'O')        % Orthodromic propagation
            modprmtr.antidromic = 0;
            mod_str =  [mod_str,', orthodromic'];
        elseif strcmp(modprmtr.model_name(7:end),'A') % Antidromic propagation
            modprmtr.antidromic = 1;
            mod_str =  [mod_str,', antidromic'];
        end
        modprmtr.mod = 2;  
        % Parameters: 600 =  24 * 25;
        alpha_vec = [0:15:60,70:10:90,95:5:105,108:115,117,120:5:135]*pi/180;
        % various, 24
    end
    [ALPHA,PPWW] = ndgrid(alpha_vec,PW_vec);
    
    Alpha = ALPHA(modprmtr.id);
    PW = PPWW(modprmtr.id);
    log_str = sprintf('Field orientation in degrees:\t%2.1f',Alpha*180/pi);
    title_str = ['$$ \alpha=', num2str(Alpha*180/pi,'%2.1f'),'\: \rm{^{\circ}}$$; '];
    
    dt = 2e-3;              % Time step 2 us, in ms
else                                            % electrode stimulation
    if strcmp(modprmtr.model_name,'PE_A')                % Point Electrode
        modprmtr.isdisk = 0;
        mod_str = 'Point electrode';
    elseif strcmp(modprmtr.model_name,'DE_A')            % Disk Electrode
        modprmtr.isdisk = 1;
        R_disk = 100e-4;        % 200 um diameter disk 
        mod_str =  'Disk electrode';
    end
    modprmtr.mod = 3;
    % Parameters: 475 =  19 * 25;
    H_vec = [ [3.1,4.5,6.8]*1e-3, kron([1e-2,1e-1],[1,1.4,2.1,3.1,4.5,6.8]),[1,1.4,2.1,3.1] ];     
    % Electrode-axon distance in cm, 31 um to 3 cm; 19 points (3 decades)
    [HHH,PPWW] = ndgrid(H_vec,PW_vec);
    H = HHH(modprmtr.id);
    PW = PPWW(modprmtr.id);
    log_str = sprintf('Axon-Electrode distance in mm:\t%2.4f',H*10);
    title_str = ['$$ H =',  num2str(H*10,'%3.4f'),'\: \rm{mm}$$; '];
    
    dt = 5e-3;              % Time step 5 us, in ms
end

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
sigma_e = 10;               % Extracellular conductivity, in mS/cm;

R_a = 3e-4;              % Axon radius, in cm; 3 um radius
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

% Calculated parameters
freq = 100e-3;                                          % 100 Hz in kHz
lambda_100 = sqrt(R_a * sigma_i / (2*pi*freq*c_m) /2);  % in cm; NOTE: NEURON has lambda_f a sqrt(2) times larger
d_lambda = 0.1;                                         % Interval is determined by d_lambda rule
dz = lambda_100 * d_lambda;                             % For R = 3um, lambda100 = 820.5 um
                                                        % d_lambda = 0.3 -> dz = 246.2 um; d_lambda = 0.1 -> dz = 82.1 um
[m,~,h,~,n,~] = HH(V_rest, T, V_rest);                  % Resting ion channel parameters
g_m_rest = g_Na * m^3 * h + g_K * n^4 + g_L;            % g_m is 0.6773 mS/cm^2 at rest
lambda_DC = sqrt(R_a * sigma_i / g_m_rest /2);          % DC length constant, for R = 3 um, lambda_DC = 0.79 mm

% Emperical parameter from simulation
v = 0.2335;                   % cm/ms; conduction speed 2.335 mm/ms

%% Specify cable
N_theta = 15;               % Discretization of azimuthal angle
d_theta =  pi / N_theta;    % Interval for integration, in radian
theta = linspace( d_theta/2, pi - d_theta/2, N_theta );         % Integration points between 0 and pi, row vector

switch modprmtr.mod
    case 1
        axon_length = L_straight;       % Axon length, 1 cm in cm
        N_comp = ceil(axon_length / dz);
        if (round(N_comp/2) == (N_comp/2))
            N_comp = N_comp + 1;                                % make N_comp odd to have a compartment at center (z = 0)
        end
        axon_length = N_comp * dz;
        z = linspace( -axon_length/2+dz/2, axon_length/2-dz/2, N_comp )';   % axial coordinates, in cm. Row vector
        dz = dz * ones(N_comp,1);
    case 2
        N_straight_comp = ceil(L_straight / dz);
        N_bend_comp = 30;
        d_beta = pi/2/N_bend_comp;
        beta = linspace(-pi/4 + d_beta/2, pi/4 - d_beta/2, N_bend_comp)';
        
        dz = [dz * ones(N_straight_comp,1); d_beta * R_bend * ones(size(beta)); dz * ones(N_straight_comp,1)];
        z = dz(1)/2 + [0; cumsum((dz(1:end-1) + dz(2:end))/2)];

        ind_bend = N_straight_comp + (1 :  N_bend_comp);
        N_comp = N_straight_comp * 2 + N_bend_comp;
        z = z - median(z); 	% Shift origin of coordinates to middle        
    case 3
        axon_min_length = lambda_DC * 5;
        % minimum axon length at least 5 DC length constant,
        % for R = 3 um, axon_min_length = 4.0 mm
        if (5*H >= axon_min_length) 	% If 5*H exceeds minimum axon length (4 mm), axon length set to 5*H: H > 790 um (H >= 1 mm)
            % dz interval is determined by d_lambda rule, as H/5>= 158(200) um > d_lambda (82.1 um)
            dz = lambda_100 * d_lambda;
            axon_min_length = dz * ceil(5*H/dz);
        else                            % Otherwise (H < 790 um (H <= 680 um)), axon length is set to 5 DC length constant
            if H/5 > (lambda_100 * d_lambda)    % If H/5 is more than d_lambda(82.1 um), then set dz to d_lambda  (H>=450 um)
                dz = lambda_100 * d_lambda;
            else
                dz = H/5;                       % Else set dz to H/5 (min: ~4 um )
            end
        end
        N_comp_uni = ceil(axon_min_length/dz);
        axon_length = dz * N_comp_uni;
        N_comp = N_comp_uni * 2 + 1;
        z = linspace( -axon_length, axon_length, N_comp )';     % axial coordinates, in cm. Row vector
        dz = dz * ones(N_comp,1);
end
          
R_i = dz / (sigma_i*pi*R_a^2);       % Axial resistance between nodes, in kOhm
Area = 2 * pi * R_a * dz;        % Element area; cm^2
C_m = c_m * Area;                   % Node membrane capacitance, in uF
ones_A = ones(N_comp, 1);                     % Column vector; empty array
TP_weight = ones(N_comp, N_theta) * d_theta / pi;
cable = struct( 'N_comp',N_comp,...         % Number of compartments
                'z',z,...                   % Center coordinates of compartment, cm (along local longitudinal axis)
                'dz',dz,...                 % Compartment length, cm
                'R', R_a * ones_A,...       % Compartment radius, cm
                'R_i',R_i,...               % Compartment axial resistance, kOhm
                'R_i_left',R_i,...          % Axial resistance to left neighbor, kOhm
                'R_i_right',R_i,...         % Axial resistance to right neighbor, kOhm
                'Area',Area,...             % Compartment Area, cm^2
                'C_m',C_m,...               % Compartment capacitance, uF
                'TP_dim', ones_A,...        % Dimension: 1 for cylindrical; 2 for sphercial
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


switch modprmtr.mod
    case 1
        cable.z_ind_no_act = 1;                                     % Force activation term to zero at antidronic terminal
        z_ind_AP = 1;                                               % Location for AP to reach antidronic terminal for detection
        z_prop = abs(cable.z(N_comp) - cable.z(z_ind_AP));          % Propagation distance from AP initiation point (terminal) to AP detection point
    case 2
        cable.z_ind_no_act = [1,N_comp];                            % Disabling the terminals, including one internode and node inward from terminal
        cable.Area([1,N_comp]) = cable.Area([1,N_comp]) * 100;      % Compartment's area; cm^2. Increase surface area by 100 to reduce polarization near terminal due to undulation
        cable.C_m([1,N_comp])  = c_m * cable.Area([1,N_comp]);      % Compartment's membrane capacitance, in uF

        if modprmtr.antidromic
            z_ind_AP = find(cable.z < (- L_flank - R_bend * pi/2 ), 1,'last');
        else
            z_ind_AP = find(cable.z > ( L_flank + R_bend * pi/2 ), 1,'first');
        end
        z_prop = abs(cable.z(N_comp) - cable.z(1));        
    case 3
        cable.z_ind_no_act = [1,N_comp];                             	% Disabling the terminals
        z_ind_AP = min( [N_comp, max([find(cable.z >= sqrt(3/2)*H*2.5,1,'first'),find(cable.z >= (5 * lambda_DC),1,'first')]) ] );
        % Distance for AP to pass; minimum of 10 mm or 2.5 times distance where
        % hyperpolarization is strongest (at sqrt(3/2) * H = 1.2 * H): 3*H
        z_prop = min(abs(median(cable.z) - cable.z(z_ind_AP))); % Assume AP initiation somewhere close to center of axon
end


%% E-field at cable

switch modprmtr.mod
    case 1
        E = 1;                      %  1 mV/cm
        Ex = E * sin(Alpha) * ones_A;      % Tranverse E-field in mV/cm
        Ez = E * cos(Alpha) * ones_A;      % Axial field in mV/cm
        d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2;
    case 2
       beta = cable.z(ind_bend)/R_bend;
       beta = [-1/4*pi * ones(N_straight_comp,1) ; beta ; 1/4*pi * ones(N_straight_comp,1)];
       
       E = 1;                      %  1 mV/cm
       Ez  = E * cos(1/2*pi - Alpha + beta);
       Ex  = E * sin(1/2*pi - Alpha + beta);
       
       d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2;
    case 3
        I = 1;                                      % 1 uA; scaling performed by waveform
        
        if modprmtr.isdisk
            % Cooridinate system:
            % Electrode surface defines horizontal x-z plane; y axis at
            % electrode center and perpendicular to surface; Electrode-axon
            % distance H measured as y coordinate of axon.
            % x axis perpendicular to axon, positive direction inward
            % y axis positive direction upward
            % z axis with positive direction to right, axon aligns with z
            x = 0;      % axon goes through center of electrode, no lateral offset
            
            [PHI, RR] = cart2pol(cable.z,x);    % cylindrical coordinates with azimuthal angle measured from the z axis
            
            H_hat = H/R_disk;                   % H_hat = eta * xi
            R_hat = RR/R_disk;                  % R_hat = (1 + xi^2)*(1 - eta^2)
            
            sqrt_minus = sqrt( (R_hat-1).^2 + H_hat^2  );
            sqrt_plus  = sqrt( (R_hat+1).^2 + H_hat^2  );
            sqrt_big   = sqrt( (sqrt_minus + sqrt_plus).^2 - 4 ); % 2 * xi
            
            % Wiley and Webster 1982
            Imp_disk = 1 / (4 * R_disk * sigma_e);          % Resistance in kOhm
            V_0 = I * Imp_disk;                             % Primary voltage in mV
            % Wiley and Webster 1982, Eq.10:
            phi_e = V_0 * 2/pi*asin( 2 ./ ( sqrt_plus + sqrt_minus ) );  % sqrt(1+xi^2) = ( sqrt_plus + sqrt_minus )/2 
            
            Er = V_0 *4 / (pi*R_disk) * ( (R_hat+1)./sqrt_plus + (R_hat-1)./sqrt_minus ) ./ ( (sqrt_plus + sqrt_minus) .* sqrt_big );
            Ey = V_0 *4 / (pi*R_disk) * ( (H_hat  )./sqrt_plus + (H_hat  )./sqrt_minus ) ./ ( (sqrt_plus + sqrt_minus) .* sqrt_big );
            
            Ex = sqrt( (Er.* sin(PHI)).^2 + Ey.^2 );    
%             Ez = Er.* cos(PHI);
        else
            I = 1;                                          % 1 uA; scaling performed by waveform
            r = sqrt( H^2 + cable.z.^2 );                   % Distance to point source, in cm
            Ex = I * H/(4*pi*sigma_e) ./ (r.^3);            % Tranverse E-field in mV/cm
            
            phi_e = I / (4*pi*sigma_e) ./ r;                % Extracellular potenial, in mV            
%             Ez = I * cable.z/(4*pi*sigma_e).*r.^(-3);       % Axial field in mV/cm
            % d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2  ;
        end
        d_phi_e = diff(phi_e);
end

stimulation.d_phi_e_left =  [ 0 ; d_phi_e ];    % Finite difference in potenial, in mV
stimulation.d_phi_e_right = [ d_phi_e ; 0 ];    % Finite difference in potenial, in mV

stimulation.ER_TP = kron(   (1+1./cable.TP_dim) .* abs( Ex ) .*...
                             cable.R , cos(theta)  );



%% Time and stimulation waveform

n_bf_start = 5;
t_start = - n_bf_start * dt;                                % Pulse on-set delay, in ms
t_end =   ceil( ( PW  + 5 + z_prop/v)/ dt ) * dt;           % Simulation end time, PW + ~ 5 ms for AP initiation + propagation time
t_vec = ( t_start : dt : t_end);                            % Time vector
if abs(round(PW/dt) - (PW/dt)) > 1e-3                       % If pulse width is not a multiple of dt
    t_vec = sort([t_vec,PW]);                                % Include PW in time vector
end    

stimulation.pulse_shape = zeros(size(t_vec));
stimulation.pulse_shape( (t_vec > 1e-6) & t_vec <= (PW + 1e-6) ) = 1;
stimulation.PW = PW;

%% Solver related parameters
solver.n_theta = N_theta;
solver.t_vec = t_vec;
solver.Temp = T;
solver.V_init = V_rest;
solver.h_func = str2func('simulate_cable_HH');

switch modprmtr.mod
    case 1
        E_rh = 10;          % mV/cm
        E_rh_TP = 5/R_a;	% mV/cm
        t_ch = 0.3;         % ms
        amp_init = E_rh / ( 1 - 2^ ( -PW/t_ch ) );        % mV/cm
        amp_init_TP = E_rh_TP / ( 1 - 2^ ( -PW/t_ch ) );  % mV/cm
        
        if Alpha < pi/2
            cos_alpha = cos(Alpha);
        else
            cos_alpha = cos(pi/2 * 0.95);
        end
        
        amp_init = min( amp_init / cos_alpha  , amp_init_TP) * 0.5 *(1 + randn(1) * 0.02 );
        
        solver.plot_t_intv = [0.25:0.25:solver.t_vec(end)];      % 250 us interval for plotting; 4 points per ms
        
        solver.thresh_find.unit_str = 'V/m';
        solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m
        
        factor = sqrt(2);
    case 2
        E_rh = 100;         % mV/cm
        E_rh_TP = 5/R_a;    % mV/cm
        t_ch = 0.3;        	% ms
        amp_init = E_rh / ( 1 - 2^ ( -PW/t_ch ) );        % mV/cm
        amp_init_TP = E_rh_TP / ( 1 - 2^ ( -PW/t_ch ) );  % mV/cm
        
        if Alpha < pi/2
            cos_alpha = cos(Alpha);
        else
            cos_alpha = cos(pi/2 * 0.95);
        end
        amp_init = min( amp_init / cos_alpha  , amp_init_TP) * 0.5 *(1 + randn(1) * 0.02 );
        
        solver.plot_t_intv = [0.25:0.5:solver.t_vec(end)];      % 250 us interval for plotting; 4 points per ms
        
        solver.thresh_find.unit_str = 'V/m';
        solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m        
        factor = (2)^(1/8);
    case 3
        % set initial search amplitude (based on RMG paper)
        t_ch = 0.36; % chronaxie from RMG; 102+-8 us
        k_rh = 1e5; % uA/cm^2; k is about 100 uA/mm^2 for 0.1 ms in RMG
        k_PW = k_rh / ( 1- 2^ (-PW/t_ch)); %uA/cm^2
        I_rh = 10;  %uA; I_0 is 25 uA for 0.1 ms in RMG
        I_0_PW = I_rh / ( 1- 2^ (-PW/t_ch)); %uA
        
        amp_init = - ( I_0_PW + k_PW * H^2 ); % uA, cathodic current
        amp_init = amp_init * 0.2 * (1 + randn(1) * 0.02 );  % start at reduce amplitude and add some variation
        
        if H <= 0.5
            solver.plot_t_intv = [0.1:0.1:solver.t_vec(end)];               % 00 us interval for plotting; 50 points per ms
        else
            solver.plot_t_intv = [0.25:0.25:4.5,5:0.5:solver.t_vec(end)];	% 100 us interval for plotting; 10 points per ms
        end
        
        solver.thresh_find.unit_str = 'A';
        solver.thresh_find.unit_amp = 1e-6; % 1 uA = 1e-6 A
        factor = sqrt(2);
end

solver.thresh_find.amp_init = amp_init;
solver.thresh_find.amp_th_acc = 0.5e-2;             % 0.5%, accuracy of threshold finding
solver.thresh_find.factor = factor;
solver.thresh_find.range = 10^4; 
solver.thresh_find.phi_AP = 0;                      % mV; threshold definition, phi_m to cross 0 mV
solver.thresh_find.z_ind_AP = z_ind_AP;

end