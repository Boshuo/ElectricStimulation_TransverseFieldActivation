function [solver, stimulation, cable] = specify_model_ES_RMG(modprmtr)

PW_vec = [ kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),1e1];
% PW in ms, 1 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane); 25 points (4 decades)

if strncmp(modprmtr.model_name,'UF_',3)                  % Uniform Field
    if strncmp(modprmtr.model_name(4:end),'Trm',3)       % Axon TeRMinal
        mod_str = 'Axon terminal';
        if strcmp(modprmtr.model_name(7),'B')            % Bare terminal
            modprmtr.has_synapse = 0;
        elseif strcmp(modprmtr.model_name(7),'S')        % with Synaptic bouton
            modprmtr.has_synapse = 1;
            mod_str =  [mod_str,' with synapse'];
        end
        modprmtr.mod = 1;
        % Parameters: 600 = 24 * 25
        alpha_vec = [0:5:75,78:3:87,88,89,89.5,90]*pi/180;  
        % various interval, 24
    elseif strncmp(modprmtr.model_name(4:end),'Bnd',3)  % Axonal BeND
        mod_str = 'Bend axon';
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
    [HH,PPWW] = ndgrid(H_vec,PW_vec);
    H = HH(modprmtr.id);
    PW = PPWW(modprmtr.id);
    log_str = sprintf('Axon-Electrode distance in mm:\t%2.4f',H*10);
    title_str = ['$$ H =',  num2str(H*10,'%3.4f'),'\: \rm{mm}$$; '];
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
% Richardson, McIntyre & Grill 2000 (RMG)
T = 37;                     % Temperature, Degree Celcius

% Conductivities of intra- and extra-cellular spaces
sigma_i = 1/0.07;           % Intracellular conductivity, in mS/cm; resistivity: 70 Ohm*cm = 0.07 kOhm*cm
sigma_e = 1/0.5;            % Extracellular conductivity, in mS/cm; resistivity: 500 Ohm*cm = 0.5 kOhm*cm

% Nodal membrane parameters
R_n = 1.65e-4;              % Nodal radius, in cm; 3.3 um diameter
c_n = 2;                    % Nodal membrane capacitance, in uF/cm^2
l_n = 1e-4;                 % Nodal length, in cm; 1 um
% "Cellular" time constant for cylindrical nodal segment
% tau_c = R_n * c_n * (sigma_i^-1 + sigma_e^-1) = 188.1 ns
% 3 * tau_c = 0.56 us  ->  reaching 95.0% of steady state TP
% 5 * tau_c = 0.94 us  ->  reaching 99.3% of steady state TP

% Ion channel parameters for nodes, in mV & mS/cm^2;
V_rest = -82;                       % in MRG 2002: -80;
E_Na = 50;        g_Na = 3000;      g_Nap = 5;  % fast & persistent Na
E_K  = -84;       g_K  = 80;        % Slow K
E_L  = -83.38;    g_L  = 80;        % Leakage
% g_bar = 104.3 mS/cm^2 at rest -> r_m = 0.0096 kOhm*cm^2
% tau_n = 19.2 us

% Internodal membrane parameters
R_in = 3e-4;                % Internodal axon radius, in cm; 6 um fiber diameter (10 um diameter including myelin)
N_lamella = 120;            % 120 myelin lamella; 2 membranes per lamella
c_in = 0.1/2/N_lamella;     % Internodal total myelin capacitance, in uF/cm^2
g_in = 1/2/N_lamella;       % Internodal total myelin conductance, in mS*cm^2
L_in = 1150e-4;             % Internodal length, in cm; 1150 um
N_in = 10;             	% 10 compartments per internode (default for straight fiber)
l_in = L_in / N_in;         % Internodal compartment length, in cm

% Emperical parameter from simulation/RMG paper
v = 6.1;                    % cm/ms; conduction speed 61 m/s = 6.1 cm/ms from RMG paper

%% Specify cable
N_theta = 15;               % Discretization of azimuthal angle
d_theta =  pi / N_theta;    % Interval for integration, in radian
theta = linspace( d_theta/2, pi - d_theta/2, N_theta );         % Integration points between 0 and pi, row vector

switch modprmtr.mod
    case 1
        min_node_uni = 20;                  % increased to 20 
        % Minimum number of nodes in either direction (excluding center node); 41 nodes total
    case 2
        min_node_uni = 50;                  % increased to 50 
        % Minimum number of nodes in either direction (excluding center node);
        % 121 nodes total. Correction to JNE 15 026003: the flanking section for the RMG cable are 6
        % cm long, instead of 2 cm
    case 3
        min_length = H * 5;                 % cm, minimum axon length in either direction
        min_node_uni = max(10, ceil( (min_length - l_n/2)/(l_n + L_in) ) );
        % Minimum number of nodes in either direction (excluding center node); at least 10 nodes (11.5 mm) as given in RMG/MRG    
end
N_node = min_node_uni * 2 + 1;              % Number of nodal compartments
N_comp = N_node + (N_node - 1)* N_in;       % Number of all compartments

switch modprmtr.mod
    case 1
        if modprmtr.has_synapse
            % Synaptic bouton
            R_synps = 4.5e-4;             	% 4.5 um in cm
            N_node = N_node + 1;         	% Additional node for synaptic bouton
            N_comp = N_comp + 1;            % Additional compartment for synaptic bouton
        end
        N_myl_comp = N_in * 10;
        N_comp = N_comp + N_myl_comp * 1; 
    case 2
        % Bent cable of 121 nodes and 120 internodes; flanking straight section has
        % 58 nodes and 58 internode each, bend has 5 nodes and 4 internode
        % Correction to JNE 15 026003: the flanking section for the RMG cable are 6
        % cm long, instead of 2 cm
        N_bend_node = 4;
        N_bend_comp = (N_in + 1) * N_bend_node + 1;
        
        N_myl_comp = N_in * 10;
        
        N_straight_node_l = round((N_node - N_bend_node -1 )/2);
        N_straight_comp_l = (N_in + 1) * N_straight_node_l + N_myl_comp;
        
        N_straight_node_r = N_node - N_bend_node - 1 - N_straight_node_l;
        N_straight_comp_r = (N_in + 1) * N_straight_node_r + N_myl_comp;
        
        N_comp = N_comp + N_myl_comp * 2; 
        R_bend = (l_n + L_in) * N_bend_node / (pi/2);
        ind_bend = N_straight_comp_l + (1 :  N_bend_comp);
    case 3
        N_myl_comp = 0;        
end

Empty_N = zeros(N_node,1);              % Column vector; empty array
Empty_th = zeros(N_node,N_theta);       % Column vector expanded; empty array
node = struct(  'cable_ID',Empty_N,...  % the index of the node within the entire cable
                'TP_dim',Empty_N,...    % Dimension: 1 for cylindrical; 2 for sphercial
                'TP_weight',Empty_th... % Integration weights, row vector for each node
              );
Empty_C = zeros(N_comp, 1);             % Column vector; empty array
cable = struct( 'N_comp',N_comp,...     % Number of compartments
                'z',Empty_C,...         % Center coordinates of compartment, cm (along local longitudinal axis)
                'dz',Empty_C,...        % Compartment length, cm
                'R',Empty_C,...         % Compartment radius, cm
                'R_i',Empty_C,...       % Compartment axial resistance, kOhm
                'R_i_left',Empty_C,...  % Axial resistance to left neighbor, kOhm
                'R_i_right',Empty_C,... % Axial resistance to right neighbor, kOhm
                'Area',Empty_C,...      % Compartment Area, cm^2
                'C_m',Empty_C,...       % Compartment capacitance, uF
                'is_node',Empty_C,...   % Is compartment a node? logical 0/1
                'N_nodes',N_node,...    % Number of nodes
                'node',node...          % Nodal data
              );

          
z_in_rel = l_in * (-0.5 + (1 : N_myl_comp));
ind_in = (1:N_myl_comp) ;
cable.z(ind_in) = z_in_rel;
cable.dz(ind_in) = l_in;                                % Compartments' length, in cm
cable.R(ind_in) = R_in;                                 % Compartments' radius, in cm
cable.R_i(ind_in) = l_in / (sigma_i * pi * R_in^2);     % Compartments' axial resistance, in kOhm
cable.Area(ind_in) = 2 * pi * R_in * l_in;              % Compartments' area; cm^2
cable.C_m(ind_in) = c_in * 2 * pi * R_in * l_in;        % Compartments' membrane capacitance, in uF
cable.is_node(ind_in) = false;                          % Compartments are internode

z_in_rel = l_in * ( (1:N_in) - 0.5 );
% Relavtive axial coordinates of the center of N_in internodal compartments
% within one internode from its left end

for ii = 1 : min_node_uni *2 + 1      % Add node and internode
    % 1 nodal compartments
    ind_n = N_myl_comp + (1 + N_in) * (ii-1) + 1;           % Index of compartment in the cable: (ii-1) nodes and internodes to the left
    cable.z(ind_n) = N_myl_comp * l_in + (l_n + L_in) * (ii-1) + l_n/2;         % Center coordinates of compartment, cm
    cable.dz(ind_n) = l_n;                                  % Compartment's length, in cm
    cable.R(ind_n) = R_n;                                   % Compartment's radius, in cm
    cable.R_i(ind_n) = l_n / (sigma_i * pi * R_n^2);        % Compartment's axial resistance, in kOhm
    cable.Area(ind_n) = 2 * pi * R_n * l_n;                 % Compartment's area; cm^2
    cable.C_m(ind_n) = c_n * 2 * pi * R_n * l_n;            % Compartment's membrane capacitance, in uF
    cable.is_node(ind_n) = true;                            % Compartment is node
        cable.node.cable_ID(ii) = ind_n;                    % The index of the node within the entire cable
        cable.node.TP_dim(ii) = 1;                          % Dimension parameter for cylindrical compartments is 1
        cable.node.TP_weight(ii,:) = d_theta / pi;          % Each patch of the membrane is only a fraction of the total area
    if ii == min_node_uni *2 + 1 
        break                                               % No internode after last nodal compartment
    end

    % N_in internodal compartments
    ind_in = N_myl_comp + ii + N_in * (ii-1) + (1:N_in);	% Index of compartments in the cable: ii nodes and  (ii-1) internodes to the left
    cable.z(ind_in) = N_myl_comp * l_in + l_n * ii + L_in * (ii-1) + z_in_rel;  % Center coordinates of compartments, cm
    cable.dz(ind_in) = l_in;                                % Compartments' length, in cm
    cable.R(ind_in) = R_in;                                 % Compartments' radius, in cm
    cable.R_i(ind_in) = l_in / (sigma_i * pi * R_in^2);     % Compartments' axial resistance, in kOhm
    cable.Area(ind_in) = 2 * pi * R_in * l_in;              % Compartments' area; cm^2
    cable.C_m(ind_in) = c_in * 2 * pi * R_in * l_in;        % Compartments' membrane capacitance, in uF
    cable.is_node(ind_in) = false;                          % Compartments are internode
end

switch modprmtr.mod
    case 1
        if modprmtr.has_synapse
            ii = ii+1;
            ind_n = ind_n + 1;
            cable.z(ind_n) = cable.z(ind_n-1) + cable.dz(ind_n-1)/2 + R_synps;
            cable.dz(ind_n) = R_synps * 2;
            cable.R(ind_n) = R_synps;
            cable.R_i(ind_n) = 0;                           
            % Spherical compartment inside is equipotential, has no axial
            % resistance itself (axial resistance to connected
            % cylindrical compartments is determined those compartments
            cable.Area(ind_n) = 4* pi * R_synps^2;
            cable.C_m(ind_n) = c_n * 4* pi * R_synps^2;
            cable.is_node(ind_n) = true;
            cable.node.cable_ID(ii) = ind_n;
            cable.node.TP_dim(ii) = 2;
            cable.node.TP_weight(ii,:) = sin(theta) * d_theta / 2;  % Different weight for spherical structure
        end
        cable.z_ind_no_act = 1;                                     % Force activation term to zero at antidronic terminal
        z_ind_AP = 1 + N_myl_comp + (1 + N_in) * 5;                 % Location for AP to reach antidronic terminal for detection
        z_prop = abs(cable.z(N_comp) - cable.z(z_ind_AP));          % Propagation distance from AP initiation point (terminal) to AP detection point
        cable.z = cable.z - N_myl_comp * l_in - (min_node_uni * (l_n + L_in) + l_n/2); 	% Shift origin of coordinates to middle
    case 2
        z_in_rel = l_in * (-0.5 + (1 : N_myl_comp));
        ind_in = (1:N_myl_comp) ;
        cable.z(ind_in) = z_in_rel;

        ind_in = N_myl_comp + 1 + (1:N_myl_comp)+ (N_in + 1) * (ii-1) ;
        cable.z(ind_in) = z_in_rel + N_myl_comp * l_in + l_n * ii + L_in * (ii-1);
        cable.dz(ind_in) = l_in;                                % Compartments' length, in cm
        cable.R(ind_in) = R_in;                                 % Compartments' radius, in cm
        cable.R_i(ind_in) = l_in / (sigma_i * pi * R_in^2);     % Compartments' axial resistance, in kOhm
        cable.Area(ind_in) = 2 * pi * R_in * l_in;              % Compartments' area; cm^2
        cable.C_m(ind_in) = c_in * 2 * pi * R_in * l_in;        % Compartments' membrane capacitance, in uF
        cable.is_node(ind_in) = false;                          % Compartments are internode
        
        cable.z_ind_no_act = [ (1 : (N_myl_comp + (N_in + 1) ) + 1) , ...
                               N_comp + (-(N_myl_comp + (N_in + 1) ) -1 : 0 ) ]; % Disabling the terminals, including 1 internodes and nodes inward from terminal
        cable.Area([1,N_comp]) = cable.Area([1,N_comp]) * 100;      % Compartment's area; cm^2. Increase surface area by 100 to reduce polarization near terminal due to undulation
        cable.C_m([1,N_comp])  = c_n * cable.Area([1,N_comp]);      % Compartment's membrane capacitance, in uF
    
        if modprmtr.antidromic
            z_ind_AP = 1 + N_myl_comp + (1 + N_in) * 15;
        else
            z_ind_AP = N_comp - N_myl_comp - ( 1+ N_in ) * 15;
        end
        z_prop = abs(cable.z(N_comp) - cable.z(1));
        cable.z = cable.z - N_myl_comp * l_in - (min_node_uni * (l_n + L_in) + l_n/2); 	% Shift origin of coordinates to middle
        
    case 3
        cable.z_ind_no_act = [1,N_comp];                             	% Disabling the terminals
        cable.z = cable.z - (min_node_uni * (l_n + L_in) + l_n/2); 	% Shift origin of coordinates to middle
        z_ind_AP = min( [N_comp, max([find(cable.z >= sqrt(3/2)*H*2.5,1,'first'),find(cable.z >= (9.5*l_n+9*L_in),1,'first')]) ] );
        % Distance for AP to pass; minimum of 10 mm or 2.5 times distance where
        % hyperpolarization is strongest (at sqrt(3/2) * H = 1.2 * H): 3*H
        z_prop = min(abs(median(cable.z) - cable.z(z_ind_AP))); % Assume AP initiation somewhere close to center of axon
end


% Axial resistance; sealed ends are reflected by d_phi_e_left/right = 0 at terminals
cable.R_i_left(1) = inf; 
cable.R_i_left(2:end)    = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(1:end-1) = ( cable.R_i(1:end-1) + cable.R_i(2:end))/2;
cable.R_i_right(end) = inf; 

% Biophysics of cable
cable.V_rest = V_rest;
cable.g_in = g_in;

% Nodal ion channels
cable.node.E_Na = E_Na;
cable.node.E_K  = E_K;
cable.node.E_L  = E_L;

cable.node.g_Na = g_Na;
cable.node.g_Nap = g_Nap;
cable.node.g_K  = g_K;
cable.node.g_L  = g_L;

%% E-field at cable

switch modprmtr.mod
    case 1
        E = 1;                      %  1 mV/cm
        Ex = E * sin(Alpha) + Empty_C;      % Tranverse E-field in mV/cm
        Ez = E * cos(Alpha) + Empty_C;      % Axial field in mV/cm
        d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2;
    case 2
       beta = cable.z(ind_bend)/R_bend;
       beta = [-1/4*pi * ones(N_straight_comp_l,1) ; beta ; 1/4*pi * ones(N_straight_comp_r,1)];
       
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
            Ez = Er.* cos(PHI);
        else
            I = 1;                                          % 1 uA; scaling performed by waveform
            r = sqrt( H^2 + cable.z.^2 );                   % Distance to point source, in cm
            Ex = I * H/(4*pi*sigma_e) ./ (r.^3);            % Tranverse E-field in mV/cm
            
            phi_e = I / (4*pi*sigma_e) ./ r;                % Extracellular potenial, in mV            
            Ez = I * cable.z/(4*pi*sigma_e).*r.^(-3);           % Axial field in mV/cm
            % d_phi_e = - Ez(1:end-1).*cable.dz(1:end-1)/2 - Ez(2:end).*cable.dz(2:end)/2  ;
        end
        d_phi_e = diff(phi_e);
end

E = sqrt(Ex.^2 + Ez.^2);
ind_sph = cable.node.cable_ID(cable.node.TP_dim == 2);
Ex(ind_sph) = E(ind_sph);

stimulation.d_phi_e_left =  [ 0 ; d_phi_e ];    % Finite difference in potenial, in mV
stimulation.d_phi_e_right = [ d_phi_e ; 0 ];    % Finite difference in potenial, in mV

stimulation.ER_TP = kron(   (1+1./cable.node.TP_dim) .* abs( Ex(cable.node.cable_ID) ) .*...
                             cable.R(cable.node.cable_ID) , cos(theta)  );



%% Time and stimulation waveform
dt = 2e-3;              % Time step 2 us, in ms

n_bf_start = 5;
t_start = - n_bf_start * dt;                                % Pulse on-set delay, in ms
t_end =   ceil( ( PW  + 1 + z_prop/v)/ dt ) * dt;           % Simulation end time, PW + ~ 1 ms for AP initiation + propagation time
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
solver.h_func = str2func('simulate_cable_RMG');

switch modprmtr.mod
    case {1,2}
        E_rh_TP = 0.2/R_n;	% mV/cm
        E_rh = 5;           % mV/cm
        t_ch = 0.3;         % ms
        amp_init = E_rh / ( 1- 2^ (-PW/t_ch));        % mV/cm
        amp_init_TP = E_rh_TP / ( 1- 2^ (-PW/t_ch));  % mV/cm
                
        if Alpha < pi/2 
            cos_alpha = cos(Alpha);
        else
            cos_alpha = cos(pi/2 *0.95);
        end
        
        amp_init = min( amp_init / cos_alpha  , amp_init_TP) * 0.5 *(1 + randn(1) * 0.02 );
        
        solver.plot_t_intv = [0.05:0.05:solver.t_vec(end)];      % 50 us interval for plotting; 20 points per ms
        
        solver.thresh_find.unit_str = 'V/m';
        solver.thresh_find.unit_amp = 0.1; % 1 mV/cm = 0.1 V/m
        factor = (2)^(1/8);
    case 3
        % set initial search amplitude (based on RMG paper)
        t_ch = 0.1;     % chronaxie from RMG; 102+-8 us
        k_rh = 2000;	% uA/cm^2; k is about 100 uA/mm^2 for 0.1 ms in RMG
        k_PW = k_rh / ( 1 - 2^ (-PW/t_ch));     %uA/cm^2
        I_rh = 2;       %uA; I_0 is 25 uA for 0.1 ms in RMG
        I_0_PW = I_rh / ( 1 - 2^ (-PW/t_ch));   %uA
        
        amp_init = - ( I_0_PW + k_PW * H^2 ); % uA, cathodic current
        amp_init = amp_init * 0.2 *(1 + randn(1) * 0.02 );  % start at reduce amplitude and add some variation
        
        if H <= 0.5
            solver.plot_t_intv = [0.02:0.02:solver.t_vec(end)];             % 20 us interval for plotting; 50 points per ms
        else
            solver.plot_t_intv = [0.05:0.05:0.9,1:0.1:solver.t_vec(end)];	% 100 us interval for plotting; 10 points per ms
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