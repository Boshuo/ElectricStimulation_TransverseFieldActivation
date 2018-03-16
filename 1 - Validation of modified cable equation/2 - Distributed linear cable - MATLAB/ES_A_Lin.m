function [results,PW] = ES_A_Lin(electrode_type)
% electrode_type: 1, point source; 2 disk electrode

%% Parameters
%%% Cellular
sigma_e = 10;               % Extracellular conductivity, in mS/cm; 
sigma_i = 28.2 ;            % Intracellular conductivity, in mS/cm; 

R = 3 * 1e-4;               % Axon radius 3 um, in cm;

c_m = 1;                    % Specific membrane capacitance, in uF/cm^2
r_m = 1.5;                  % Specific membrane resistance, in kOhm*cm^2 
g_m = 1 / r_m;              % Specific membrane conductance, in mS/cm^2
tau = c_m * r_m;            % Membrane time constant, in ms

E_L = 0;                    % Leakage reversal potential (reduced), in mV

%%% Electrode
I = -2e3;                   % Electrode current amplitude, in uA
H = [4.5e-4,6.8e-4, kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),10];     
% Axon-electrode distance in cm, 4.5 um to 10 cm; 6 point per decade
PW = [kron([1e-4,1e-3,1e-2,1e-1,1e0,1e1],[1,1.4,2.1,3.1,4.5,6.8]),1e2];
% PW in ms, 0.1 us to 100 ms; 
% 6 point per decade: 
% 1         1.4       2.1       3.1       4.5       6.8
% 1.0000    1.4678    2.1544    3.1623    4.6416    6.8129

%%% Other
d_lambda = 0.1;             % d_lambda rule of NEURON
freq = 100e-3;              % Frequency for calculating d_lambda, 100 Hz in kHz
lambda_f = sqrt(R *sigma_i / (2*pi*freq*c_m) /2); % Length constant at 100 Hz, in cm
% For R = 3um: lambda100 = 820.5 um; d_lambda = 0.1 -> compartment length = 82.1 um 
% NOTE: NEURON has lambda_f sqrt(2) times larger
lambda_DC = sqrt(R * sigma_i / g_m /2);
% DC length constant, for R = 3 um, lambda_DC = 0.79 mm
axon_min_length = lambda_DC * 5;
% Minimum axon length at least 5 DC length constant. for R = 3 um, axon_min_length = 4.0 mm
    
%%% Simulation time
dt_default = 1e-3;                          % default time step in ms, 1 us
n_bf_start = 3;                             % number of time steps before pulse start
t_start = - n_bf_start * dt_default;        % Simulation time start in ms, with pulse onset at 0
t_end = max(PW);                            % Simulation end time

t_vec = t_start:dt_default:t_end;           % Time vector of simulation, in ms
t_vec = unique(round([t_vec,PW]/1e-5))*1e-5;% Add time points when pulses end, the scaling of values to intergers is to allow "unique" work properly 
dt = [0,diff(t_vec)];                       % Actual time steps, in ms

%% Pre-allocation
emptymat = zeros(size(PW));
emptycell = cell(size(H));
results = struct(	'H',num2cell(H),...     % Axon-electrode distance
                  	'ExR',emptycell,...     % Maximum of E * R, at center of axon
                    't_ss',emptycell,...    % Time to reach steady state polarization
                    'phi_m_bar',emptymat... % Maximum average membrane potential
                    );

%% Loops for all parameters    
for ii = 1:length(H)
    tic_H = tic;
    disp(['Axon-electrode distance: ',num2str(H(ii)/1e-1,'%2.5f'),' mm. ']);
    %%% Axon length and compartment length
    if (10*H(ii) >= axon_min_length)                % H > 400 um (H >= 450 um)
        % If 10*H exceeds minimum axon length, axon length set to 10*H
        axon_length = 10 * H(ii);
    else                                            % H < 400 um (H <= 310 um)
        axon_length = axon_min_length;
        % Otherwise axon length is set to 5 DC length constant (rounded)
    end
    if ( H(ii)/10 <= lambda_f * d_lambda )          % H < 821 um (H <= 680 um)
        dz = H(ii) / 10;                            % Compartment length
        % If H/10 is less than inteval specified by d_lambda rule, then set interval to distance /10
    else                                            % H > 821 um (H >= 1 mm)
        dz = lambda_f * d_lambda;
        % Otherwise interval is determined by d_lambda rule
    end
    axon_length = dz * ceil(axon_length/dz);
    
    %%% Geometric parameters and E-field
    z = linspace( -axon_length, axon_length, 2*round(axon_length/dz)+1 )';
    % Axial coordinates of axon, in cm
    display(['Number of compartments: ',num2str(length(z))]);
    switch electrode_type
        case 1
            r = sqrt(H(ii)^2+z.^2);                     % Distance between compartment center to point source, in cm
            phi_e = I / (4*pi*sigma_e) ./ r;            % Extracellular potenial at compartment center, in mV
            delta_phi_e = [diff(phi_e(1:2)); diff(phi_e,2); -diff(phi_e(end-1:end))];
            % Differential of longitudinal field, first order at terminals,
            % otherwise second order
            Ex = I * H(ii)/(4*pi*sigma_e) ./ (r.^3);	% Extracellular tranverse E-field in mV/cm
            % Ez = I * z/(4*pi*sigma_e).*r.^(-3);      	% Extracellular axial field in mV/cm
        case 2
            x = 0;
            R_disk = 100e-4;        % 200 um diameter disk
            [PHI, RR] = cart2pol(z,x);  % cylindrical coordinates with azimuthal angle measured from the z axis
            
            H_hat = H(ii)/R_disk;
            R_hat = RR/R_disk;
            
            sqrt_minus = sqrt( (R_hat-1).^2 + H_hat^2  );
            sqrt_plus  = sqrt( (R_hat+1).^2 + H_hat^2  );
            sqrt_big   = sqrt(  (sqrt_minus + sqrt_plus).^2 - 4 );
            
            Imp_disk = 1 / (4*R_disk*sigma_e);	% Resistance in kOhm
            V_0 = I * Imp_disk;                 % Primary voltage in mV
            phi_e = V_0 * 2/pi*asin( 2 ./ (sqrt_plus + sqrt_minus) );  % Extracellular potenial at compartment center, in mV
            delta_phi_e = [diff(phi_e(1:2)); diff(phi_e,2); -diff(phi_e(end-1:end))];
            % Differential of longitudinal field, first order at terminals,
            % otherwise second order
            Er = V_0 *4/ (pi*R_disk) * ( (R_hat-1)./sqrt_minus + (R_hat+1)./sqrt_plus ) ./ ( (sqrt_plus + sqrt_minus) .* sqrt_big );
            Ey = V_0 *4/ (pi*R_disk) * ( (H_hat  )./sqrt_minus + (H_hat  )./sqrt_plus ) ./ ( (sqrt_plus + sqrt_minus) .* sqrt_big );
            
            Ex = sqrt((Er.* sin(PHI)).^2 + Ey.^2);      % Extracellular tranverse E-field in mV/cm
    end
    
    weight = 2 * ones(size(z));                 % Weights
    weight(1) = 1; 
    weight(end) = 1;
    
    results(ii).ExR =  max(abs(Ex)) * R;             % in mV
             
    %%% Neuronal parameters     
    phi_m_bar = zeros(size(z));     % Average membrane potential, mV
    A = 2 * pi * R * dz;            % Compartment area; cm^2
    R_i = dz / (sigma_i*pi*R^2);   	% Axial resistance between compartments, in kOhm
    Lambda2 = ones(size(z)) ./ (R_i * g_m * A);  % Local length constant, same along axon for linear membrane
    
    %%% Simulations
    tic_PW = tic;
    for t_ind= 2 : length(t_vec)   % Backward Euler steps
        
        % Generate tridiagonal matrix and right-hand-side for backward Euler
        F = tau .* phi_m_bar + dt(t_ind) * ( E_L + Lambda2 .* delta_phi_e * (t_vec(t_ind)>0) );
        A = tau + dt(t_ind) + dt(t_ind) * Lambda2 .* weight;
        B = -Lambda2 * dt(t_ind);
        
        phi_m_bar_old = phi_m_bar;          %% Save membrane potential of previous time
        phi_m_bar = tridiag(A, B, B, F);    %% Backward Euler step
        
        ind_PW = find(abs (t_vec(t_ind) - PW) < 1e-6);      %% Record results if the time step is equal to one of the pulse durations
        if ~isempty(ind_PW)
            results(ii).phi_m_bar(ind_PW) = max(phi_m_bar);
            t_PW = toc(tic_PW);tic_PW = tic;
            disp(['PW: ', num2str(PW(ind_PW),'%1.4f'),' ms. Computation time: ', num2str(t_PW,'%2.5f'),' s.']);
        end
        
        if max(abs (phi_m_bar ./ phi_m_bar_old -1)) <= 1e-4         % Skip remaining simulation time if steady-state reached
            results(ii).t_ss = t_vec(t_ind);               % Steady state time
            ind_PW = find(t_vec(t_ind)<PW);
            for jj = ind_PW
                results(ii).phi_m_bar(jj) = max(phi_m_bar);
            end
            disp(['Steady state reached at t=',num2str(t_vec(t_ind),'%2.3f'),' ms.']);
            break
        end
    end
    
    t_H = toc(tic_H);
    disp(['Total computation time: ', num2str(t_H,'%2.5f'),' s.']);
    disp(' ');
    
end

end