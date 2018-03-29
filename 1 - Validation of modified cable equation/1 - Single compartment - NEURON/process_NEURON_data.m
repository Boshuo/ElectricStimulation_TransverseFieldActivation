clearvars;clc;
addpath('MATLAB functions');        % Path of shared functions and data

if ~exist('Processed data and figures','dir')
    mkdir('Processed data and figures');
end
%% Set folder and R_vec according to radii used in the NEURON
folder = 'NEURON results';

PW_vec = [kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),10]; % ms, 25 pulse duration

%% Axon
R_vec = 0.25 * 2.^(0:6-1);       % Compartment radius, in um; 

threshold.E = zeros(length(PW_vec),length(R_vec));
threshold.ER = zeros(length(PW_vec),length(R_vec));

for ii = 1:numel(threshold.E)
    filename  = sprintf('%s\\axon\\%d.txt',folder,ii);
    FID = fopen(filename,'r');
    [str] = fgetl(FID);
    while ~strncmp(str,'Threshold search',16)
        [str] = fgetl(FID);
    end
    data = fscanf(FID,'\tE-field: %f mV/mm\n\tE*R: %f mV\n');
    threshold.E(ii) = data(1);
    threshold.ER(ii) = data(2);
    fclose(FID);
end

filename = fullfile('Processed data and figures','UF_Axon_HH_compiled_result_NEURON.mat');
save(filename,'threshold','R_vec','PW_vec');

%% Soma
R_vec = 0.5 * 2.^(0:6-1);       % Compartment radius, in um; 

threshold.E = zeros(length(PW_vec),length(R_vec));
threshold.ER = zeros(length(PW_vec),length(R_vec));

for ii = 1:numel(threshold.E)
    filename  = sprintf('%s\\soma\\%d.txt',folder,ii);
    FID = fopen(filename,'r');
    [str] = fgetl(FID);
    while ~strncmp(str,'Threshold search',16)
        [str] = fgetl(FID);
    end
    data = fscanf(FID,'\tE-field: %f mV/mm\n\tE*R: %f mV\n');
    threshold.E(ii) = data(1);
    threshold.ER(ii) = data(2);
    fclose(FID);
end

filename = fullfile('Processed data and figures','UF_Soma_HH_compiled_result_NEURON.mat');
save(filename,'threshold','R_vec','PW_vec');

plot_UF_SC_HH_NEURON;                 % plot NEURON results
plot_UF_SC_HH_NEURON_and_MATLAB;      % plot NEURON and MATLAB results

rmpath('MATLAB functions');
