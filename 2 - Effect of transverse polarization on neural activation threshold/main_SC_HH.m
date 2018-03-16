function results = main_SC_HH(	mod_prmtr, out_ctrl )
%
% mod_prmtr: structure specifying model parameter:
% model_name                Model: 'UF_Axon', 'UF_Soma'
% id                        Parameter ID of test case (integer)
%
% out_ctrl: structure specifying output control:
% if_save_data   	Whether to save results in a .mat file for each simulation (logical or 0/1)
% if_write_log  	Whether to write threshold finding process in a .txt log (logical or 0/1)
% if_plot         	Whether to plot threshold finding process (logical or 0/1)

%% Logistics
addpath('Shared functions and data');        % Path of shared functions and data

folder_name = [ mod_prmtr.model_name,'_HH'];               % Folder for model, add membrane type to folder name          
create_folders(out_ctrl,folder_name);                       % Create subfolders for data, logs, and figures

out_ctrl.log_fid = 0;                                       % Defaults is to display in MATLAB command window
if out_ctrl.if_write_log                                    % Write in .txt file
    logfilename = fullfile(folder_name,'Logs',['log_',num2str(mod_prmtr.id),'.txt']);	% Log filename
    write_fun( out_ctrl.log_fid,{' ',['Simulation process will be written in: ',logfilename]});
    out_ctrl.log_fid = fopen(logfilename,'w');              % Open log file
end

%% Set parameters for simulation loop
MCE = 1;                                                    % 1 using conventional cable equation (CE), 1 using modified CE
CE_str = {'MCE'};                                             % Description of CE used

results = struct(   'th_MCE', NaN   );  

%%  Main loop for simulations
t_main = tic;
[solver, stimulation, cable] = specify_model_SC_HH( mod_prmtr );	% Specify parameters for solver, stimulation, and cable

t_th = tic;
    
solver.is_MCE =  MCE;                                               % Whether to use MCE
ii = 1;

if ii == 1
    write_fun( out_ctrl.log_fid, solver.txt.log_txt);              	% Output model specification related text
end

write_fun(out_ctrl.log_fid, {'-----------------------------------------------------------',...
    ['Solver: ',CE_str{ii}],' '});

if out_ctrl.if_plot         % Set up figure
    out_ctrl.h_fig = figure('Position',[00 00 1400 800],'Color',[1,1,1]);
    axes('position',[0.05,0.95,0.9,0.00]);box off; axis off;
    title(  [solver.txt.fig_title, CE_str{ii}], 'Interpreter','latex','FontSize',14);
end

results.(['th_',CE_str{ii}]) = threshold_finding( solver, stimulation, cable, out_ctrl );

if out_ctrl.if_plot                                         % Save and close figures
    figure_filename = fullfile(folder_name,'Figures',['Fig_',num2str(mod_prmtr.id),'_',CE_str{ii},'.fig']);
    saveas(out_ctrl.h_fig,figure_filename,'fig');
    close(out_ctrl.h_fig);
end

write_fun(out_ctrl.log_fid, {' ',sprintf('Run time for search: %3.2f min.',toc(t_th)/60),' '});

%% Saving results and closing files
write_fun( out_ctrl.log_fid,    {'-----------------------------------------------------------',...
                                 sprintf('Total run time: %3.3f min.', toc(t_main)/60)});    % Output run time 
if out_ctrl.if_save_data
    filename = fullfile(folder_name,'Results',['result_',num2str(mod_prmtr.id),'.mat']);
    save(filename,'mod_prmtr', 'results');
    write_fun(out_ctrl.log_fid, {sprintf('Results saved in %s.', filename)});
end
if out_ctrl.if_write_log
    out_ctrl.log_fid = fclose(out_ctrl.log_fid);
end
rmpath('Shared functions and data');
end
