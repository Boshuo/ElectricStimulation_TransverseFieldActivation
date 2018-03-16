model = {   'UF_Axon_HH',   'UF_Soma_HH',...
            'UF_TrmB_HH',   'UF_TrmB_RMG',  'UF_TrmS_RMG',...
            'UF_BndO_HH',   'UF_BndA_HH',   'UF_BndO_RMG',	'UF_BndA_RMG',...
            'PE_A_HH',      'DE_A_HH',      'PE_A_RMG',     'DE_A_RMG'};
        

for ii = 1 : length(model)
    model_name = model{ii};
    compile(model_name);
end


function compile(model_name)
PW_vec = [ kron([1e-3,1e-2,1e-1,1e0],[1,1.4,2.1,3.1,4.5,6.8]),1e1];
% PW in ms, 1 us to 10 ms; >10 ms rheobase (steady state for linear
% membrane); 25 points (4 decades)

if strncmp(model_name,'UF_',3)                  % Uniform Field
    if strncmp(model_name(4:end),'Axon',4)      % Axon
        mod = 0;
        R_vec = [0.25,0.5,1,2,4,8]*1e-4;        % Compartment radius, in cm;
        [RR,PPWW] =  ndgrid(R_vec,PW_vec);
    elseif strncmp(model_name(4:end),'Soma',4)  % Soma
        mod = 0;
        R_vec = [     0.5,1,2,4,8,16]*1e-4; 	% Compartment radius, in cm;
        [RR,PPWW] =  ndgrid(R_vec,PW_vec);
    elseif strncmp(model_name(4:end),'Trm',3)       % Axon TeRMinal
        mod = 1;
        % Parameters: 600 = 24 * 25
        alpha_vec = [0:5:75,78:3:87,88,89,89.5,90]*pi/180;  
        % various interval, 24
         [ALPHA,PPWW] = ndgrid(alpha_vec,PW_vec);    
    elseif strncmp(model_name(4:end),'Bnd',3)  % Axonal BeND
        mod = 2;  
        % Parameters: 600 =  24 * 25;
        alpha_vec = [0:15:60,70:10:90,95:5:105,108:115,117,120:5:135]*pi/180;
        % various, 24
         [ALPHA,PPWW] = ndgrid(alpha_vec,PW_vec);    
    end   
else                                            % Electrode stimulation
    mod = 3;
    % Parameters: 475 =  19 * 25;
    H_vec = [ [3.1,4.5,6.8]*1e-3, kron([1e-2,1e-1],[1,1.4,2.1,3.1,4.5,6.8]),[1,1.4,2.1,3.1] ];     
    % Electrode-axon distance in cm, 31 um to 3 cm; 19 points (3 decades)
    [HH,PPWW] = ndgrid(H_vec,PW_vec);    
end

NaN_matrix = NaN(size(PPWW));

switch mod
    case 0
        compiled_results = struct(  'model',model_name,...
                                    'R',RR,'PW',PPWW,...
                                    'parameter_id',NaN_matrix,...
                                    'th_MCE',NaN_matrix);
    case {1,2}
        compiled_results = struct(  'model',model_name,...
                                    'ALPHA',ALPHA,'PPWW',PPWW,...
                                    'parameter_id',NaN_matrix,...
                                    'th_CE',NaN_matrix,...
                                    'th_MCE',NaN_matrix,...
                                    'th_per_diff_MCE',NaN_matrix);
    case 3
        compiled_results = struct(  'model',model_name,...
                                    'HH',HH,'PPWW',PPWW,...
                                    'parameter_id',NaN_matrix,...
                                    'th_CE',NaN_matrix,...
                                    'th_MCE',NaN_matrix,...
                                    'th_per_diff_MCE',NaN_matrix);
end

if exist([model_name,'/Results'],'dir') ==  0
    error('LoadConvertDataFromCluster:NoDataFromClusterParameters','Data folder does not exist for given parameters.');
end

for ii = 1: numel(PPWW)
    filename = fullfile(model_name,'Results',['result_',num2str(ii),'.mat']);
    if exist(filename,'file') > 0
        load(filename,'mod_prmtr', 'results');
        
        parameter_id = mod_prmtr.id;
        compiled_results.parameter_id(parameter_id)	= parameter_id;
        if mod
            compiled_results.th_CE(parameter_id) = results.th_CE;
            compiled_results.th_per_diff_MCE(parameter_id) = results.th_per_diff_MCE;
        end
        compiled_results.th_MCE(parameter_id) = results.th_MCE;
        
    else
    disp(['File does not exist: ',filename]);
    end
end

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
save(filename,'compiled_results');
disp(['Compiled ', model_name,'.']);
end
