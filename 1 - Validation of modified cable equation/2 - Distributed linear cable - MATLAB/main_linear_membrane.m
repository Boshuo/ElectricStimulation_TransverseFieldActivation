clearvars;
figure_format;
%%
if ~exist('Results','dir')
    mkdir('Results');
end
mkdir('Figures');

electrode_str = {'PE','DE'};

for electrode_type = 1 : 2
    [results,PW] = ES_A_Lin(electrode_type);
    
    filename = fullfile('Results',[electrode_str{electrode_type},'_A_Lin.mat']);
    save(filename,'results','PW');
    
    if electrode_type == 1
        plot_trans = 1;
        plot_Schnabel2001;        
    else 
        plot_trans = 0;
    end
    plot_Neu2016;
end