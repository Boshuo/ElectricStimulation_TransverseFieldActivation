figure_format;

h_f = figure;
set(h_f,'Position',[00 00 720 600],'Color',format_figure.Color);
box on; hold on;

%% Axon
model_name = 'UF_Axon_HH';

filename = fullfile('..','..','2 - Effect of transverse polarization on neural activation threshold',model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

threshold_axon_M = compiled_results.th_MCE;       % convert from mV/mm to mV/cm; 1 mV/cm = 10 mV/mm

PW_axon = compiled_results.PW;
R_axon  = compiled_results.R;    %R_vec = [0.25,0.5,1,2,3,4,8,16]*1e-4;   

R_vec_axon_M =  R_axon(:,1);       % Compartment radius, in cm;
PW_vec_axon_M = PW_axon(1,:);

filename = fullfile('Processed data and figures','UF_Axon_HH_compiled_result_NEURON.mat');
load(filename,'threshold','PW_vec','R_vec');

R_vec_axon_N =  R_vec*1e-4;       % Compartment radius, in cm;
PW_vec_axon_N = PW_vec;

ind_R_N = length(R_vec_axon_N);
ind_R_M = length(R_vec_axon_M);

threshold_axon_N = threshold.E'*10;
threshold_axon_N(~threshold_axon_N) = NaN;

%% Soma
model_name = 'UF_Soma_HH';

filename = fullfile('..','..','2 - Effect of transverse polarization on neural activation threshold',model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

threshold_soma_M = compiled_results.th_MCE;       % convert from mV/mm to mV/cm; 1 mV/cm = 10 mV/mm

PW_soma = compiled_results.PW;
R_soma  = compiled_results.R;   %R_vec = [1,2,3,4,8,16,32]*1e-4;  
R_vec_soma_M =  R_soma(:,1);       % Compartment radius, in cm;
PW_vec_soma_M = PW_soma(1,:);


filename = fullfile('Processed data and figures','UF_Soma_HH_compiled_result_NEURON.mat');
load(filename,'threshold','PW_vec','R_vec');

R_vec_soma_N =  R_vec*1e-4;       % Compartment radius, in cm;
PW_vec_soma_N = PW_vec;

threshold_soma_N = threshold.E'*10;   % from V/m to mV/cm
threshold_soma_N(~threshold_soma_N) = NaN;

%%
plot(log10(PW_vec_axon_M),log10(abs(threshold_axon_M(ind_R_M,:))*R_vec_axon_M(ind_R_M)),'-s','Color',[1,1,1]*0.6,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[1,1,1]*0.6,'MarkerEdgeColor','none');
plot(log10(PW_vec_axon_N),log10(abs(threshold_axon_N(ind_R_N,:))*R_vec_axon_N(ind_R_N)),'-s','Color',[1,1,1]*0.0,'LineWidth',1.5,'MarkerSize',6, 'MarkerFaceColor','w','MarkerEdgeColor',[1,1,1]*0.0);

plot(log10(PW_vec_soma_M),log10(abs(threshold_soma_M(ind_R_M,:))*R_vec_soma_M(ind_R_M)),'-o','Color',[1,1,1]*0.6,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[1,1,1]*0.6,'MarkerEdgeColor','none');
plot(log10(PW_vec_soma_N),log10(abs(threshold_soma_N(ind_R_N,:))*R_vec_soma_N(ind_R_N)),'-o','Color',[1,1,1]*0.0,'LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',[1,1,1]*0.0);

%%
axis equal;

axis(log10([0.8e-3,1.25e1,10^0.5,10^3.5]));

legend_text = {'Axon - MATLAB','Axon - NEURON','Soma - MATLAB','Soma - NEURON'};
h_legend = legend(legend_text,'Location','Northeast');
set(h_legend,'Box','off','Interpreter','latex','FontSize',16);

xticks = log10(kron(10.^(-3:1),(1:9)));
xticklabels = cell(length(xticks),1);
for ii = 1 : length(xticklabels)
    if round(xticks(ii)) == xticks(ii)
        xticklabels{ii}= ['$$10^{',num2str(xticks(ii),'%d'),'}$$'];
    end
end

yticks = log10(kron(10.^(0:4),(1:9)));
yticklabels = cell(length(yticks),1);
for ii = 1 : length(yticklabels)
     if round(yticks(ii)) == yticks(ii)
        yticklabels{ii}= ['$$10^{',num2str(yticks(ii),'%d'),'}$$'];
     end
end
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels);

xlabel('Pulse width $$PW \: \rm{(ms)}$$','Interpreter','latex');
ylabel('Normalized threshold $$ E''_{\mathrm{x}}  R \: \mathrm{(mV)}$$','Interpreter','latex');
title({'Strength{\textendash}duration curve for transverse activation'},'Interpreter','latex');

set([get(gca,'XLabel'),get(gca,'YLabel')], format_axis_label);
set(get(gca,'Title'), format_title);
set(gca,format_axis);
set(gca,'XMinorTick','off','YMinorTick','off');

saveas(gca,fullfile('Processed data and figures','UF_SC_HH_combined.fig'));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Processed data and figures','UF_SC_HH_combined.tif'),...
            'tif','WriteMode','overwrite', 'Resolution',300);
close(h_f);