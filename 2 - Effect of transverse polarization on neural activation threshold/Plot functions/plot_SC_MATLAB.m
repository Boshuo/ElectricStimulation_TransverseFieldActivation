figure_format;
figure;
set(gcf,'Position',[00 00 1500 600],'Color',format_figure.Color);

amp_log_max = 6.5;        % Upper bound of threshold amplitudes, 10^amp_log_max, mV/cm
amp_log_min = 2.5;        % Lower bound of threshold amplitudes, 10^amp_log_min, mV/cm

PW_sep = 15e-2;             % ms, 
% separation between two regions of pulse duration 
% long PW region PW>=PW_sep: Lapicque-like behaviour
% short PW region PW<=PW_sep : log-log slope > -1

opts = optimset('MaxFunEvals',1500,'MaxIter',1500,'TolF',1e-12,'TolX',1e-6);

%% Axon
% Parameters reloaded so that script can be run without being called by
% compile_NUERON_data.m
model_name = 'UF_Axon_HH';

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

threshold_axon = compiled_results.th_MCE;       % convert from mV/mm to mV/cm; 1 mV/cm = 10 mV/mm

R_vec_axon =  compiled_results.R(:,1);
PW_vec = compiled_results.PW(1,:);

ind_PW_long = find(PW_vec >= PW_sep,1,'first')-1 : length(PW_vec);
ind_PW_short = 1: find(PW_vec <= PW_sep,1,'last')+ 1;

% long PW region
t_ch_axon = zeros(size(R_vec_axon));    % chronaxies
rh_axon = zeros(size(R_vec_axon));      % rheobase
% short PW region
slope_axon = zeros(size(R_vec_axon));   % log-log slope of fitted line
interc_axon = zeros(size(R_vec_axon));  % intercept of fitted line

color = parula(length(R_vec_axon));
legend_text = cell(size(R_vec_axon));

h_ax(1) = axes('position',[0.075,0.1,0.25,0.70]);
hold on;box on;
axis(log10([1e-3,1e1,10^amp_log_min,10^amp_log_max]));


for ii = 1 : length(R_vec_axon)
    plot(log10(PW_vec),log10(abs(threshold_axon(ii,:))/1e1),'-','LineWidth',2,'Color',color(ii,:));
    
    para_init = [PW_vec(end),abs(threshold_axon(ii,end))]; % t_ch, I_rh
    [para,Jmin,ExitFlag] = fminsearch(@(x) JError_Lapicque(x,PW_vec(ind_PW_long),abs(threshold_axon(ii,ind_PW_long))), para_init,opts);
    t_ch_axon(ii) = para(1);
    rh_axon(ii) = para(2);
    
    stats = regstats(log10(abs(threshold_axon(ii,ind_PW_short))),log10(PW_vec(ind_PW_short)),'linear',{'beta','rsquare'});
    r_square = stats.rsquare;
    %  disp(r_square);
    slope_axon(ii) = stats.beta(2);
    interc_axon(ii) = stats.beta(1);
    % disp(stats.beta)
    
    if R_vec_axon(ii) >= 1e-4 
        legend_text{ii} = ['$$R =',num2str(R_vec_axon(ii)/1e-4,'%1.1f'),' \: \rm{ \mu m}$$ '];
    else
        legend_text{ii} = ['$$R =',num2str(R_vec_axon(ii)/1e-4,'%1.2f'),' \: \rm{ \mu m}$$ '];
    end 

end

for ii = 1 : length(R_vec_axon)
    E_fit = rh_axon(ii) ./ (1 - 2.^ ( - PW_vec(ind_PW_long)./t_ch_axon(ii)));
    plot(log10(PW_vec(ind_PW_long)),log10(E_fit/1e1),'--','LineWidth',1.5,'Color','k');
    plot(log10(t_ch_axon(ii)),log10(rh_axon(ii)*2/1e1),'k*');
    plot(log10(PW_vec(ind_PW_short)),log10(PW_vec(ind_PW_short))*slope_axon(ii) + interc_axon(ii) -1,':','LineWidth',2,'Color','k');
end

h_legend(1) = legend(legend_text,'Location','Northeast');
h_xlabel(1) = xlabel('Pulse duration $$PW \: \rm{(ms)}$$','Interpreter','latex');
h_ylabel(1) = ylabel('Threshold $$ {E''} \: \rm{(mV/mm)}$$','Interpreter','latex');
h_title(1) = title({['Axon: log-log slope: $$',num2str(slope_axon(end),'%1.3f'),'$$'],...
                    ['$$t_{\rm{ch}} =',num2str(t_ch_axon(end)/1e-3,'%1.1f'),' \: \rm{ \mu s}$$, ',...
                    '$${E''}_{\rm{rh}}\cdot R =',num2str(rh_axon(end)*R_vec_axon(end),'%1.3f'),' \: \rm{mV}$$']},'Interpreter','latex');

%% Soma
filename = fullfile('Processed data and figures','UF_Soma_HH_compiled_result_NEURON.mat');
load(filename,'threshold','PW_vec','R_vec');

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

threshold_soma = compiled_results.th_MCE;       

R_vec_soma =  compiled_results.R(:,1);
PW_vec = compiled_results.PW(1,:);

% long PW region
t_ch_soma = zeros(size(R_vec_soma));    % chronaxies
rh_soma = zeros(size(R_vec_soma));      % rheobase
% short PW region
slope_soma = zeros(size(R_vec_soma));   % log-log slope of fitted line
interc_soma = zeros(size(R_vec_soma));  % intercept of fitted line

h_ax(2) = axes('position',[0.375,0.1,0.25,0.70]);
hold on;box on;
axis(log10([1e-3,1e1,10^amp_log_min,10^amp_log_max]));

for ii = 1 : length(R_vec_soma)
    plot(log10(PW_vec(:)),log10(abs(threshold_soma(ii,:))/1e1),'-','LineWidth',2,'Color',color(ii,:));
    
    para_init = [PW_vec(end),abs(threshold_soma(ii,end))]; % t_ch, I_rh
    [para,Jmin,ExitFlag] = fminsearch(@(x) JError_Lapicque(x,PW_vec(ind_PW_long),abs(threshold_soma(ii,ind_PW_long))), para_init,opts);
    t_ch_soma(ii) = para(1);
    rh_soma(ii) = para(2);
    
    stats = regstats(log10(abs(threshold_soma(ii,ind_PW_short))),log10(PW_vec(ind_PW_short)),'linear',{'beta','rsquare'});
    r_square = stats.rsquare;
    % disp(r_square);
    slope_soma(ii) = stats.beta(2);
    interc_soma(ii) = stats.beta(1);
    % disp(stats.beta)
    if R_vec_soma(ii) >= 1e-4
        legend_text{ii} = ['$$R =',num2str(R_vec_soma(ii)/1e-4,'%1.1f'),' \: \rm{ \mu m}$$ '];
    else
        legend_text{ii} = ['$$R =',num2str(R_vec_soma(ii)/1e-4,'%1.2f'),' \: \rm{ \mu m}$$ '];
    end

end
 
for ii = 1 : length(R_vec_soma)
    E_fit = rh_soma(ii) ./ (1 - 2.^ ( - PW_vec(ind_PW_long)./t_ch_soma(ii)));
    plot(log10(PW_vec(ind_PW_long)),log10(E_fit/1e1),'--','LineWidth',1.5,'Color','k');
    plot(log10(t_ch_soma(ii)),log10(rh_soma(ii)*2/1e1),'k*');
    plot(log10(PW_vec(ind_PW_short)),log10(PW_vec(ind_PW_short))*slope_soma(ii) + interc_soma(ii) -1,':','LineWidth',2,'Color','k');
end

h_legend(2) = legend(legend_text,'Location','Northeast');
h_xlabel(2) = xlabel('Pulse duration $$PW \: \rm{(ms)}$$','Interpreter','latex');
h_title(2) = title({['Soma: log-log slope: $$',num2str(slope_soma(end),'%1.3f'),'$$'],...
                    ['$$t_{\rm{ch}} =',num2str(t_ch_soma(end)/1e-3,'%1.1f'),' \: \rm{ \mu s}$$, ',...
                    '$${E''}_{\rm{rh}}\cdot R =',num2str(rh_soma(end)*R_vec_soma(end),'%1.3f'),' \: \rm{mV}$$']},'Interpreter','latex');
%%

xticks = -3:1;
xticklabels = cell(length(xticks),1);
for ii = 1 : length(xticklabels)
    xticklabels{ii} = ['$$10^{',num2str(xticks(ii),'%d'),'}$$'];
end

yticks = ceil(amp_log_min):floor(amp_log_max);
set(h_ax(1:2),'YTick',yticks);
yticklabels = cell(length(yticks),1);
for ii = 1 : length(yticklabels)
    yticklabels{ii} = ['$$10^{',num2str(yticks(ii),'%d'),'}$$'];
end
set(h_ax(1:2),'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels);

%%
h_ax(3) = axes('position',[0.7,0.1,0.25,0.70]);
hold on;box on;

h_line(1) = plot(log10(PW_vec),log10(abs(threshold_axon(end,:)) * R_vec_axon(end)),'-sk','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','none');
h_line(2) = plot(log10(PW_vec),log10(abs(threshold_soma(end,:)) * R_vec_soma(end)),'-ok','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','none');
legend_text = {'Axon','Soma'};

for ii = 1 : length(R_vec_axon)-1
    plot(log10(PW_vec(:)),log10(abs(threshold_axon(ii,:)) * R_vec_axon(ii)),'-sk','LineWidth',0.5,'MarkerSize',6,'MarkerFaceColor','none');    
end
for ii = 1 : length(R_vec_soma)-1
    plot(log10(PW_vec(:)),log10(abs(threshold_soma(ii,:)) * R_vec_soma(ii)),'-ok','LineWidth',0.5,'MarkerSize',6,'MarkerFaceColor','none');
end

axis(log10([1e-3,1e1,10^0,10^4]));

set(gca,'YTick',0:4);
yticks = 0:4;
yticklabels = cell(length(yticks),1);
for ii = 1 : length(yticklabels)
    yticklabels{ii} = ['$$10^{',num2str(yticks(ii),'%d'),'}$$'];
end
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels);

h_legend(3) = legend(h_line,legend_text,'Location','Northeast');
h_xlabel(3) = xlabel('Pulse duration $$PW \: \rm{(ms)}$$','Interpreter','latex');
h_xlabel(3) = ylabel('Normalized threshold $$ {E}''_{\rm{x}} \cdot R \: \rm{(mV)}$$','Interpreter','latex');
h_title(3) = title({'Normalized threshold',''},'Interpreter','latex');

%%
set(h_ylabel, format_axis_label);
set(h_ylabel, format_axis_label);
set(h_title, format_title);
set(h_ax, format_axis);
set(h_legend,'Box','off','Interpreter','latex','FontSize',12);
`

saveas(gcf,'UF_HH_MATLAB.fig');
[imind,cm] = rgb2ind(frame2im(getframe(gcf)),256);
imwrite(    imind,cm,'UF_HH_MATLAB.tif',...
            'tif','WriteMode','overwrite', 'Resolution',300);
close(gcf);