clearvars;close all;
addpath('Plot functions');
figure_format;
rmpath('Plot functions');

h_f = figure;
set(h_f,'Position',[0,0,1500,800],format_figure);

xtick = 0 : 5 : 90;
xticklabel = cell(length(xtick),1);
for ii = 1 : length(xticklabel)
    if abs(round(xtick(ii)/15)-xtick(ii)/15) < 1e-2
        xticklabel{ii} = ['$$',num2str((xtick(ii)),'%d'),'^{\circ}$$'];
    else
        xticklabel{ii} = ' ';
    end
end

E_log_max = 5.5;    E_log_cnt_max = fix(E_log_max);
E_log_min = 0.5;    E_log_cnt_min = fix(E_log_min);
E_log_lvl = E_log_min:0.01:E_log_max;
E_log_range = [E_log_min,E_log_max];

max_diff = -100;
per_diff_lvl = [max_diff:5:-30,-25:2.5:0];

%
axes('Position',[0.08,0.55,0.88,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{E-field threshold}','modified CE'},'Interpreter','latex','Color','k','FontSize',20);

caxis(gca,E_log_range);
colormap(gca,cm1);

ticks = E_log_cnt_min : E_log_cnt_max;
ticklabels = cell(size(ticks));
for ii = 1 : length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',ticks(ii));
end

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.915,0.55,0.015,0.4],'Ticks',ticks,'TickLabels',ticklabels);
cb_title.String = '$$\mathrm{V/m}$$';
set(hcb.Title,cb_title);

%
axes('Position',[0.08,0.1,0.88,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{Percentage change}','vs. conventional CE'},'Interpreter','latex','Color','k','FontSize', 20);

caxis(gca,[max_diff,0]);
colormap(gca,cm2);

ticks = max_diff : 25 : 0;
ticklabels = cell(size(ticks));

for ii = 1 : length(ticks)
    ticklabels{ii} = sprintf('$$%d \\%%$$',ticks(ii));
end

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.915,0.1,0.015,0.4],'Ticks',ticks,'TickLabels',ticklabels);
cb_title.String = '';
set(hcb.Title,cb_title);

% 
axes('Position',[0.10,0.1,0.01,0.85]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'Pulse width $$ PW \: \mathrm{(ms)}$$'},'Interpreter','latex','Color','k','FontSize', 18);

%%  HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH
model_name='UF_TrmB_HH';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PPWW = compiled_results.PPWW;
% th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
th_per_diff_MCE(end,:) = max_diff;

%--------------------------------------------------------------------------
h_ax(1) = axes('Position',[0.11,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm1);
caxis(gca,E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind = find(ALPHA(:,1)<=(85*pi/180));
[C,h]=contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^(E_log_cnt_min:E_log_cnt_max-1),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',500);
ind = find(ALPHA(:,1)>=(80*pi/180));
contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^(E_log_cnt_min:E_log_cnt_max),'-','LineWidth',1.5,'Color',[1,1,1]*0);

title(gca, '\textbf{HH}','Interpreter','latex');

%--------------------------------------------------------------------------
h_ax(2) = axes('Position',[0.11,0.1,0.25,0.40]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind = find(PPWW(1,:)>=(0.02));
contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-45,-15,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind = find(PPWW(1,:)<=(0.021));
[C,h]=contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-45,-15,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex');

%%  RMG str RMG str RMG str RMG str RMG str RMG str RMG str RMG str RMG str
model_name='UF_TrmB_RMG';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PPWW = compiled_results.PPWW;
% th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
th_per_diff_MCE(end,:) = max_diff;

%--------------------------------------------------------------------------
h_ax(3) = axes('position',[0.38,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm1);
caxis(E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind = find(ALPHA(:,1)<=(87*pi/180));
[C,h]=contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^(E_log_cnt_min:E_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',500);
ind = find(ALPHA(:,1)>=(87*pi/180));
contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^(E_log_cnt_min:E_log_cnt_max-1),'-','LineWidth',1.5,'Color',[1,1,1]*0);

title(gca, '\textbf{RMG}','Interpreter','latex');

%--------------------------------------------------------------------------
h_ax(4) = axes('position',[0.38,0.1,0.25,0.4]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind = find(PPWW(1,:)>=(0.01));
contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-3,-3],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind = find(PPWW(1,:)<=(0.01));
[C,h]=contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-3,-3],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex');

xlabel({'E-field orientation $$ \alpha$$'},'Interpreter','latex','FontSize',18); 

%%  RMG syn RMG syn RMG syn RMG syn RMG syn RMG syn RMG syn RMG syn RMG syn
model_name='UF_TrmS_RMG';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PPWW = compiled_results.PPWW;
% th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
th_per_diff_MCE(end,:) = max_diff;

%--------------------------------------------------------------------------
h_ax(5) = axes('position',[0.65,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm1);
caxis(E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind = find(ALPHA(:,1)<=(90*pi/180));
[C,h]=contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^(E_log_cnt_min:E_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',500);
[C,h]=contour(ALPHA(ind,:)/pi*180,log10(PPWW(ind,:)),abs(th_MCE(ind,:))/10,10.^([E_log_cnt_max,E_log_cnt_max]-1),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex');

title(gca, '\textbf{RMG} with synapse','Interpreter','latex');

%--------------------------------------------------------------------------
h_ax(6) = axes('position',[0.65,0.1,0.25,0.4]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind = find(PPWW(1,:)>=(0.021));
contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-50,-20,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind = find(PPWW(1,:)<=(0.021));
[C,h]=contour(ALPHA(:,ind)/pi*180,log10(PPWW(:,ind)),th_per_diff_MCE(:,ind),[-50,-20,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,[-50,-20,-5],'FontSize',14,'Interpreter','latex');

%--------------------------------------------------------------------------
%%
set(h_ax,'XTick', xtick,'YTick', ytick);
set(h_ax(1:2:5),'XTickLabel',{})
set(h_ax(2:2:6),'XTickLabel',xticklabel);
set(h_ax(1:2),'YTickLabel',yticklabel);
set(h_ax(3:6),'YTickLabel',{});
axis(h_ax,[-0.2,90.2,-3.01,1.01]); %axis equal;

set(h_ax,format_axis);
set([h_ax.XLabel,h_ax.YLabel ], format_axis_label);
set([h_ax.Title], format_title);

%%
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(imind,cm,'UF_Trm_combined.tif','tif','WriteMode','overwrite', 'Resolution',300,'Compression','none');
saveas(h_f,'UF_Trm_combined.fig');
close(h_f);