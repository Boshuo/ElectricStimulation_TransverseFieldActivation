clearvars;close all;

figure_format;

h_f = figure;
set(h_f,'Position',[00 00 1500 600],'Color',format_figure.Color);

xtick = 0 : 15 : 135;
xticklabel = cell(length(xtick),1);
for ii = 1 : length(xticklabel)
    if abs(round(xtick(ii)/15)-xtick(ii)/15) < 1e-2
        xticklabel{ii} = ['$$',num2str((xtick(ii)),'%d'),'^{\circ}$$'];
    else
        xticklabel{ii} = ' ';
    end
end

E_log_max = 6.5;	E_log_cnt_max = fix(E_log_max);
E_log_min = 0.75;	E_log_cnt_min = fix(E_log_min);
E_log_lvl = E_log_min:0.01:E_log_max;
E_log_range = [E_log_min,E_log_max];

ticks = E_log_cnt_min : E_log_cnt_max;
ticklabels = cell(size(ticks));
for ii = 1 : length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',ticks(ii));
end
ticks0 = [ticks,E_log_max];
ticklabels0 = [ticklabels,'$$+\infty$$'];

max_diff = -100;
per_diff_lvl = [max_diff:5:-30,-25:2.5:0];

per_ticks = max_diff : 25 : 0;
per_ticklabels = cell(size(per_ticks));
for ii = 1 : length(per_ticks)
    per_ticklabels{ii} = sprintf('$$%d \\%%$$',per_ticks(ii));
end

%%
model_name='UF_BndO_RMG';
filename = fullfile('..',model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PPWW = compiled_results.PPWW;
th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
ind = isnan(th_CE);
th_CE(ind) = 10^(E_log_max)*10;
th_per_diff_MCE(ind) = max_diff;
ind = isnan(th_MCE);
th_MCE(ind) = 10^(E_log_max)*10;
th_per_diff_MCE(ind) = max_diff;

%--------------------------------------------------------------------------
h_ax(1) = axes('Position',[0.1,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm0);
caxis(gca,E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_CE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70/180*pi));
ind_PW =find(PPWW(1,:)<=1);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);
ind_PW =find(PPWW(1,:)>=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(70/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-3)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW =find(PPWW(1,:)<=5);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-3)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(130/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_alpha = find(ALPHA(:,1)<=(114/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
text(114 ,-0.58 , num2str(10.^(E_log_cnt_max-2)),'FontSize',13,'Interpreter','latex','HorizontalAlign','Left','Rotation',-28);

xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
ylabel({'Pulse duration $$ PW \: \mathrm{(ms)}$$'},'Interpreter','latex');
title({'Conventional cable equation','Threshold $${E}$$ in $$\mathrm{V/m}$$'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',ticks0,'TickLabels',ticklabels0);
set(hcb,format_cb);

%--------------------------------------------------------------------------
h_ax(2) = axes('Position',[0.4,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm0);
caxis(gca,E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70/180*pi));
ind_PW =find(PPWW(1,:)<=1);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);
ind_PW =find(PPWW(1,:)>=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(70/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW =find(PPWW(1,:)<=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-3)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(130/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_alpha = find(ALPHA(:,1)<=(114/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
text(114 ,-1.48 , num2str(10.^(E_log_cnt_max-2)),'FontSize',13,'Interpreter','latex','HorizontalAlign','Left','Rotation',-4)

xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
title({'Modified Cable equation','Threshold $${E}$$ in $$\mathrm{V/m}$$'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',ticks0,'TickLabels',ticklabels0);
set(hcb,format_cb);

%--------------------------------------------------------------------------
h_ax(3) = axes('Position',[0.7,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind_PW = find(PPWW(1,:)>=(0.1));
ind_alpha= find(ALPHA(:,1)<=113*pi/180);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-5,-30,-50],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW = find(PPWW(1,:)>=(0.06));
ind_alpha= find(ALPHA(:,1)>=113*pi/180);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-5,-30,-50],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex');

ind_PW = find(PPWW(1,:)<=(0.14));
ind_alpha= find(ALPHA(:,1)<=120*pi/180);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-30,-30],'-','LineWidth',2,'Color',[1,1,1]*0);

plot([45,45],[-3,1],'k:','LineWidth',1.5);
plot([90,90],[-3,1],'k:','LineWidth',1.5);


xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
title({'Percentage difference of threshold','modified vs. conventional'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',per_ticks,'TickLabels',per_ticklabels);
set(hcb,format_cb);


%%
set(h_ax,'XTick', xtick,'YTick', ytick,'XTickLabel',xticklabel,'YTickLabel',yticklabel);
axis(h_ax,[-0.4,135.4,-3.01,1.01]); 
set(h_ax,format_axis);
set([h_ax.XLabel,h_ax.YLabel ], format_axis_label);
set([h_ax.Title], format_title);

%%
filename = fullfile([model_name,'_compiled_result']);
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(imind,cm,[filename,'.tif'],'tif','WriteMode','overwrite', 'Resolution',300);
saveas(h_f,[filename,'.fig']);
close(h_f);
 