clearvars;close all;

figure_format;

h_f = figure;
set(h_f,'Position',[00 00 1500 600],'Color',format_figure.Color);

xtick = log10([kron([1e-3,1e-2,1e-1,1,10],[1,2,3,4,5,6,7,8,9]),100]);
xticklabel = cell(length(xtick),1);
for ii = 1 : length(xticklabel)
    if abs(round(xtick(ii)) - xtick(ii)) < 1e-2
        xticklabel{ii} = ['$$10^{',num2str((xtick(ii)),'%d'),'}$$'];
    else
        xticklabel{ii} = ' ';
    end
end

I_log_max = 6.05;   I_log_cnt_max = fix(I_log_max);
I_log_min = -2.05;  I_log_cnt_min = fix(I_log_min);
I_log_lvl = I_log_min:0.05:I_log_max;
I_log_range = [I_log_min,I_log_max];

ticks = I_log_cnt_min : I_log_cnt_max;
ticklabels = cell(size(ticks));
for ii = 1 : length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',ticks(ii));
end

max_diff = -100;
per_diff_lvl = (max_diff:0.2:0);
c_1 = -50;

per_ticks = max_diff : 25 : 0;
per_ticklabels = cell(size(per_ticks));

for ii = 1 : length(per_ticks)
    per_ticklabels{ii} = sprintf('$$%d \\%%$$',per_ticks(ii));
end

%%
model_name = 'PE_A_HH';
filename = fullfile('..',model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

HH = compiled_results.HH;
PPWW = compiled_results.PPWW;
th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
%--------------------------------------------------------------------------
h_ax(1) = axes('Position',[0.1,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm1);  caxis(I_log_range);

contourf(log10(HH/1e-1),log10(PPWW),log10(abs(th_CE)/1e3),I_log_lvl,'LineStyle','none');

[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_CE)/1e3,10.^(I_log_cnt_min+1:I_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',350);
[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_CE)/1e3,10.^((I_log_cnt_max-1)*[1,1]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex','LabelSpacing',350);

xlabel('Axon-electrode distance $$ h \: \rm{ (mm)}$$','Interpreter','latex');
ylabel('Pulse duration $$PW \: \rm{(ms)}$$','Interpreter','latex');
title({'Conventional cable equation','Threshold $$|I_{\rm{th}}|$$ in $$\rm{mA}$$'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',ticks,'TickLabels',ticklabels);
set(hcb,format_cb);

%--------------------------------------------------------------------------
h_ax(2) = axes('Position',[0.4,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm1);  caxis(I_log_range);

contourf(log10(HH/1e-1),log10(PPWW),log10(abs(th_MCE)/1e3),I_log_lvl,'LineStyle','none');

[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^(I_log_cnt_min+1:I_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',350);
[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^((I_log_cnt_max-1)*[1,1]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex','LabelSpacing',350);

xlabel('Axon-electrode distance $$ h \: \rm{ (mm)}$$','Interpreter','latex');
title({'Modified cable equation','Threshold $$|I_{\rm{th}}|$$ in $$\rm{mA}$$'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',ticks,'TickLabels',ticklabels);
set(hcb,format_cb);

%--------------------------------------------------------------------------
h_ax(3) = axes('Position',[0.7,0.1,0.25,0.75]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(log10(HH/1e-1),log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind_PW = find(PPWW(1,:)<=0.07 & PPWW(1,:)>=0.004);
[C,h]=contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex');
ind_PW = find(PPWW(1,:)>=0.06 );
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW = find( PPWW(1,:)<=0.005);
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_PW = find( PPWW(1,:)<=0.04  & PPWW(1,:)>=0.0020);
[C,h]=contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-60,-60],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex');
ind_PW = find( PPWW(1,:)<=0.0022);
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-60,-60],'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_PW = find( PPWW(1,:)<=0.04 );
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-70,-70],'-','LineWidth',1.5,'Color',[1,1,1]*0);

xlabel('Axon-electrode distance $$ h \: \rm{ (mm)}$$','Interpreter','latex');
title({'Percentage difference of threshold','modified vs. conventional'},'Interpreter','latex');

hcb = colorbar(gca,'SouthOutside'); 
set(hcb,'Ticks',per_ticks,'TickLabels',per_ticklabels);
set(hcb,format_cb);

%%
set(h_ax,'XTick', xtick,'YTick', ytick,'XTickLabel',xticklabel,'YTickLabel',yticklabel);
axis(h_ax,log10([31e-3*0.98,31*1.02,1e-3*0.98,10*1.02])); 
set(h_ax,format_axis);
set([h_ax.XLabel,h_ax.YLabel ], format_axis_label);
set([h_ax.Title], format_title);

%%
filename = fullfile([model_name,'_compiled_result']);
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(imind,cm,[filename,'.tif'],'tif','WriteMode','overwrite', 'Resolution',300);
saveas(h_f,[filename,'.fig']);
close(h_f);
