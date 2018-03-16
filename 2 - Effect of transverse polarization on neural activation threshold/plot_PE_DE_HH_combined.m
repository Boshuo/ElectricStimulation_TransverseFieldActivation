clearvars;close all;
addpath('Plot functions');
figure_format;
rmpath('Plot functions');

h_f = figure;
set(h_f,'Position',[0,0,1500,800],format_figure);

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

max_diff = -100;
per_diff_lvl = (max_diff:0.2:0);

%
axes('Position',[0.08,0.55,0.88-0.27,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{Current threshold}','modified CE'},'Interpreter','latex','Color','k','FontSize', 20);

caxis(gca,I_log_range);
colormap(gca,cm1);

ticks = I_log_cnt_min : I_log_cnt_max;
ticklabels = cell(size(ticks));
for ii = 1 : length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',ticks(ii));
end

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.915-0.27,0.55,0.015,0.4],'Ticks',ticks,'TickLabels',ticklabels);
cb_title.String = '$$\mathrm{mA}$$';
set(hcb.Title,cb_title);

%
axes('Position',[0.08,0.1,0.88,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{Percentage change}','vs. conventional CE'},'Interpreter','latex','Color','k','FontSize', 20);

caxis(gca,[max_diff,0]);
colormap(gca,cm2);
c_1 = -50;

ticks = max_diff : 25 : 0;
ticklabels = cell(size(ticks));

for ii = 1 : length(ticks)
    ticklabels{ii} = sprintf('$$%d \\%%$$',ticks(ii));
end

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.915-0.27,0.1,0.015,0.4],'Ticks',ticks,'TickLabels',ticklabels);
cb_title.String = '';
set(hcb.Title,cb_title);

% 
axes('Position',[0.10,0.1,0.01,0.85]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'Pulse width $$ PW \: \mathrm{(ms)}$$'},'Interpreter','latex','Color','k','FontSize', 18);

axes('Position',[0.105,0.08,0.525,0.01]);
set(gca,'XColor','w','YColor','w');
xlabel({'Axon\textendash electrode distance $$ H \: \mathrm{(mm)}$$'},'Interpreter','latex','Color','k','FontSize', 18); 

overlay_linear = 1;
%%  PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE PE 
model_name = 'PE_A_HH';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

HH = compiled_results.HH;
PPWW = compiled_results.PPWW;
% th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
%--------------------------------------------------------------------------
h_ax(1) = axes('Position',[0.11,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm1);  caxis(I_log_range);

contourf(log10(HH/1e-1),log10(PPWW),log10(abs(th_MCE)/1e3),I_log_lvl,'LineStyle','none');

[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^(I_log_cnt_min+1:I_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',350);
[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^((I_log_cnt_max-1)*[1,1]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex','LabelSpacing',350);

title(gca, 'Point source - \textbf{HH}','Interpreter','latex');

%--------------------------------------------------------------------------
h_ax(2) = axes('Position',[0.11,0.1,0.25,0.40]);
box on; hold on;
colormap(gca,cm2);  caxis(gca,[max_diff,0]);

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

%%  DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE DE
model_name = 'DE_A_HH';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

HH = compiled_results.HH;
PPWW = compiled_results.PPWW;
% th_CE = compiled_results.th_CE;
th_MCE = compiled_results.th_MCE;
th_per_diff_MCE = compiled_results.th_per_diff_MCE;
%--------------------------------------------------------------------------
h_ax(3) = axes('position',[0.38,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm1);  caxis(I_log_range);

contourf(log10(HH/1e-1),log10(PPWW),log10(abs(th_MCE)/1e3),I_log_lvl,'LineStyle','none');

[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^(I_log_cnt_min*[1,1]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',450);
[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^(I_log_cnt_min+1:I_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',350);
[C,h]=contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^(I_log_cnt_max-1:I_log_cnt_max),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex','LabelSpacing',350);

contour(log10(HH/1e-1),log10(PPWW),abs(th_MCE)/1e3,10.^ceil([I_log_max,I_log_max]),'-','LineWidth',1.5,'Color',[1,1,1]*0);

title(gca, 'Disk electrode - \textbf{HH}','Interpreter','latex');

%--------------------------------------------------------------------------
h_ax(4) = axes('position',[0.38,0.1,0.25,0.4]);
box on; hold on;
colormap(gca,cm2);  caxis(gca,[max_diff,0]);

contourf(log10(HH/1e-1),log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind_PW = find(PPWW(1,:)<=0.07 & PPWW(1,:)>=0.004);
[C,h]=contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',14,'Interpreter','latex');
ind_PW = find(PPWW(1,:)>=0.06 );
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW = find( PPWW(1,:)<=0.005);
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[c_1:10:-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_PW = find( PPWW(1,:)<=0.04  & PPWW(1,:)>=0.003);
[C,h]=contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-60,-60],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex');
ind_PW = find( PPWW(1,:)<=0.004);
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-60,-60],'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_PW = find( PPWW(1,:)<=0.04 );
contour(log10(HH(:,ind_PW)/1e-1),log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-70,-70],'-','LineWidth',1.5,'Color',[1,1,1]*0);

%%
if overlay_linear
    for kk = 1: 2
        if kk == 1
            h_temp = h_ax(2);
            filename = fullfile('..','1 - Validation of modified cable equation','2 - Distributed linear cable - MATLAB','Results','PE_A_Lin.mat');
            load(filename,'PW','results');
        else
            h_temp = h_ax(4);
            filename = fullfile('..','1 - Validation of modified cable equation','2 - Distributed linear cable - MATLAB','Results','DE_A_Lin.mat');
            load(filename,'PW','results');
        end
        
        H = [results.H];
        per_diff = zeros(length(PW),length(H));
        
        for ii = 1 : length(H)  % Metric used in Neu 2016
            per_diff(:,ii)= - (2 * results(ii).ExR ) ./ (results(ii).phi_m_bar + 2 * results(ii).ExR ) *100;
        end
        
        ind_PW = find((PW>=1e-3)&(PW<=10));
        ind_H = find(H>30e-4 & H<0.1);
        contour(h_temp, log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.05,0.1:0.1:0.3,0.7,0.9]*100,'--','LineWidth',1,'Color','r');
        contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.5,0.5]*100,'--','LineWidth',2,'Color','r');
        
        ind_H = find(H<0.4 & H>0.05);
        [C,h]=contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.050,0.1:0.1:0.3,0.7,0.9]*100,'--','LineWidth',1,'Color','r');
        clabel(C,h,'FontSize',14,'Interpreter','latex','Labelspacing',400,'Color','r');
        [C,h]=contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.5,0.5]*100,'--','LineWidth',2,'Color','r');
        clabel(C,h,'FontSize',14,'Interpreter','latex','LabelSpacing',400,'Color','r');
        
        ind_H = find(H>0.3 & H < 3.2);
        contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.050,0.1:0.1:0.3,0.7,0.9]*100,'--','LineWidth',1,'Color','r');
        contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),-[0.5,0.5]*100,'--','LineWidth',2,'Color','r');
        [C,h]=contour(h_temp,log10(H(ind_H)/1e-1),log10(PW(ind_PW)),per_diff(ind_PW,ind_H),[-0.99,-0.99]*100,'--','LineWidth',1,'Color','r');
        clabel(C,h,'FontSize',14,'Interpreter','latex','Labelspacing',100,'Color','r');

    end
end
%--------------------------------------------------------------------------
%%
set(h_ax,'XTick', xtick,'YTick', ytick);
set(h_ax(1:2:4),'XTickLabel',{})
set(h_ax(2:2:4),'XTickLabel',xticklabel);
set(h_ax(1:2),'YTickLabel',yticklabel);
set(h_ax(3:4),'YTickLabel',{});
axis(h_ax,log10([31e-3*0.98,31*1.02,1e-3*0.98,10*1.02])); %axis equal;

set(h_ax,format_axis);
set([h_ax.XLabel,h_ax.YLabel ], format_axis_label);
set([h_ax.Title], format_title);

%%
im = frame2im(getframe(h_f));
[imind,cm] = rgb2ind(im(:,1:1100*2,:),256);
imwrite(    imind,cm,'PE_DE_HH_combined.tif','tif','WriteMode','overwrite', 'Resolution',150,'Compression','none');
saveas(h_f,'PE_DE_HH_combined.fig');
close(h_f);