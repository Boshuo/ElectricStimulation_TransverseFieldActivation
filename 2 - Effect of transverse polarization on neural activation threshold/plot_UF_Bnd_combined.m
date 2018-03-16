clearvars;close all;
addpath('Plot functions');
figure_format;
rmpath('Plot functions');

xtick = 0 : 15 : 135;
xticklabel = cell(length(xtick),1);
for ii = 1 : length(xticklabel)
    if abs(round(xtick(ii)/15)-xtick(ii)/15) < 1e-2
        xticklabel{ii} = ['$$',num2str((xtick(ii)),'%d'),'^{\circ}$$'];
    else
        xticklabel{ii} = ' ';
    end
end

E_log_max = 6.25;	E_log_cnt_max = fix(E_log_max);
E_log_min = 0.75;	E_log_cnt_min = fix(E_log_min);
E_log_lvl = E_log_min:0.01:E_log_max;
E_log_range = [E_log_min,E_log_max];

max_diff = -100;
per_diff_lvl = [max_diff:5:-30,-25:2.5:0];

h_f1 = figure;
set(h_f1,'Position',[0,0,1500,800],format_figure);

axes('Position',[0.08,0.55,0.1,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{E-field threshold}','conventional CE'},'Interpreter','latex','Color','k','FontSize', 20);

axes('Position',[0.08,0.1,0.1,0.4]);
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{E-field threshold}','modified CE'},'Interpreter','latex','Color','k','FontSize', 20);

axes('Position',[0.10,0.1,0.1,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'Pulse width $$ PW \: \mathrm{(ms)}$$'},'Interpreter','latex','Color','k','FontSize', 18);

caxis(gca,E_log_range);
colormap(gca,cm0);

ticks = E_log_cnt_min : E_log_cnt_max;
ticklabels = cell(size(ticks));
for ii = 1 : length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',ticks(ii));
end
ticks = [ticks,E_log_max];
ticklabels = [ticklabels,{'$$+\infty$$'}];

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.94-0.275,0.1,0.015,0.85],'Ticks',ticks,'TickLabels',ticklabels,'TickLength',0.01);
cb_title.String = '$$\mathrm{V/m}$$';
set(hcb.Title,cb_title);

%%
h_f2 = figure;
set(h_f2,'Position',[0,0,1500,800],'Color','w');

axes('Position',[0.08,0.55,0.1,0.4]); 
set(gca,format_axis);set(gca,'XColor','w','YColor','w');
ylabel({'\textbf{Percentage change}','modified vs. conventional CE'},'Interpreter','latex','Color','k','FontSize', 20);

caxis(gca,[max_diff,0]);
colormap(gca,cm2);

ticks = max_diff : 25 : 0;
ticklabels = cell(size(ticks));

for ii = 1 : length(ticks)
    ticklabels{ii} = sprintf('$$%d \\%%$$',ticks(ii));
end

hcb = colorbar(gca); 
set(hcb,format_cb);
set(hcb,'position',[0.94-0.275,0.55,0.015,0.4],'Ticks',ticks,'TickLabels',ticklabels);
cb_title.String = '';
set(hcb.Title,cb_title);


%%  HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH HH
model_name = 'UF_BndO_HH';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PPWW = compiled_results.PPWW;
th_CE_ortho = compiled_results.th_CE;
th_MCE_ortho = compiled_results.th_MCE;

model_name = 'UF_BndA_HH';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

th_CE_anti = compiled_results.th_CE;
th_MCE_anti = compiled_results.th_MCE;


th_CE   = min(th_CE_ortho,th_CE_anti);
th_MCE  = min(th_MCE_ortho,th_MCE_anti);
th_per_diff_MCE = ( th_MCE ./ th_CE - 1) * 100;

th_CE(isnan(th_CE)) = 1e8;
th_per_diff_MCE(isnan(th_per_diff_MCE)) = -100;
th_per_diff_MCE(abs(th_per_diff_MCE) <= 0.5) =0;

%--------------------------------------------------------------------------
figure(h_f1);
h_ax(1) = axes('Position',[0.125,0.55,0.25,0.40]);
box on; hold on;
colormap(gca,cm0);
caxis(gca,E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_CE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70*pi/180));
[C,h]=contour(ALPHA(ind_alpha,:)/pi*180,log10(PPWW(ind_alpha,:)),abs(th_CE(ind_alpha,:))/10,10.^(E_log_cnt_min:(E_log_cnt_max-2)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);
ind_alpha = find(ALPHA(:,1)>=(70*pi/180));
ind_PW = find(PPWW(1,:)>=0.002);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW = find(PPWW(1,:)<0.015);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PPWW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^([1,1]*(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',00);


ind_dir = zeros(size(ALPHA));
for ii = 1:numel(ALPHA)
    if isnan(th_CE_ortho(ii)) && isnan(th_CE_anti(ii))      % Neither has threshold
        ind_dir(ii) = 2;
    elseif isnan(th_CE_anti(ii))                            % Only orthodromic has threshold
        ind_dir(ii) = 1;
    elseif  abs(th_CE_ortho(ii)/th_CE_anti(ii)-1) <= 1e-2    % Same threshold
        ind_dir(ii) = 0;
    elseif abs(th_CE_ortho(ii)) * 0.995 < abs(th_CE_anti(ii)) * 1.005   % Smaller orthodromic threshold
        ind_dir(ii) = 0.5;
    end
end
contour(ALPHA/pi*180,log10(PPWW),ind_dir,[1.5,1.5],'-','LineWidth',2,'Color',[1,1,1]*1);
contour(ALPHA/pi*180,log10(PPWW),ind_dir,[0.52,0.52],'--','LineWidth',2,'Color',[1,1,1]*1);
ind_PW =find(PPWW(1,:) > 3);
contour(ALPHA(:,ind_PW)/pi*180,log10(PPWW(:,ind_PW)),ind_dir(:,ind_PW),[0.48,0.48],'--','LineWidth',2,'Color',[1,1,1]*0.8);

title(gca, '\textbf{HH}','Interpreter','latex');

plot([116,125], [-2,-2],'-','LineWidth',1.5,'Color','w');
patch([120,120,116],[-1.92,-2.08,-2],'w','LineStyle','none');
text(126,-2,'$$\alpha_{2}$$','Interpreter','latex','Color','w','FontSize',16, 'VerticalAlignment','Middle');

%--------------------------------------------------------------------------
h_ax(2) = axes('Position',[0.125,0.1,0.25,0.4]);
box on; hold on;
colormap(gca,cm0);
caxis(gca,E_log_range);

contourf(ALPHA/pi*180,log10(PPWW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70/180*pi));
[C,h]=contour(ALPHA(ind_alpha,:)/pi*180,log10(PPWW(ind_alpha,:)),abs(th_MCE(ind_alpha,:))/10,10.^(E_log_cnt_min:E_log_cnt_max),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);

ind_alpha = find(ALPHA(:,1)>=(70/180*pi));
contour(ALPHA(ind_alpha,:)/pi*180,log10(PPWW(ind_alpha,:)),abs(th_MCE(ind_alpha,:))/10,10.^(E_log_cnt_min:E_log_cnt_max-2),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(95/180*pi));
[C,h]=contour(ALPHA(ind_alpha,:)/pi*180,log10(PPWW(ind_alpha,:)),abs(th_MCE(ind_alpha,:))/10,10.^((E_log_cnt_max-1)*[1,1]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',12,'Interpreter','latex','LabelSpacing',250);

ind_dir = zeros(size(ALPHA))-4;
for ii = 1:numel(ALPHA)
    if isnan(th_MCE_ortho(ii)) && isnan(th_MCE_anti(ii))    % Neither has threshold
        ind_dir(ii) = 10;
    elseif isnan(th_MCE_ortho(ii))                          % Only antidromic has threshold
        ind_dir(ii) = 3;
    elseif  abs(th_MCE_ortho(ii)/th_MCE_anti(ii)-1) <= 1e-2  % Same threshold
        ind_dir(ii) = 0;
    elseif abs(th_MCE_anti(ii)) * 0.995 < abs(th_MCE_ortho(ii)) * 1.005   % Larger orthodromic threshold
        ind_dir(ii) = 2;
    elseif abs(th_MCE_ortho(ii)) * 0.995 < abs(th_MCE_anti(ii)) * 1.005   % Smaller orthodromic threshold
        ind_dir(ii) = 1;    
    end
end

diff_th_MCE = diff(th_MCE,1);

for ii = 1 : size(th_MCE,2)
    ind_neg = find(diff_th_MCE(:,ii) < 0,1,'first');
    ind_pos = find(ind_dir(1:ind_neg,ii) == 1 ,1 ,'last');
    ind_dir(ind_pos + 1 : ind_neg - 1,ii) = 1;
    ind_dir(ind_neg:end,ii) = 2;
end

ind_PW = find(PPWW(1,:) >= 0.04);
contour(ALPHA(:,ind_PW)/pi*180,log10(PPWW(:,ind_PW)),ind_dir(:,ind_PW),[0.5,0.5],'--','LineWidth',2,'Color',[1,1,1]*0.8);
contour(ALPHA/pi*180,log10(PPWW),ind_dir,[1.5,1.5],'--','LineWidth',2,'Color',[1,1,1]);

plot([116,125]-7, [-2,-2],'-','LineWidth',1.5,'Color','w');
patch([120,120,116]-7,[-1.92,-2.08,-2],'w','LineStyle','none');
text(126-7,-2,'$$\alpha_{1}$$','Interpreter','latex','Color','w','FontSize',16, 'VerticalAlignment','Middle');

%% --------------------------------------------------------------------------
figure(h_f2);
axes('Position',[0.125,0.525,0.525,0.01]);
set(gca,'XColor','w','YColor','w');
xlabel({'E-field orientation $$ \alpha$$'},'Interpreter','latex','Color','k','FontSize',18); 

h_ax(3) = axes('Position',[0.125,0.55,0.25,0.40]);
box on; hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PPWW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

plot([45,45],[-3,1],'k:','LineWidth',1.5);
plot([90,90],[-3,1],'k:','LineWidth',1.5);

ind_PW = find(PPWW(1,:)>=(0.03));
contour(ALPHA(:,ind_PW)/pi*180,log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-50,-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_PW = find(PPWW(1,:)<=(0.04));
[C,h]=contour(ALPHA(:,ind_PW)/pi*180,log10(PPWW(:,ind_PW)),th_per_diff_MCE(:,ind_PW),[-50,-10,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',300);

%%  RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG RMG
model_name = 'UF_BndO_RMG';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

ALPHA = compiled_results.ALPHA;
PW = compiled_results.PPWW;
th_CE_ortho = compiled_results.th_CE;
th_MCE_ortho = compiled_results.th_MCE;

model_name = 'UF_BndA_RMG';
filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

th_CE_anti = compiled_results.th_CE;
th_MCE_anti = compiled_results.th_MCE;

th_CE = min(th_CE_ortho,th_CE_anti);
th_MCE = min(th_MCE_ortho,th_MCE_anti);
th_per_diff_MCE = ( th_MCE ./ th_CE - 1) *100;

th_CE(isnan(th_CE)) = 1e8;
th_per_diff_MCE(isnan(th_per_diff_MCE)) = -100;
th_per_diff_MCE(abs(th_per_diff_MCE) < 1) =0;

%--------------------------------------------------------------------------
figure(h_f1);
h_ax(4) = axes('position',[0.4,0.55,0.25,0.4]);
box on; hold on;
colormap(gca,cm0);
caxis(E_log_range);

contourf(ALPHA/pi*180,log10(PW),log10(abs(th_CE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70/180*pi));
ind_PW =find(PW(1,:)<=1);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);
ind_PW =find(PW(1,:)>=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(70/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW =find(PW(1,:)<=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_CE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-2)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

title(gca, '\textbf{RMG}','Interpreter','latex');

ind_dir = zeros(size(ALPHA));
for ii = 1:numel(ALPHA)
    if isnan(th_CE_ortho(ii)) && isnan(th_CE_anti(ii))      % Neither has threshold
        ind_dir(ii) = 2;
    elseif  abs(th_CE_ortho(ii)/th_CE_anti(ii)-1)<1e-2      % Same threshold
        ind_dir(ii) = 0;
    elseif abs(th_CE_ortho(ii))*0.995 < abs(th_CE_anti(ii)) * 1.005     % Smaller orthodromic threshold
        ind_dir(ii) = 1;
    end
end

for ii = 1 : size(th_CE,2)
    ind = find(ind_dir(:,ii) == 1 ,1 ,'last');
    ind_dir(ind + 1 : end,ii) = 3;
end
for ii = 1 : size(th_CE,1)
    ind = find(ind_dir(ii,:) == 2 ,1 ,'last');
    ind_dir(ii,ind + 1: end) = 3;
end

ind_0 = (ind_dir == 0);
ind_1 = (ind_dir == 1);
ind_2 = (ind_dir == 2);
ind_3 = (ind_dir == 3);

ind_dir(ind_2) = 1;
ind_dir(~ind_2) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.5,0.5],'-','LineWidth',2,'Color',[1,1,1]*1);
ind_dir(ind_1) = 1;
ind_dir(~ind_1) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.75,0.75],'--','LineWidth',2,'Color',[1,1,1]*0.8);
ind_dir(ind_3) = 1;
ind_dir(~ind_3) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.75,0.75],':','LineWidth',2,'Color',[1,1,1]*1);

plot([116,125]-2.5, [-2,-2],'-','LineWidth',1.5,'Color','w');
patch([120,120,116]-2.5,[-1.92,-2.08,-2],'w','LineStyle','none');
text(125-2.5,-2.05,'$$\alpha_{2}$$','Interpreter','latex','Color','w','FontSize',16, 'VerticalAlignment','Top');

%% --------------------------------------------------------------------------
h_ax(5) = axes('position',[0.4,0.1,0.25,0.4]);
box on; hold on;
colormap(gca,cm0);
caxis(E_log_range);

contourf(ALPHA/pi*180,log10(PW),log10(abs(th_MCE)/10),E_log_lvl,'LineStyle','none');

ind_alpha = find(ALPHA(:,1)<=(70/180*pi));
ind_PW =find(PW(1,:)<=1);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','LabelSpacing',500);
ind_PW =find(PW(1,:)>=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(70/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW =find(PW(1,:)<=1);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^(E_log_cnt_min:(E_log_cnt_max-3)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

ind_alpha = find(ALPHA(:,1)>=(130/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_alpha = find(ALPHA(:,1)<=(114/180*pi));
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),abs(th_MCE(ind_alpha,ind_PW))/10,10.^([(E_log_cnt_max-2),(E_log_cnt_max-2)]),'-','LineWidth',1.5,'Color',[1,1,1]*0);
text(114 ,-1.48 , num2str(10.^(E_log_cnt_max-2)),'FontSize',13,'Interpreter','latex','HorizontalAlign','Left','Rotation',-4)

ind_dir(ind_0) = 0;
ind_dir(ind_1) = 1;
ind_dir(ind_2) = 2;
ind_dir(ind_3) = 3;

ind_dir(th_per_diff_MCE < -5) = 2;
ind_0 = (ind_dir ==0);
ind_1 = (ind_dir ==1);
ind_2 = (ind_dir ==2);
ind_3 = (ind_dir ==3);

ind_dir(ind_2) = 1;
ind_dir(~ind_2) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.5,0.5],'--','LineWidth',2,'Color',[1,1,1]*1);
ind_dir(ind_1) = 1;
ind_dir(~ind_1) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.75,0.75],'--','LineWidth',2,'Color',[1,1,1]*0.8);
ind_dir(ind_3) = 1;
ind_dir(~ind_3) = 0;
contour(ALPHA/pi*180,log10(PW),ind_dir,[0.75,0.75],':','LineWidth',2,'Color',[1,1,1]*1);

plot([116,125]-2.5, [-2,-2],'-','LineWidth',1.5,'Color','w');
patch([120,120,116]-2.5,[-1.92,-2.08,-2],'w','LineStyle','none');
text(125-2.5,-2.05,'$$\alpha_{1}$$','Interpreter','latex','Color','w','FontSize',16, 'VerticalAlignment','Top');

%%
figure(h_f2);
h_ax(6) = axes('Position',[0.4,0.55,0.25,0.40]);
box on;hold on;
colormap(gca,cm2);
caxis(gca,[max_diff,0]);

contourf(ALPHA/pi*180,log10(PW),th_per_diff_MCE,per_diff_lvl,'LineStyle','none');

ind_PW = find(PW(1,:)>=(0.1));
ind_alpha= find(ALPHA(:,1)<=113*pi/180);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-5,-30,-50],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind_PW = find(PW(1,:)>=(0.06));
ind_alpha= find(ALPHA(:,1)>=113*pi/180);
[C,h]=contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-5,-30,-50],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex');

ind_PW = find(PW(1,:)<=(0.14));
ind_alpha= find(ALPHA(:,1)<=120*pi/180);
contour(ALPHA(ind_alpha,ind_PW)/pi*180,log10(PW(ind_alpha,ind_PW)),th_per_diff_MCE(ind_alpha,ind_PW),[-30,-30],'-','LineWidth',2,'Color',[1,1,1]*0);

plot([45,45],[-3,1],'k:','LineWidth',1.5);
plot([90,90],[-3,1],'k:','LineWidth',1.5);

%%
set(h_ax,'XTick', xtick,'YTick', ytick);
set(h_ax([1,2,4,5]),'XTickLabel',{})
set(h_ax([3,6]),'XTickLabel',xticklabel);
set(h_ax(1:3),'YTickLabel',yticklabel);
set(h_ax(4:6),'YTickLabel',{});
axis(h_ax,[-0.4,135.4,-3.015,1.015]); %axis equal;

set(h_ax,format_axis);
set([h_ax.XLabel,h_ax.YLabel ], format_axis_label);
set([h_ax.Title], format_title);


%%
im1 = frame2im(getframe(h_f1));
im2 = frame2im(getframe(h_f2));

im = cat(1,im1(1:725*2,1:1100*2,:),im2(1:450*2,1:1100*2,:));
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,'UF_Bnd_combined.tif','tif','WriteMode','overwrite', 'Resolution',300,'Compression','none');

% saveas(h_f1,'UF_Bnd_combined_th.fig');
% saveas(h_f2,'UF_Bnd_combined_per.fig');
close(h_f1);
close(h_f2);