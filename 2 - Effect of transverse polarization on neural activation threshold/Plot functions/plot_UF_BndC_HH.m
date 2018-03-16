% close all;clc;%clear all;

%% plotting
%%%% figure format
journal_axis.LineWidth = 1.5;
journal_axis.FontSize = 14;
% journal_axis.FontName = 'Times New Roman';
journal_axis.Color = [1 1 1]; % background color of the plot area
journal_axis.TickLabelInterpreter = 'latex';
journal_axis.TickDir = 'out';

journal_axis_label.FontSize = 14;
journal_axis_label.FontName = 'Times New Roman';
journal_axis_label.Interpreter = 'latex';

journal_figure.Color = [1 1 1]; % background color of figure window

journal_trace.LineWidth = 1.5;
% journal_trace.Color = [0 0 0];

journal_title.FontSize = 14;
journal_title.FontName = 'Times New Roman';
journal_title.FontWeight = 'Normal';
journal_title.Interpreter = 'latex';
%%%% end figure format

cm1 = flipud(parula(256));
cm2 = parula(256);
% cm2(1,:) = [1,1,1]*(cm2(1,1) * 0.2126+ cm2(1,2)*0.7158 + cm2(1,3)* 0.0722);

for ii = 256:-1:129
    cm2(ii,:)= brighten(cm2(ii,:),(ii-129)/128);
end

cm2 = flipud(cm2);
%%
model_name='UF_BndO_HH';

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');

alpha_all = compiled_results.ALPHA;
PW_all = compiled_results.PPWW;
th_CE_ortho = compiled_results.th_CE;
th_MCE_ortho = compiled_results.th_MCE;
th_per_diff_MCE_ortho = compiled_results.th_per_diff_MCE;
% 
max_noIP = max(th_CE_ortho(:));
% th_CE_ortho(isnan(th_CE_ortho)) = 1e7;
th_per_diff_MCE_ortho(isnan(th_per_diff_MCE_ortho)) = -100;
th_per_diff_MCE_ortho(abs(th_per_diff_MCE_ortho)<0.05) =0;
%%

model_name='UF_BndA_HH';

filename = fullfile(model_name,[model_name,'_compiled_result.mat']);
load(filename,'compiled_results');


th_CE_anti = compiled_results.th_CE;
th_MCE_anti = compiled_results.th_MCE;
th_per_diff_MCE_anti = compiled_results.th_per_diff_MCE;
% 
max_noIP = max(th_CE_anti(:));
% th_CE_anti(isnan(th_CE_anti)) = 1e7;
th_per_diff_MCE_anti(isnan(th_per_diff_MCE_anti)) = -100;
th_per_diff_MCE_anti(abs(th_per_diff_MCE_anti)<0.05) =0;
%%

th_CE = min(th_CE_ortho,th_CE_anti);
th_MCE = min(th_MCE_ortho,th_MCE_anti);
th_per_diff_MCE = ( th_MCE ./ th_CE - 1) *100;

th_CE(isnan(th_CE)) = 1e8;
th_per_diff_MCE(isnan(th_per_diff_MCE)) = -100;
th_per_diff_MCE(abs(th_per_diff_MCE)<0.05) =0;

%%
figure;
set(gcf,'Position',[00 00 1500 600],'Color',journal_figure.Color);

%%
axes('position',[0.1,0.1,0.25,0.8]);
box on; hold on;

amp_log_max = 6;
amp_log_min = 1;
cm = cm1;
cm(end-8:end,:) = cm(end-8:end,:).*kron((0.95:-0.1:0.15)',[1,1,1]);
colormap(gca,cm);
caxis([amp_log_min,amp_log_max]);
% surf(log10(H_all/1e-1),YY_all/1e-1,0* H_all,log10(abs(th_CE)/10),'LineStyle','-');
contourf(alpha_all/pi*180,log10(PW_all),log10(abs(th_CE)/10),amp_log_min:0.02:amp_log_max,'LineStyle','none');

% plot([109,109],[-3,1],'w:','LineWidth',2);
% plot([45,45],[-3,1],'w:','LineWidth',2);

ind = find(alpha_all(:,1)<=(115*pi/180));
[C,h]=contour(alpha_all(ind,:)/pi*180,log10(PW_all(ind,:)),abs(th_CE(ind,:))/10,10.^round(amp_log_min:(amp_log_max-2)),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',10,'FontWeight','Bold');
ind = find(alpha_all(:,1)<=(115*pi/180));
ind_PW =find(PW_all(1,:)<0.01);
[C,h]=contour(alpha_all(ind,ind_PW)/pi*180,log10(PW_all(ind,ind_PW)),abs(th_CE(ind,ind_PW))/10,10.^([1,1]*round(amp_log_max-1)),'-','LineWidth',1.5,'Color',[1,1,1]*0);

xtick=0:5:135;
xticklabel=cell(length(xtick),1);
for ii=1:length(xticklabel)
    if abs(round(xtick(ii)/15)-xtick(ii)/15) < 1e-2
        xticklabel{ii}= ['$$',num2str(xtick(ii),'%d'),'^{\circ}$$'];
    else
         xticklabel{ii}='';
    end
end

ytick=log10([kron([1e-3,1e-2,1e-1,1],[1,2,3,4,5,6,7,8,9]),10]);
yticklabel=cell(length(ytick),1);
for ii=1:length(yticklabel)
    if abs(round((ytick(ii)))-(ytick(ii))) < 1e-6
        yticklabel{ii}= ['$$10^{',num2str((ytick(ii)),'%d'),'}$$'];
    else
        yticklabel{ii}='';
    end
end

set(gca,'XTick', xtick);
set(gca,'XTickLabel',xticklabel);
set(gca,'YTick',ytick);
set(gca,'YTickLabel',yticklabel);

axis([0,135,-3,1]); %axis equal;

xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
ylabel({'Pulse duration $$ PW \: \mathrm{(ms)}$$'},'Interpreter','latex');
title({'Conventional \textbf{1D} cable equation','Threshold $$\overrightarrow{E}$$ in $$\mathrm{V/m}$$'},'Interpreter','latex');
set(get(gca,'XLabel'), journal_axis_label);
set(get(gca,'YLabel'),journal_axis_label);
set(get(gca,'Title'), journal_title);
set(gca,journal_axis);

% amp_log_max = amp_log_max+0.5;
% amp_log_min = amp_log_min+0.5;
caxis(gca,[amp_log_min,amp_log_max]);
ticks=round(amp_log_min:amp_log_max);
ticks(end) = amp_log_max;
ticklabels=cell(size(ticks));
for ii=1:length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',round(ticks(ii)));
end
ticklabels{end} = '$$+\infty$$';
hcb=colorbar(gca,'SouthOutside'); 
set(hcb,'TickDirection','out');
set(hcb,'TickLabelInterpreter','latex');
set(hcb,'FontSize',14);
set(hcb,'Ticks',ticks);
set(hcb,'TickLabels',ticklabels);

ind_dir = zeros(size(alpha_all));
for ii = 1:numel(alpha_all)
    if isnan(th_CE_ortho(ii)) && isnan(th_CE_anti(ii))
        ind_dir(ii) = 2;
    elseif isnan(th_CE_anti(ii))
        ind_dir(ii) = 1;
    elseif  abs(th_CE_ortho(ii)/th_CE_anti(ii)-1)<1e-2
        ind_dir(ii) = 0;
    elseif abs(th_CE_ortho(ii))<abs(th_CE_anti(ii))
        ind_dir(ii) = 0.5;
    end
end
contour(alpha_all/pi*180,log10(PW_all),ind_dir,[0.55,0.55],'--','LineWidth',2,'Color',[1,1,1]*1);
contour(alpha_all/pi*180,log10(PW_all),ind_dir,[1.5,1.5],'-','LineWidth',2,'Color',[1,1,1]*1);
ind_PW =find(PW_all(1,:)>0.7);
contour(alpha_all(:,ind_PW)/pi*180,log10(PW_all(:,ind_PW)),ind_dir(:,ind_PW),[0.5,0.5],'--','LineWidth',2,'Color',[1,1,1]*0.8);


% return
%%
axes('position',[0.4,0.1,0.25,0.8]);
box on; hold on;

amp_log_max = 6;
amp_log_min = 1;

colormap(gca,cm1);

contourf(alpha_all/pi*180,log10(PW_all),log10(abs(th_MCE)/10),amp_log_min:0.01:amp_log_max,'LineStyle','none');

[C,h]=contour(alpha_all/pi*180,log10(PW_all),abs(th_MCE)/10,10.^round(amp_log_min:amp_log_max),'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',10,'FontWeight','Bold');

set(gca,'XTick', xtick);
set(gca,'XTickLabel',xticklabel);
set(gca,'YTick',ytick);
set(gca,'YTickLabel',yticklabel);

axis([0,135,-3,1]); %axis equal;

xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
% ylabel('Axon vertical placement $$Y \: \mathrm{(mm)}$$');
title({'Cable equation with \textbf{Initial Polarization}','Threshold $$\overrightarrow{E}$$ in $$\mathrm{V/m}$$'},'Interpreter','latex');
set(get(gca,'XLabel'), journal_axis_label);
set(get(gca,'YLabel'),journal_axis_label);
set(get(gca,'Title'), journal_title);
set(gca,journal_axis);

% amp_log_max = amp_log_max+0.5;
% amp_log_min = amp_log_min+0.5;
caxis(gca,[amp_log_min,amp_log_max]);
ticks=round(amp_log_min:amp_log_max);
ticklabels=cell(size(ticks));
for ii=1:length(ticks)
    ticklabels{ii}=sprintf('$$10^{%1.0f}$$',round(ticks(ii)));
end
hcb=colorbar(gca,'SouthOutside'); 
set(hcb,'TickDirection','out');
set(hcb,'TickLabelInterpreter','latex');
set(hcb,'FontSize',14);
set(hcb,'TickLabels',ticklabels);
set(hcb,'Ticks',ticks);

ind_dir = zeros(size(alpha_all))-4;
for ii = 1:numel(alpha_all)
    if isnan(th_MCE_ortho(ii)) && isnan(th_MCE_anti(ii))
        ind_dir(ii) = 10;
    elseif isnan(th_MCE_ortho(ii))
        ind_dir(ii) = 2;
    elseif  abs(th_MCE_ortho(ii)/th_MCE_anti(ii)-1)<1e-2
        ind_dir(ii) = 0;
    elseif abs(th_MCE_ortho(ii))*1.005<abs(th_MCE_anti(ii)*0.995)
        ind_dir(ii) = 1;
    else
        ind_dir(ii) = 2;
    end
end

% ind_dir(396) = 0;
% ind_dir(413) = 0;
% % ind_dir([390,391,411,412])
% ind_dir([307:21:370,390,411]) = 1;
ind_dir(alpha_all>=(110*pi/180) & PW_all>0.00) = 2;
% ind_dir([391,412]) = 2;
% 
ind_2 = (ind_dir ==2);
ind_dir(ind_2) =0;
contour(alpha_all/pi*180,log10(PW_all),ind_dir,[0.5,0.5],'--','LineWidth',2,'Color',[1,1,1]*0.8);
% ind = find(alpha_all(:,1)<=(111*pi/180));
% ind_PW =find(PW_all(1,:)<=1.5);
% contour(alpha_all(ind ,ind_PW)/pi*180,log10(PW_all(ind ,ind_PW)),ind_dir(ind ,ind_PW),[0.5,0.5],'-','LineWidth',2,'Color',[1,1,1]*1);
% ind = find(alpha_all(:,1)<=(100*pi/180));
% ind_PW =find(PW_all(1,:)>0.9);
% contour(alpha_all(ind ,ind_PW)/pi*180,log10(PW_all(ind ,ind_PW)),ind_dir(ind ,ind_PW),[0.5,1.4],'-','LineWidth',2,'Color',[1,1,1]*1);
% ind = find(alpha_all(:,1)<=(110*pi/180) & alpha_all(:,1)>=(104*pi/180));
% ind_PW =find(PW_all(1,:)>=1.1);
% contour(alpha_all(ind ,ind_PW)/pi*180,log10(PW_all(ind ,ind_PW)),ind_dir(ind ,ind_PW),[1.2,1.2],'-','LineWidth',2,'Color',[1,1,1]*1);
% plot([106,106],log10([1,10]),'-','LineWidth',2,'Color',[1,1,1]*1);
ind_dir(ind_2) =2;
ind_dir(~ind_2) =1;
contour(alpha_all/pi*180,log10(PW_all),ind_dir,[1.5,1.5],'--','LineWidth',2,'Color',[1,1,1]);

% ind_PW =find(PW_all(1,:)>0.7);
% contour(alpha_all(:,ind_PW)/pi*180,log10(PW_all(:,ind_PW)),ind_dir(:,ind_PW),[0.5,0.55],'-','LineWidth',2,'Color',[1,1,1]*0.8);

% return
%%
axes('position',[0.7,0.1,0.25,0.8]);
box on; hold on;

% cm = jet(128);
% cm(1:2,:)=cm(1:2,:).*kron((0.6:0.3:0.95)',[1,1,1]);
colormap(gca, cm2);
caxis(gca,-fliplr([-100,0]));

c_lvl = -100:.1:0;
% surf(log10(H_all/1e-1),YY_all/1e-1,H_all*0,th_per_diff_MCE,'LineStyle','none');
contourf(alpha_all/pi*180,log10(PW_all),-th_per_diff_MCE,-c_lvl,'LineStyle','none');



plot([45,45],[-3,1],'k:','LineWidth',1.5);
plot([90,90],[-3,1],'k:','LineWidth',1.5);
% plot([109,109],[-3,1],'w--','LineWidth',2);
% plot([116,116],[-3,1],'w--','LineWidth',2);

ind = find(PW_all(1,:)>=(.01));
[C,h]=contour(alpha_all(:,ind)/pi*180,log10(PW_all(:,ind)),th_per_diff_MCE(:,ind),[-60,-15,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
% clabel(C,h,'FontSize',10,'FontWeight','Bold');
% ind = find(PW_all(1,:)<=(0.1) & PW_all(1,:)>=(0.01));
% [C,h]=contour(alpha_all(:,ind)/pi*180,log10(PW_all(:,ind)),th_per_diff_MCE(:,ind),[-15,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
ind = find(PW_all(1,:)<=(0.01));
[C,h]=contour(alpha_all(:,ind)/pi*180,log10(PW_all(:,ind)),th_per_diff_MCE(:,ind),[-60,-15,-5],'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',10,'FontWeight','Bold');

set(gca,'XTick', xtick);
set(gca,'XTickLabel',xticklabel);
set(gca,'YTick',ytick);
set(gca,'YTickLabel',yticklabel);

axis([0,135,-3,1]); %axis equal;

xlabel({'Axon''s angle with E-field $$ \alpha$$'},'Interpreter','latex'); 
% ylabel('Axon vertical placement $$Y \: \mathrm{(msm)}$$','Interpreter','latex');
title({'','Percentage difference of threshold: \textbf{IP} vs. \textbf{1D}'},'Interpreter','latex');
set(get(gca,'XLabel'), journal_axis_label);
set(get(gca,'YLabel'),journal_axis_label);
set(get(gca,'Title'), journal_title);
set(gca,journal_axis);

hcb=colorbar(gca,'SouthOutside'); 
set(hcb,'TickDirection','out');
set(hcb,'TickLabelInterpreter','latex');
set(hcb,'FontSize',14);
ticks=-fliplr([-100:25:0]);
ticklabels=cell(size(ticks));
ticklabels{end}='$$-100\%$$';
for ii=1:length(ticks)-1
    ticklabels{ii}=sprintf('$$%1.1f  \\%%$$',-ticks(ii));
end
set(hcb,'Ticks',ticks);
set(hcb,'TickLabels',ticklabels);

ind_dir(ind_2) =2;
ind_dir(~ind_2) =1;
contour(alpha_all/pi*180,log10(PW_all),ind_dir,[1.5,1.5],'r--','LineWidth',2);
plot([116,116],[-3,1],'r--','LineWidth',2);
% return

model_name='UF_Bnd_HH';
filename = fullfile([model_name,'_compiled_result']);
saveas(gca,[filename,'.fig']);
frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(    imind,cm,[filename,'.tif'],...
            'tif','WriteMode','overwrite', 'Resolution',300);
  