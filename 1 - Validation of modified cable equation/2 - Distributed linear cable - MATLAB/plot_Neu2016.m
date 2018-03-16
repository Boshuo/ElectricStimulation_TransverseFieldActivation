H = [results.H];
t_ss = [results.t_ss];
per_diff=zeros(length(PW),length(H));

for ii = 1 : length(H)  % Metric used in Neu 2016
    per_diff(:,ii)= - (2 * results(ii).ExR ) ./ (results(ii).phi_m_bar + 2 * results(ii).ExR ) *100;
end

figure;
hold on;

% Breaking countour lines into smaller regions to control the location of
% contour labels
h_ind = find(H<0.1);
contour(log10(H(h_ind)/1e-1),log10(PW),per_diff(:,h_ind),-[0.05,0.1:0.1:0.9,0.95,0.99,0.999]*100,'-','LineWidth',1.5,'Color',[1,1,1]*0);
h_ind = find(H<0.4 & H>0.05);
[C,h]=contour(log10(H(h_ind)/1e-1),log10(PW),per_diff(:,h_ind),-[0.05,0.1:0.1:0.9,0.95,0.99]*100,'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','Labelspacing',400);
h_ind = find(H>0.3);
contour(log10(H(h_ind)/1e-1),log10(PW),per_diff(:,h_ind),-[0.05,0.1:0.1:0.9,0.95,0.99]*100,'-','LineWidth',1.5,'Color',[1,1,1]*0);
[C,h]=contour(log10(H(h_ind)/1e-1),log10(PW),per_diff(:,h_ind),[-0.999,-0.999]*100,'-','LineWidth',1.5,'Color',[1,1,1]*0);
clabel(C,h,'FontSize',13,'Interpreter','latex','Labelspacing',100);

plot(log10(H/1e-1),ones(size(H))*log10(mean(t_ss)),'k--','LineWidth',2);    % Average of steady-state time

% Plot shaded regions
colormap(gray(64));
caxis([0,1]);
S.vertices = [log10([4.5,4.5,20,29,29,20]')-3,[-4,2,2,2,-4,-4]'];
S.faces = (1:6);
S.FaceVertexCData = [0.8,0.8,0.9,1,1,0.9]'-0.3;
S.FaceColor= 'interp'; S.EdgeColor = 'none';
S.LineStyle = 'none'; S.FaceAlpha  = 0.5;
patch(S);

S.vertices = [log10([31e-3*[10,10^0.25,1],1*[10,10^0.75,1]]'),[0,0.75,1,1,0.75,0]'];
S.FaceVertexCData = [1,0.9,0.8,0.8,0.9,1]'-0.25;
S.FaceAlpha  = 0.3;
patch(S);
S.vertices = [log10([1,10^0.75,10,10,10^0.75,1]'),[-1,-1.75,-2,1,0.75,0]'];
patch(S);
S.vertices = [log10([[1,10^0.75,10],31e-3*[1,10^0.25,10]]'),[-1,-1.75,-2,-2,-1.75,-1]'];
patch(S);
S.vertices = [log10(31e-3*[10,10^0.25,1,1,10^0.25,10]'),[0,0.75,1,-2,-1.75,-1]'];
patch(S);
S.vertices = [log10([31e-3*[10,10],1*[1,1]]'),[-1,0,0,-1]'];
S.faces = (1:4);
S.FaceVertexCData = [1,1,1,1]'-0.25;
patch(S);

%% Figure formatting
set(gcf,'Position',[50 50 750 600],'Color',format_figure.Color);
set(gca, format_axis);
axis(log10([4.5e-3,100,1e-4,100]));

set(get(gca,'XLabel'), format_axis_label);
set(get(gca,'YLabel'), format_axis_label);
set(get(gca,'Title'), format_title);
title({'Percentage difference of maximum depolarization: ','\textbf{Conventional CE} vs. \textbf{modified CE}'});
xlabel('Axon \textendash electrode distance $$ H \: \rm{ (mm)}$$');
ylabel('Pulse width $$PW \: \rm{(ms)}$$');

xtick = log10([1e-3*[5,6,7,8,9],kron([1e-2,1e-1,1,10],[1,2,3,4,5,6,7,8,9]),100]);
set(gca,'XTick', xtick);
xticklabel = cell(length(xtick),1);
for ii = 1 : length(xticklabel)
    if abs(round(xtick(ii))-xtick(ii)) < 1e-2
        xticklabel{ii} = ['$$10^{',num2str(xtick(ii),'%d'),'}$$'];
    else
        xticklabel{ii} = '';
    end
end
set(gca,'XTickLabel',xticklabel);

ytick = log10([kron([1e-4,1e-3,1e-2,1e-1,1,10],[1,2,3,4,5,6,7,8,9]),100]);
set(gca,'YTick',ytick);
yticklabel=cell(length(ytick),1);
for ii = 1 : length(yticklabel)
    if abs(round((ytick(ii)))-(ytick(ii))) < 1e-6
        yticklabel{ii} = ['$$10^{',num2str((ytick(ii)),'%d'),'}$$'];
    else
        yticklabel{ii} = '';
    end
end
set(gca,'YTickLabel',yticklabel);

set(gca,format_axis);
figurename = fullfile('Figures',[electrode_str{electrode_type},'_A_Lin_Neu']);
saveas(gcf,[figurename,'.fig']);

[imind,cm] = rgb2ind(frame2im(getframe(gcf)),256);
imwrite(    imind,cm, [figurename,'.tif'],'tif','WriteMode','overwrite', 'Resolution',600,'Compression','none');

if plot_trans
    % Only retain curve and box for overlaying on Fig.3 of Neu 2016
    set(gca,'DataAspectRatio',[1,1.156*6/(5-log10(4.5)),1]);
    xtick = log10([1e-2,1e-1,1,10,100]);
    set(gca,'XTick', xtick);
    ytick = log10([1e-4,1e-3,1e-2,1e-1,1,10,100]);
    set(gca,'YTick',ytick);
    
    xlabel('');ylabel('');title('');
    delete(findobj(get(gca,'Children'),'Type','Patch'));
    delete(findobj(get(gca,'Children'),'Type','Line'));
    set(gca,'XTickLabel',{},'YTickLabel',{},'TickDir','out','XColor','r','YColor','r');
    set(get(gca,'Children'),'Color','r','LineWidth',1.5','LineStyle','--','Showtext','off')
    [imind,cm] = rgb2ind(frame2im(getframe(gcf)),256);
    imwrite(    imind,cm,[figurename,'_outline.tif'],'tif','WriteMode','overwrite', 'Resolution',150,'Compression','none');
end

close(gcf);