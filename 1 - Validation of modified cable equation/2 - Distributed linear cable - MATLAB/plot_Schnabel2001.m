% Schnabel data:        R = 3 um; PW = max (100 ms)

phi_m_ratio = zeros(length(results),1);    % Metric used in Schnabel & Struijk 2001
for ii = 1:length(results)
    phi_m_ratio(ii)= 2 * results(ii).ExR / results(ii).phi_m_bar(end) +1;
end

figure;
h_line = plot(  log10([results.H]/1e-1),phi_m_ratio,...    % semi-log x
                '-.ok','MarkerSize',6,'LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','w');

%% Figure formatting
set(gcf,'Position',[50,50,600*1.5,600],format_figure);
set(gca,format_axis);
axis([-3,3,0,4]);
set(gca,'DataAspectRatio',[1,1.27*4/6,1]);                          % Aspect ratio adjusted to fit Fig. 2 of Schnabel & Struijk 2001

set(get(gca,'Title'), format_title);
set(get(gca,'XLabel'),format_axis_label);
set(get(gca,'YLabel'),format_axis_label);
title( 'Ratio of maximum depolarization at steady state');
xlabel('Axon-electrode distance $$ h \: \rm{ (mm)}$$');
ylabel('$$\mathrm{max}_{z,\theta}(\varphi_{\mathrm{m}(z,\theta)}) / \mathrm{max}_{z}(\overline{\varphi}_{\mathrm{m}(z)})$$');

% Change x-axis labels to log scale
xticks = get(gca,'XTick');
xticklabels = cell(size(xticks));
for ii = 1: length(xticks)
   xticklabels{ii} = ['$$10^{',num2str(xticks(ii),'%d'),'}$$']; 
end
set(gca,'XTickLabel',xticklabels);

figurename = fullfile('Figures',[electrode_str{electrode_type},'_A_Lin_Schnabel']);
saveas(gcf,[figurename,'.fig']);
[imind,cm] = rgb2ind(frame2im(getframe(gcf)),256);
imwrite(    imind,cm,[figurename,'.tif'],'tif','WriteMode','overwrite', 'Resolution',600,'Compression','none');

if plot_trans
    % Only retain curve and box for overlaying on Fig.2 of Schnabel & Struijk 2001
    set(gca,'XTickLabel',{},'YTickLabel',{},'TickDir','out','XColor','b','YColor','b');
    set(h_line,'Color','b','MarkerEdgeColor','b');
    xlabel('');ylabel('');title('');
    
    [imind,cm] = rgb2ind(frame2im(getframe(gcf)),256);
    imwrite(    imind,cm,[figurename,'_outline.tif'],'tif','WriteMode','overwrite', 'Resolution',150,'Compression','none');
end

close(gcf);