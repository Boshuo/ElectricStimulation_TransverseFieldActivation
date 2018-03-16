%% figure format
format_axis.LineWidth = 1;
format_axis.FontSize = 16;
% format_axis.FontName = 'Times New Roman';
format_axis.Color = [1 1 1]; % background color of the plot area
format_axis.TickLabelInterpreter = 'latex';
format_axis.TickDir = 'out';

format_axis_label.FontSize = 18;
format_axis_label.FontName = 'Times New Roman';
format_axis_label.Interpreter = 'latex';

format_figure.Color = [1 1 1]; % background color of figure window

format_trace.LineWidth = 1.5;

format_title.FontSize = 20;
format_title.Interpreter = 'latex';

%% Color map
cm0 = flipud(parula(256));
cm1 = cm0;
cm2 = parula(256);

cm0(end-8:end,:) = cm0(end-8:end,:).*kron((0.95:-0.1:0.15)',[1,1,1]);
for ii = 256:-1:129
    cm2(ii,:)= brighten(cm2(ii,:),(ii-129)/128);
end

cb_title.Interpreter = 'latex';
cb_title.FontSize = 16;
cb_title.FontWeight ='Normal';
cb_title.Units = 'Normalized';
cb_title.HorizontalAlignment = 'Center';
cb_title.VerticalAlignment = 'Bottom';
cb_title.LineWidth = 1;
cb_title.Color = 'k';

format_cb.TickDirection = 'out';
format_cb.TickLength = 0.02;
format_cb.TickLabelInterpreter = 'latex';
format_cb.FontSize = 16;
format_cb.AxisLocation = 'out';

%%

ytick = log10([kron([1e-3,1e-2,1e-1,1],[1,2,3,4,5,6,7,8,9]),10]);
yticklabel = cell(length(ytick),1);
for ii = 1 : length(yticklabel)
    if abs(round(ytick(ii)) - ytick(ii)) < 1e-2
        yticklabel{ii} = ['$$10^{',num2str((ytick(ii)),'%d'),'}$$'];
    else
        yticklabel{ii} =' ';
    end
end