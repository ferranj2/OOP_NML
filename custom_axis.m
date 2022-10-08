function ax = custom_axis(f)
if nargin == 0
    f = figure;
end
ax = axes('Parent',f);
hold(ax,'on');
ax.BoxStyle = 'full';
ax.Box = 'on';
ax.NextPlot = 'Add'; %Normies use the "hold on" command.
%ax.XLim = [-1,+1];
%ax.YLim = [-1,+1];
axis(ax,'square');

%X-labels
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontName = 'TimesNewRoman';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.XAxis.Label.String = 'X'; %Default string to display.

%Y-Labels
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontName = 'TimesNewRoman';
ax.YAxis.TickLabelInterpreter = 'latex';
ax.YAxis.Label.String = 'Y'; %Default string to display.

%Grid-related.
ax.XGrid = 'on';
ax.YGrid = 'on';

end