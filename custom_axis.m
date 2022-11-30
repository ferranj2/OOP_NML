function ax = custom_axis(f)
if nargin == 0
    f = figure(...
        'Units','pixels',...
        'Position',[34 127.5000 560 420]);
end%if
ax = axes('Parent',f);
hold(ax,'on');
ax.BoxStyle = 'full';
ax.Box = 'on';
ax.NextPlot = 'Add'; %Normies use the "hold on" command.
%ax.XLim = [-1,+1];
%ax.YLim = [-1,+1];
axis(ax,'equal');

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

%Z-labels
ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.FontName = 'TimesNewRoman';
ax.ZAxis.TickLabelInterpreter = 'latex';
ax.ZAxis.Label.String = 'Z'; %Default string to display.

%Grid-related.
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';

end