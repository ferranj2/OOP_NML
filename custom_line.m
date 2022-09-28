function L = custom_line(varargin)
L = line;
%Default presets.
L.DisplayName = '?';
L.Color = [0,0,0];
L.LineStyle = '-';
L.LineWidth = 1;

L.Marker = 'none';
L.MarkerEdgeColor = [0,0,0];
L.MarkerFaceColor = [0,0,1];

if nargin == 0
    return;
end%if

 ii = 0;
 valid_input = false;
 while ii < nargin
     ii = ii + 1;
     iip1 = ii + 1;
     %Assign XY Data in bulk.
     if strcmpi(varargin{ii},'Parent') == 1
         L.Parent = varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'XY') == 1
         L.XData = varargin{iip1}(:,1);
         L.YData = varargin{iip1}(:,2);
         valid_input = true;
     end%if
     %Assign only X-Data
     if strcmpi(varargin{ii},'X') == 1
         L.XData = varargin{iip1};
         valid_input = true;
     end%if
     %Assign only Y-Data
     if strcmpi(varargin{ii},'Y') == 1
         L.YData = varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'Z') == 1
         L.ZData = varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'DisplayName') == 1
         L.DisplayName = varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'Color') == 1
         L.Color = varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'Marker') == 1
         L.Marker= varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'MarkerEdgeColor') == 1
         L.MarkerEdgeColor= varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'MarkerFaceColor') == 1
         L.MarkerFaceColor= varargin{iip1};
         valid_input = true;
     end%if
     if strcmpi(varargin{ii},'LineStyle') == 1
         L.LineStyle = varargin{iip1};
         valid_input = true;
     end%if
     
     if ~valid_input
         error('One of the inputs is unrecognizable!');
     end%if
     ii = iip1;
     valid_input = false;
 end%while


end%function