function col = TmazeColors(colormode)
% col = TMAZECOLORS(colormode)
% Generates a struct with rgb values for plotting Tmaze figures. This is a
% way of easily saving and trying out different colormodes.
%
% output has the following organization
%       col.all.f = [r g b]; % or 'r' depending on the mode
%       col.R042.f = [r g b]; 
%          ~~~
%       col.R064.w = [r g b];
%
%       f is food color, w is water color
%       f color is always darker than w color
%
% colormode = 
%    'inventory1'; as in Tmaze data inventory, with combined rats in grey
%    'inventory2'; as in inventory, with combined rats in red and blue
%    'inventory3'; " " with red and blue less saturated
%    'inventory4'; as with inventory 3, but with individual rat colors less
%    saturated
%    'single'; coordinating colors when not using food AND water data in
%    same plot; in this case, use only the "f" field ex: col.R042.f, but not
%    col.R042.w. For all data, col.all.w = 'k' and col.all.f = 'w'.
%    'rb'; food is always red and water is always blue, regardless of rat
%    'grey'; food is dark grey and water is light grey
%
% aacarey Sept 2015

%% Example usage with dynamic struct fields during plotting:
% 
% rats = {'R042','R044','R050','R064'};
% for iRat = 1:length(rats)
%    colors = TmazeColors('inventory2');
%    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
%    d = [val1 val2 val3 val4]; % bar graph data values
%    %*** plot stuff here **
%    for iBar = 1:length(d)
%        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
%        hold on
%    end
% end

switch colormode
    case 'inventory1' % with grey for combined data
        
        col.R042.f = [0 100/255 0];  %dark green
        col.R042.w = [143/255 188/255 143/255];  % lighter green
        col.R044.f = [54/255 100/255 139/255];  %dark blue
        col.R044.w = [176/255 196/255 222/255];  % lighter blue
        col.R050.f = [104/255 34/255 139/255];  %dark purple
        col.R050.w = [216/255 191/255 216/255];  % lighter purple
        col.R064.f = [198/255 113/255 113/255];  %dark salmon
        col.R064.w = [238/255 213/255 210/255];  % lighter salmon/dusty pink
        col.all.f =  [105/255 105/255 105/255];  %dark grey
        col.all.w =  [211/255 211/255 211/255];  % lighter grey
        
    case 'inventory2' % with rb for combined data
        
        col.R042.f = [0 100/255 0];  %dark green
        col.R042.w = [143/255 188/255 143/255];  % lighter green
        col.R044.f = [54/255 100/255 139/255];  %dark blue
        col.R044.w = [176/255 196/255 222/255];  % lighter blue
        col.R050.f = [104/255 34/255 139/255];  %dark purple
        col.R050.w = [216/255 191/255 216/255];  % lighter purple
        col.R064.f = [198/255 113/255 113/255];  %dark salmon
        col.R064.w = [238/255 213/255 210/255];  % lighter salmon/dusty pink
        col.all.f =  'r';
        col.all.w =  'b';
        
    case 'inventory3' % with rb for combined data
        
        col.R042.f = [0 100/255 0];  %dark green
        col.R042.w = [143/255 188/255 143/255];  % lighter green
        col.R044.f = [54/255 100/255 139/255];  %dark blue
        col.R044.w = [176/255 196/255 222/255];  % lighter blue
        col.R050.f = [104/255 34/255 139/255];  %dark purple
        col.R050.w = [216/255 191/255 216/255];  % lighter purple
        col.R064.f = [198/255 113/255 113/255];  %dark salmon
        col.R064.w = [238/255 213/255 210/255];  % lighter salmon/dusty pink
        col.all.f =  [176/255 23/255 31/255]; % deep red
        col.all.w =  [39/255 64/255 139/255]; % dark blue
        
    case 'inventory4' % with rb for combined data
        
        col.R042.f = [92/255 121/255 75/255];
        col.R042.w = [191/255 213/255 178/255];
        col.R044.f = [89/255 101/255 135/255];
        col.R044.w = [194/255 198/255 210/255];
        col.R050.f = [106/255 75/255 121/255];
        col.R050.w = [197/255 181/255 205/255];
        col.R064.f = [198/255 113/255 113/255];
        col.R064.w = [238/255 213/255 210/255];
        col.all.f =  [176/255 23/255 31/255];
        col.all.w =  [39/255 64/255 139/255];
      
    case 'rb'
        col.R042.f = 'r';
        col.R042.w = 'b';
        col.R044.f = 'r';
        col.R044.w = 'b';
        col.R050.f = 'r';
        col.R050.w = 'b';
        col.R064.f = 'r';
        col.R064.w = 'b';
        col.all.f =  'r';
        col.all.w =  'b';
        
    case 'grey'
        col.R042.f = [105/255 105/255 105/255];  %dark grey
        col.R042.w = [211/255 211/255 211/255];  % lighter grey
        col.R044.f = [105/255 105/255 105/255]; 
        col.R044.w = [211/255 211/255 211/255]; 
        col.R050.f = [105/255 105/255 105/255];  
        col.R050.w = [211/255 211/255 211/255]; 
        col.R064.f = [105/255 105/255 105/255]; 
        col.R064.w = [211/255 211/255 211/255]; 
        col.all.f =  [105/255 105/255 105/255];  
        col.all.w =  [211/255 211/255 211/255];  
        
    otherwise
        error('Unrecognized colormode. Better check that spelling ^_^')
end

