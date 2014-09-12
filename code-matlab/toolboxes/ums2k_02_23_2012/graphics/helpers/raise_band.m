function raise_band (hObject, event, hndls)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% raise_band - callback to raise clicked patch to front of visibility
%
% Usage:
%      raise_band (hObject, event, me)
%
% Description:  
% Pulls the clicked object to the top of the ui stack; useful
% for raising partially masked objects to the front of a plot.
% GUI-shortcut for 'uistack':  Left-clicking brings to the top,
% right-clicking sends to bottom.  This is accomplished by finding the
% highest (lowest) z-values in the figure and setting the current objects'
% z-value to one higher (lower).
%
% Input: 
%   hObject - handle to object
%   event  -  event identifier (ignored)
%   hndls  - handles to the "band"
%

% get top z position
[o,h] = gcbo;
pos = get( findobj(h,'Type','patch'), 'ZData' );

% only do something if pos is a cell and therefore there are many objects
if iscell(pos)
    
  % find the minimum and maximum 
  extrema = minmax(cell2mat(pos) );
   
   % normal click results in raising the pathces to make them visible
    switch(get(h, 'SelectionType')),
        case ('normal'),
            uistack(hndls, 'top');
            for j = 1:length(hndls)
              set( hndls(j),'Zdata', (extrema(2)+1)*get_z_ones( hndls(j) )  );
            end

      % right click results in lowering the patches behind other objects
        case ('alt'),            
            uistack(hndls, 'bottom');
            for j = 1:length(hndls)
              set( hndls(j),'Zdata', (extrema(1)-1)*get_z_ones( hndls(j) )  );
            end
    end

    drawnow;
end

% helper function to generate array of ones of the approiate size
function b = get_z_ones( hndl )
            b = ones( size( get(hndl,'ZData') ) );
    