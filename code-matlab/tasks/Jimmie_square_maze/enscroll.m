%% Enable scroll and zoom
function enscroll
set(gcf,'KeyPressFcn',@enscroll);

global event_number
global dataPoint
global tbin_centers
global h1
function enscroll(~,event)
if strcmp(event.Key,'rightarrow');    
    set(gca,'XLim',[tbin_centers(dataPoint(event_number+1)-5000) tbin_centers(dataPoint(event_number+1)+10000)]);    
%     if exist('h1', 'var')
%   delete(h1)
%     end
    hold on; h1 = plot(tbin_centers(dataPoint(event_number+1)),0:1:70,'color','black');
    event_number = event_number+1;
elseif strcmp(event.Key,'leftarrow');    
    set(gca,'XLim',[tbin_centers(dataPoint(event_number-1)-5000) tbin_centers(dataPoint(event_number-1)+10000)]);
%     if exist('h1', 'var')
%   delete(h1)
%     end
    hold on; h1 = plot(tbin_centers(dataPoint(event_number-1)),0:1:70,'color','black');
    event_number = event_number-1;
elseif strcmp(event.Key,'uparrow');
    zoom(2);
elseif strcmp(event.Key,'downarrow');
    zoom(.5);
end
end
end