%% Enable scroll and zoom
function enscroll_nosepokes
set(gcf,'KeyPressFcn',@enscroll_nosepokes);

global event_number_nosepokes
global dataPoint_nosepokes
global tbin_centers
global h2
function enscroll_nosepokes(~,event)
if strcmp(event.Key,'rightarrow');    
    set(gca,'XLim',[tbin_centers(dataPoint_nosepokes(event_number_nosepokes+1)-5000) tbin_centers(dataPoint_nosepokes(event_number_nosepokes+1)+10000)]);    
%     if exist('h2', 'var')
%   delete(h2)
%     end
    hold on; h2 = plot(tbin_centers(dataPoint_nosepokes(event_number_nosepokes+1)),0:1:70,'color','black');
    event_number_nosepokes = event_number_nosepokes+1;
elseif strcmp(event.Key,'leftarrow');    
    set(gca,'XLim',[tbin_centers(dataPoint_nosepokes(event_number_nosepokes-1)-5000) tbin_centers(dataPoint_nosepokes(event_number_nosepokes-1)+10000)]);
%     if exist('h2', 'var')
%   delete(h2)
%     end
    hold on; h2 = plot(tbin_centers(dataPoint_nosepokes(event_number_nosepokes-1)),0:1:70,'color','black');
    event_number_nosepokes = event_number_nosepokes-1;
elseif strcmp(event.Key,'uparrow');
    zoom(2);
elseif strcmp(event.Key,'downarrow');
    zoom(.5);
end
end
end