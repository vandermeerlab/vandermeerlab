function navigate(src,event)
% NAVIGATE Dynamically alter figure xlimits using keyboard and user input
% 
% Navigate is a callback/KeyPress function. Navigate requires that some variables 
% be set to global. Note that Navigate will overwrite figure titles during
% event navigation.
%
% Navigate's global variables:
%    evtTimes: the center times of events/iv (for event navigation)
%    time: the time vector 
%    windowSize: the desired size of the viewing window in time vector units
%
% See MultiRaster for basic usage.
%
% If I have my own keypress function that does other things, how can I use 
% navigate at the same time?  
%  -- Global variables must be properly assigned, then make your callback
%   function myKeyPressFcn(source,event)
%       if strcmp(event.Key,'return')
%           % do some stuff here
%       else
%           navigate(source,event)
%       end
%   end
%
% KEYBOARD COMMANDS:
%
%    Help
%    h: display a message box with all keyboard commands (please update
%    commandList if new features added)
%
%    Zoom 
%    w:              user inputs a new windowSize (seconds)
%    Crtl+equal:     zoom in (aka control plus); resets windowSize (/2) 
%    Ctrl+hyphen:    zoom out (aka control minus); resets windowSize (*2) 
%    Ctrl+backslash: view entire session; does not reset windowSize 
%
%    Event-Based Navigation
%    a:       1 event left
%    d:       1 event right
%    Shift+a: 10 events left
%    Shift+d: 10 events right
%    Ctrl+a:  50 events left
%    Ctrl+d:  50 events right
%    f:       first event
%    c:       center event
%    l:       last event
%
%    Move Viewing Window
%    j:           user inputs a timepoint to "jump" to (seconds)
%    leftarrow:   0.5 windows left
%    rightarrow:  0.5 windows right
%    Shift+left:  10 windows left
%    Shift+right: 10 windows right
%    Ctrl+left:   50 windows left
%    Ctrl+right:  50 windows right
%    b:           beginning of the recording
%    m:           middle of the recording
%    e:           end of the recording
%    
% aacarey. Aug 2014. 
% youkitan edit, 2015-01-20 (can display candidate score).
% aacarey edit, Feb 2015 (added j and w commands), Mar 2015 (added zoom and h)
% aacarey edit, Nov 2015 (fixed problems/"shortcomings" with event navigation)

%% ****Some notes on how this function works internally****

% Based on the current window limits and user keyboard input, this
% function generates new limits to change the view. Depending on the type 
% of movement asked for by the user, the function's if statements generate 
% different types of output:

% TYPE 1 if statement outputs are already limits and thus go directly into the
% "set(gca,'XLim',limits)" at the very end of the function. (Outputs cannot 
% enter the "permission" if statement.)

% TYPE 2 if statements must enter the "permission" if statement, which
% converts their outputs into new limits before sending the new output to set(gca,'XLim',limits)
% at the end of the function code. The function was organized
% internally in this way to reduce the amount of repetition in the code.

% There is likely a better way to do this, but I am a novice. - ACarey
%% Declare global variables

global evtTimes windowSize time usrfield


%%
limits = get(gca,'XLim'); %get current limits; exists as limits(1) and limits(2), in array [1 2]
current_location = mean(limits); % finds center of current viewing window

modifier = get(gcf,'currentmodifier'); %checks for modifiers, such as Shift or Ctrl

shiftPressed = ismember('shift',modifier);
ctrlPressed = ismember('control',modifier);
%altPressed = ismember('alt',modifier); % if alt is pressed while on a figure, the focus is moved elsewhere, so the figure has to be clicked on again

% Here's an alternative way to detect modifier keys:
%strcmp(event.Modifier{:},'control')

flank = 0.5*windowSize;

permission = 1; %this controls whether or not the outputs from the if statements can enter the  
% "permission" if statements, or whether they go directly towards changing
% the limits at the very end of the function

% evtTitle variable: this controls whether we display an event title, which 
% occurs only when we use the "find event" functionality of navigate;
% "move window" changes the current view regardless of the presence of an
% event, so we don't want to display an event title if an event is not the
% central focus of the viewing window! 

%% HELP

if strcmp(event.Key,'h')
commandList = ...
{'Make sure figure window is in focus (click window)';
' ';
'    Zoom'; 
'    w:              user inputs a new windowSize (seconds)';
'    Crtl+equal:     zoom in (aka control plus); resets windowSize (/2)'; 
'    Ctrl+hyphen:    zoom out (aka control minus); resets windowSize (*2) ';
'    Ctrl+backslash: view entire session; does not reset windowSize'; 
' ';
'    Event-Based Navigation';
'    a:       1 event left';
'    d:       1 event right';
'    Shift+a: 10 events left';
'    Shift+d: 10 events right';
'    Ctrl+a:  50 events left';
'    Ctrl+d:  50 events right';
'    f:       first event';
'    c:       center event';
'    l:       last event';
' '
'    Move Viewing Window';
'    j:           user inputs a timepoint to "jump" to (seconds)';
'    leftarrow:   0.5 windows left';
'    rightarrow:  0.5 windows right';
'    Shift+left:  10 windows left';
'    Shift+right: 10 windows right';
'    Ctrl+left:   50 windows left';
'    Ctrl+right:  50 windows right';
'    b:           beginning of the recording';
'    m:           middle of the recording';
'    e:           end of the recording'};

msgbox(commandList,'How to Navigate');
end
%% CHANGE WINDOW SIZE (ZOOM)

if strcmp(event.Key,'w') 
    win = inputdlg('New window size = ____ seconds');
    if ~isempty(win) && ~isnan(str2double(win{1})) && isscalar(str2double(win{1}))
        windowSize = str2double(win{1});
        flank = 0.5*windowSize;
        limits = [current_location-flank current_location+flank];
    end
end

if ctrlPressed && strcmp(event.Key,'backslash')
    limits = [time(1)-200 time(end)+200]; % adding the constant is for margins
    % note: this doesn't actually reset windowSize, and it's probably better this way 
end

if ctrlPressed && strcmp(event.Key,'equal') % then zoom in
    windowSize = windowSize./2;
    flank = 0.5*windowSize;
    limits = [current_location-flank current_location+flank];
end

if ctrlPressed && strcmp(event.Key,'hyphen') % then zoom out
    windowSize = windowSize.*2;
    flank = 0.5*windowSize;
    limits = [current_location-flank current_location+flank];
end

%disp(event.Key)

%% JUMP TO TIME INPUT BY USER

if strcmp(event.Key,'j')
    jumpHere = inputdlg('Jump to time = ___ seconds');
    if ~isempty(jumpHere) && ~isnan(str2double(jumpHere{1})) && isscalar(str2double(jumpHere{1}))
        next_location = str2double(jumpHere{1});
        limits = [next_location-flank next_location+flank];
    end 
end


%% MOVE WINDOW
window = (limits(2)-limits(1));
       %type 1: output is a new set of limits
        if shiftPressed == 0 && ctrlPressed == 0 && strcmp(event.Key,'leftarrow')==1 %string comparison
                %shift figure axes left
                shift = window*0.5;
                limits = limits - shift;
                permission = 0;
                evtTitle = 0;
         
        %type 1: output is a new set of limits        
        elseif shiftPressed == 0 && ctrlPressed == 0 && strcmp(event.Key,'rightarrow')==1 
                %shift figure axes right
                shift = window*0.5;
                limits = limits + shift;
                permission = 0;
                evtTitle = 0;
         
        %type 1: output is a new set of limits        
        elseif shiftPressed == 1 && strcmp(event.Key,'leftarrow') == 1
            % shift 10 windows left
                shift = window*10;
                limits = limits - shift;
                permission = 0;
                evtTitle = 0;
         
        %type 1: output is a new set of limits        
        elseif shiftPressed == 1 && strcmp(event.Key,'rightarrow') == 1
            % shift 10 windows right
                shift = window*10;
                limits = limits + shift;
                permission = 0;
                evtTitle = 0;
        
        %type 1: output is a new set of limits
        elseif ctrlPressed == 1 && strcmp(event.Key,'leftarrow') == 1
            % shift 50 windows left
                shift = window*50;
                limits = limits - shift;
                permission = 0;
                evtTitle = 0;
        
        %type 1: output is a new set of limits        
        elseif ctrlPressed == 1 && strcmp(event.Key,'rightarrow') == 1
            % shift 50 windows right
                shift = window*50;
                limits = limits + shift;
                permission = 0;
                evtTitle = 0;
         
        %type 2: output is a next_location, which later is converted to new limits        
        elseif strcmp(event.Key,'m')==1  
                %shift figure axes to the midpoint of recording
                next_location = time(ceil(end/2));
                limits = [next_location-flank next_location+flank];
                evtTitle = 0;
          
        %type 2: output is a next_location, which later is converted to new limits        
        elseif strcmp(event.Key,'b')==1
            %shift figure axes to beginning of recording
                next_location = time(1)+flank;
                limits = [next_location-flank next_location+flank];
                evtTitle = 0;

        %type 2: output is a next_location, which later is converted to new limits        
        elseif strcmp(event.Key,'e') == 1
            %shift figure axes to end of recording
                next_location = time(end)-flank;
                limits = [next_location-flank next_location+flank];  
                evtTitle = 0;
        else
            evtTitle = 1; % if a non-assigned key was pressed, don't set the title to blank
        end
        
% Do not show a title if the user is navigating based on the time axis
% (i.e. "move window" key commands) rather than based on event locations
% (i.e. "find event" key commands).
    
    if evtTitle == 0 % evtTitle will be blank if you are using the move window keys rather than the find event keys
        title(sprintf('')); 
    end
        
%% FIND EVENT
current_index = nearest_idx3(current_location,evtTimes); % find where we're located with respect to time. Take that index. 
next_idx = current_index;

%all type 2: output is a next_location, which later is converted to new limits
        if shiftPressed == 0 && ctrlPressed == 0 && strcmp(event.Key,'a')==1 %string comparison
                %shift figure axes one event left
                if any(evtTimes == current_location)
                    if current_index-1 == 0
                        next_idx = 1;
                    else
                        next_idx = current_index -1;
                    end
                    next_location = evtTimes(next_idx);
                else
                    next_location = nearval3(current_location,evtTimes,-1);
                end
                evtTitle = 1;
            
        elseif shiftPressed == 0 && ctrlPressed == 0 && strcmp(event.Key,'d')==1 %string comparison
            %shift figure axes one event right
            if any(evtTimes == current_location)
                if current_index == length(evtTimes)
                    next_idx = current_index;
                else
                    next_idx = current_index + 1;
                end
                next_location = evtTimes(next_idx);
            else
                next_location = nearval3(current_location,evtTimes,1);
            end
            evtTitle = 1;

        elseif shiftPressed == 1 && strcmp(event.Key,'a') == 1
                %shift figure axes ten events left 
                if current_index == 0 || current_index - 10 <= 0
                    next_idx = 1;
                else 
                    next_idx = current_index - 10;
                end
                next_location = evtTimes(next_idx);
                evtTitle = 1;

        elseif shiftPressed == 1 && strcmp(event.Key,'d') == 1 
            %shift figure axes ten events right
                if current_index == length(evtTimes) || current_index + 10 >= length(evtTimes)
                    next_idx = length(evtTimes);
                else 
                    next_idx = current_index + 10;
                end
                next_location = evtTimes(next_idx);
                evtTitle = 1;
                
        elseif ctrlPressed == 1 && strcmp(event.Key,'a') == 1
            %shift figure axes 50 events left
                if current_index == 0 || current_index - 50 <= 0
                    next_idx = 1;
                else 
                    next_idx = current_index - 50;
                end
                next_location = evtTimes(next_idx);
                evtTitle = 1;
            
        elseif ctrlPressed == 1 && strcmp(event.Key,'d') == 1 
            %shift figure axes 50 events right
                if current_index == length(evtTimes) || current_index + 50 >= length(evtTimes)
                    next_idx = length(evtTimes);
                else 
                    next_idx = current_index + 50;
                end
                next_location = evtTimes(next_idx);
                evtTitle = 1;
                
        elseif strcmp(event.Key,'leftarrow')==0 && strcmp(event.Key,'f') == 1
            %shift figure axes to first event
                next_idx = 1; 
                next_location = evtTimes(next_idx);
                evtTitle = 1;
            
        elseif strcmp(event.Key,'leftarrow')==0 && strcmp(event.Key,'l') == 1 
            %shift figure axes to last event
                next_idx = length(evtTimes); 
                next_location = evtTimes(next_idx);
                evtTitle = 1;
                
        elseif strcmp(event.Key,'c') == 1
            % shift figure axes to the center event
                next_idx = ceil(length(evtTimes)/2);
                next_location = evtTimes(next_idx);
                evtTitle = 1;
                
        else
            permission = 0; % can't enter the "permission" if statement if the wrong key is pressed
            evtTitle = 1; % ** if you press an unassigned key while viewing an event, the title will 
                          % disappear unless evtTitle is something other than 0
        end
  
%% LIMITS ASSIGNMENT
    % "permission" if statement: allows only type 2 to enter -- converts
    % next_location into a new set of limits, and also prints the event
    % title if we're using the "find event" keys
    if permission == 1
        leftLim = next_location - flank;
        rightLim = next_location + flank;
        limits = [leftLim rightLim];
        if false
%         if ~isempty(usrfield)
%             str_title = sprintf('event %d/%d',next_idx,length(evtTimes));
% 
%             for iU = 1:length(usrfield)  
%                 str_font = '\fontsize{8}';
%                 str_usr{iU} = [str_font,sprintf('%s: %.3f',usrfield(iU).label,usrfield(iU).data(next_idx))];
%             end
%             
%             str = [{str_title},str_usr];
%             title(str)
           
        % no usr field input for events
        elseif evtTitle == 1
            title(sprintf('event %d/%d',next_idx,length(evtTimes)));
            
        end
    end
    
    if permission == 0 && evtTitle == 0
        title(' ')
    end
            
    % this is where the new limits are assigned to the viewing window
       set(gca,'XLim',limits)