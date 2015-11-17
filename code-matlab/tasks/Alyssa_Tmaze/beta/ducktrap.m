function ducktrap(cfg_in,varargin)
% DUCKTRAP Manually identify intervals in LFP and/or spike train data.
%
% DUCKTRAP(cfg,varargin) plots S (spiketrain) and/or CSC (LFP)and allows 
% the user to manually identify intervals by clicking on the figure to produce a 
% containment window. The order that S and CSC are passed in does not
% matter. * Do not use more than one S or more than one CSC.
%
% DUCKTRAP(cfg,S) 
% DUCKTRAP(cfg,CSC) 
% DUCKTRAP(cfg,S,CSC) or DUCKTRAP(cfg,CSC,S) 
%
% Use keyboard input to navigate to different areas of the plot (type "doc
% navigate" in the command window, or the letter h on your keyboard). 
% When you see a data region you want to keep, click on the center of the
% region (if mode is 'fixed') or click on the beginning and end (if the
% mode is 'unfixed'). The interval will plot in an orange color. If you
% want to keep the interval, hit 'Enter'. The interval will now appear in
% green. Once you've hit Enter there's no way to undo the choice, so don't 
% mess up. When an interval appears in orange, hitting any button other
% than Enter will make it disappear.
% 
% Use the Spectrogram button to plot a spectrogram under the current
% viewing window. Change the frequency range and z axis using the
% uicontrols below the Spectrogram button. Note that ducktrap will refuse 
% to do spectrograms for window sizes that are large (because ducktrap's 
% author had a bad experience one time), and will warn you for intermediate 
% window sizes that the spectrogram burden is unfavorable.
%
%   INPUTS:
%       cfg: config struct with fields controlling function behavior
%         S: spiketrain struct. Output from LoadSpikes.
%       CSC: timestamped data (LFP, ouput from LoadCSC)
%
%   OUTPUT
%       Output is saved directly to a directory specified by the user.
%       The variable name is "evt"
%       evt: iv struct with fields:
%           .tstart  - [num x 1 double] start times
%           .tend    - [num x 1 double] end times
%           .label   - the CSC used (if applicable) 
%
%   CFG OPTIONS:
%
%       cfg.mode - default 'fixed'; How the manually identified intervals 
%           are defined by the user. 
%           'fixed':  all events have the same duration set by cfg.trapwin, 
%                     and the user clicks the center time of the event to
%                     define the interval
%         'unfixed':  events will not have the same duration, instead the 
%                     user must define start and end times for the events.
%
%       cfg.trapwin - default 0.06 seconds (60 milliseconds). Applies to
%                     'fixed' mode only.
%       cfg.resume  - default []; iv struct containing previously identified
%                     intervals . (If you don't want to do it all in one
%                     sitting, save the progress and continue later by
%                     inputting the intervals here). Note: does not verify
%                     previous configs or whether the resume intervals are
%                     from the correct session.
%       cfg.sidekick - []; send inputs to sidekick() in a [1 x nArg] array.
%                      ex: [cfg_sidekick,tsd1,tsd2,iv1,tsd3]. See sidekick
%                      for more information about sidekick's plotting.
%
%       cfg.clickColor - plot color for undecided intervals; [r g b] values 
%                        or short string like 'r' or 'b' (default orange-ish)
%       cfg.keepColor  - plot color for kept intervals; [r g b] values or
%                        short string (default green)
%       cfg.EnableRobot - Java robot automatically returns focus to figure
%                         window after user interaction with UI controls. 
%                         If there is something strange in the neighborhood, 
%                         disable robot by setting cfg.EnableRobot = 0.
%
% For plot appearance options, see the config specification for MultiRaster
%
% (proposed mundane name for ducktrap: ManualIV)
% aacarey, Oct 2015 (complete rewrite from original ducktrap, Jan 2015)

%% DUCKTRAP: THE STORY
%
% ~~~ Chapter 1 ~~~
%
%   There once was a lowly biologist who liked to stalk ducks. She recorded 
% the sounds of the ducks, which were neighbours to the geese and the noisy
% waterfall. Later, she wanted to count how many times the ducks quacked, 
% but more often than not her detection would only find the loudest ducks, 
% and sometimes it would also detect the pesky geese (that incidentally also 
% liked to shit everywhere with no regard for public decency). So, the lowly 
% biologist complained to the arrogant mathematician in passing. Actually
% it was more like he overheard her unintelligibly cursing at her computer
% screen, which had done nothing wrong and was otherwise minding its
% business properly.
%
% Mildy interested in the reason for her outburst, the arrogant mathematician 
% inquired, "What the hell was that all about?"
%
% "[Censored]" she explained.
% 
% And he told her, “You biologists are so dumb, you’re probably doing 
% something like filtering the signal and looking for high power in a single 
% ripple band. If you’re just thresholding over a single frequency band you’ll 
% miss all the quiet ducks and detect all the loud geese as false positives.”
% 
% And her response was something along the lines of, “I know: I’m a lowly 
% biologist. I’m worse than lowly. In fact, I’m a plebian—no, a platyhelminth! 
% An intestinal parasite!”
% 
% Somehow this activated the arrogant mathematician’s math superiority complex, 
% “Look, if you want to find the quiet ducks, too, here’s how you should do it. 
% First, find me at least 100 real quacks.”
% 
% Supremely grateful (but also a smidgen annoyed by the arrogant mathematician’s 
% arrogance), the lowly [intestinal] biologist set out to snare some ducks.
%
% ~~~ Chapter 2 ~~~
%
%   There twice was a lowly biologist who liked to stalk ducks. Having built
% a paltry ducktrap in the first place, she had to remake it. This, of
% course, led to much frustration because of the unfortunate constraint
% that she could not built the trap out of actual things. Noticing that LB
% was redder than a MATLAB error, AM asked her what's wrong.
%
% "Well my original ducktrap is so poultry I have to remake it...(-_-)...and 
% it keeps breaking and I've appealed to the internet gods but they haven't
% answered."
%
% After taking a quick look at it, the arrogant mathematician began his 
% arrogant lecture, spewing nonsense like "state machine" and "automata" 
% and "listeners" etc... By the time he was finished, the partially-conscious 
% LB knew what a state machine was, but she mostly just wanted a donut. 

%% Set cfg parameters and check inputs

% ducktrap-specific cfg options
cfg_def.mode = 'fixed'; % 'fixed' or 'unfixed'
cfg_def.trapwin = 0.06; % window size in x axis units
cfg_def.clickColor = [255/255 99/255 71/255]; % color of intervals when initially plotted
cfg_def.keepColor = [113/255 198/255 113/255]; % color of intervals that have been kept
cfg_def.resume = []; % iv struct containing previously identified intervals
cfg_def.sidekick = [];
cfg_def.EnableRobot = 1;
%cfg_def.nansub = 4; % Every nth tsd sample is replaced by nan to speed navigation. But looks bad.

% MR-specific cfg options
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'spandex';
cfg_def.spkColor = 'k';
cfg_def.lfpColor = 'k';
cfg_def.lfpHeight = 15;
cfg_def.lfpMax = 15;
cfg_def.axislabel = 'on';
cfg_def.windowSize = 1;
cfg = ProcessConfig2(cfg_def,cfg_in); % use ProcessConfig2 because there's complications with the MR fields on cfg_in

S = []; % spike trains
cfg.lfp  = []; % continuously sampled channel 
for iVarg = 1:length(varargin)   
    if isfield(varargin{iVarg},'data') && isfield(varargin{iVarg},'label') % then it's a CSC
        cfg.lfp = varargin{iVarg}; % pass it into MultiRaster
        %cfg.lfp.data(1:cfg.nansub:end) = nan;
        
    elseif isfield(varargin{iVarg},'t')
        S = varargin{iVarg};
    end
end

if isempty(S) && ~isempty(cfg.lfp) % then user wants to see the LFP only
    S = ts;
    S.t{1}(1,1) = cfg.lfp.tvec(1); % make a fake S because MultiRaster requires this as an input in order to work
    S.t{1}(2,1) = cfg.lfp.tvec(end);
    cfg.spkColor = 'w'; % make the fake S invisible, assuming figure background is white ^_^
elseif ~isempty(S) && isempty(cfg.lfp)
    cfg = rmfield(cfg,'lfp');
elseif isempty(S) && isempty(cfg.lfp)
    error('Require CSC and/or S as inputs')
end

%% initialize some things, set global variables (any variables that are 
% here and inside of the nested functions are automatically global)

MultiRaster(cfg,S); hold on;
ax_main = gca; 
if ~isempty(cfg.sidekick)
   sidekick(cfg.sidekick{:}); hold on;
end
set(ax_main,'layer','top')
hfig = gcf; 
set(hfig,'Name',mfilename,'KeyPressFcn',@keystuff,'WindowButtonDownFcn',@clickstuff,'CloseRequestFcn',@leaveme); 

% I'm pulling some things out of cfg so I can "trace" them easier if I want to
mode = cfg.mode;
trapwin = cfg.trapwin; 
clickColor = cfg.clickColor;
keepColor = cfg.keepColor;

numKept = 0; % number of user-defined intervals
numPrev = 0;% number of events from a previous bout of identification
clicktimes = []; % locations of x-axis clicks

state = 'start'; % 'start','once','twice': state changes while function is running

x1 = []; x2 = []; y1 = []; y2 = []; H = []; % click spots and plot handle

xLoc = 1.1; yLoc = 1.08; % the horizontal and vertical location of the count text
hTxt = ''; % count text handle is global

% resume from previous session, if wanted
if ~isempty(cfg.resume)
    plot([cfg.resume.tstart cfg.resume.tend],[0 0],'LineWidth',3,'Color',[30/255 144/255 255/255])
    plot([cfg.resume.tstart cfg.resume.tend],[0 0],'o','MarkerFaceColor',[30/255 144/255 255/255],'Color',[30/255 144/255 255/255],'MarkerSize',8)
    % go to the last interval (assume they went start -> end)
    lastTime = mean([cfg.resume.tstart(end) cfg.resume.tend(end)]);
    set(gca,'XLim',[lastTime-0.5 lastTime+0.5])   
    numPrev = length(cfg.resume.tstart); % count the number of events from before
end

updateCount % displays how many intervals have been identified

%% ui control buttons

uipressed = 0; % this keeps track of when the focus is removed from the 
% figure axes; if you add new callbacks that use uicontrols, make sure you 
% keep track of uipressed and call RoboDuck at the end of the function to
% return focus to the axes so the user doesn't have to do this with their
% own clicks

% save button
uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Units','normalized','Position', [0.93 0.25 0.05 0.05],...
    'FontUnits','normalized','Callback', @saveme);

% quit button
uicontrol('Style', 'pushbutton', 'String', 'Quit',...
    'Units','normalized','Position', [0.93 0.125 0.05 0.05],...
    'FontUnits','normalized','Callback', @leaveme);

if isfield(cfg,'lfp')
    % spectrogram button
    spButton = uicontrol('Style', 'pushbutton', 'String', 'Spectrogram',...
        'Units','normalized','Position', [0.01 0.88 0.1 0.05],...
        'FontUnits','normalized','Visible','off','Callback', @spectraxis);
    
    % spectrogram z scale drop down option
    zPop = uicontrol('Style','popupmenu','String',{'root','decibel-watt','raw'},...
        'Units','normalized','Position',[0.01 0.76 0.1 0.05],'Visible','off');
    
    % spectrogram frequency range 
    frange(1) = uicontrol('Style','edit','String','0',...
        'Units','normalized','Position',[0.01 0.82 0.04 0.04],'Visible','off');
    frange(2) = uicontrol('Style','edit','String','300',...
        'Units','normalized','Position',[0.07 0.82 0.04 0.04],'Visible','off');
    
    % create exes for a spectrogram. These do not move with navigate, unlike the main
    % axes, and instead are set to invisible unless a spectrogram is plotted
    ax_spec = axes('Position', get(gca, 'Position'),'Visible','off');
    set(hfig,'CurrentAxes',ax_main)
end

%% helper functions

    function updateCount
        % how many intervals have been indentified. Displays text.
        txt = ['Count: ',num2str(numKept + numPrev)];
        hTxt = text(xLoc,yLoc,txt,'Units','normalized',...
            'VerticalAlignment','top',...
            'HorizontalAlignment','right',...
            'FontSize',20,'BackgroundColor','w');       
    end % of updateCount

    function RoboDuck
        % use java robot to help return focus to the axes after interacting with a uicontrol
        % RoboDuck and uipressed work together to produce the desired behaviour
        
        if cfg.EnableRobot
            % get original location of mouse cursor
            ml_orig = get(0,'PointerLocation');
            
            % get center of figure window
            figLoc = get(hfig,'Position');
            ml_new = [figLoc(1)+figLoc(3)/2 figLoc(2)+figLoc(4)/2];
            
            % set new mouse location for the robot
            set(hfig,'Pointer','custom','PointerShapeCData',NaN(16,16)) % make it invisible
            set(0,'PointerLocation',ml_new) 
            
            % robot!
            drawnow
            robot = java.awt.Robot ;
            robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK);
            robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);
            
            % now return pointer location to where the user last had it
            set(hfig,'Pointer','arrow') % make it visible again
            set(0,'PointerLocation',ml_orig)
            
        end
    end

%% callback functions

% ~~~~~~ CLICK STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function clickstuff(~,~)
        % handles mouse input
        
        if ~uipressed
            clickpoint = get(gca,'CurrentPoint');
            x = clickpoint(1,1); y = clickpoint(1,2);
            
            switch state
                case 'start'
                    switch mode
                        case 'fixed' % interval is defined by trapwin around click location
                            state = 'twice';
                            x1 = x; y1 = y;
                            x_plot = [x1-trapwin/2 x1+trapwin/2];
                            H(1) = plot(gca,x_plot,[y1 y1],'s','MarkerFaceColor',clickColor,'Color',clickColor,'MarkerSize',8);
                            H(2) = plot(gca,x_plot,[y1 y1],'LineWidth',3,'Color',clickColor);
                            
                        case 'unfixed' % interval start and stop are defined by user
                            state = 'once';
                            x1 = x;
                            y1 = y;
                            H(1) = plot(gca,x1,y1,'s','MarkerFaceColor',clickColor,'Color',clickColor,'MarkerSize',8);
                            
                        otherwise
                            error('Unrecognized mode')
                    end
                    
                case 'once'
                    assert(strcmp(mode,'unfixed'))
                    state = 'twice';
                    x2 = x; y2 = y;
                    H(2) = plot(gca,x2,y2,'s','MarkerFaceColor',clickColor,'Color',clickColor,'MarkerSize',8);
                    H(3) = plot(gca,[x1 x2],[y1 y2],'LineWidth',3,'Color',clickColor);
                    
                case 'twice'
                    state = 'start';
                    delete(H);
            end
        end
        uipressed = 0;
    end % of clickstuff

% ~~~~~~ KEY STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function keystuff(source,event)
        % handles keyboard input
        if strcmp(event.Key,'return')
            switch state
                case 'start'
                    %Do nothing...
                    
                case 'once'
                    state = 'start';
                    delete(H)
                    
                case 'twice'
                    state = 'start';
                    switch mode
                        case 'fixed' % interval is defined by trapwin around click location
                            set(H(1),'Color',keepColor,'Marker','o','MarkerFaceColor',keepColor); set(H(2),'Color',keepColor)
                            numKept = numKept + 1;
                            clicktimes(numKept) = x1;
                            
                            if ishandle(hTxt)
                                delete(hTxt)
                            end
                            
                            % update count display
                            updateCount
                            
                        case 'unfixed' % interval start and stop are defined by user
                            
                            set(H,'Marker','o','Color',keepColor,'MarkerFaceColor',keepColor);
                            numKept = numKept + 1;
                            
                            % make sure tstart is before tend
                            if x1 < x2
                                clicktimes.tstart(numKept) = x1;
                                clicktimes.tend(numKept) = x2;
                            else
                                clicktimes.tstart(numKept) = x2;
                                clicktimes.tend(numKept) = x1;
                            end
                            
                            if ishandle(hTxt)
                                delete(hTxt)
                            end
                            
                            % update count display
                            updateCount
                            
                        otherwise
                            error('Unrecognized mode')
                    end
            end
            H = [];
        else 
            if ~isempty(cfg.lfp)
                set(spButton,'Visible','on'); set(zPop,'Visible','on'); set(frange,'Visible','on')
                mafak = findobj(ax_spec); % use findobj, otherwise can't make imagesc invisible (just axes)
                set(mafak,'Visible','off'); set(ax_main,'Color','w') 
            end
            navigate(source,event)
        end
    end % of keystuff

% ~~~~~~ LEAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function leaveme(~,~)
        % quit request function
        uipressed = 1;
        choice = questdlg('Are you sure you want to close the figure? Any unsaved data will be lost.','Quit Requested','Yes','No','No');
        switch choice
            case 'Yes'
                delete(gcf)
            case 'No'
                RoboDuck
                return
        end
    end % of leaveme

% ~~~~~~ SAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function saveme(~,~)
        % save request function
        uipressed = 1;
        choice = questdlg('Would you like to save the data?','Save Requested','Yes','No','No');
        switch choice
            case 'Yes'
                switch mode
                    case 'fixed'
                        if ~isempty(cfg.resume)
                           clicktimes = [clicktimes IVcenters(cfg.resume)]; 
                        end
                        clicktimes = sort(clicktimes);
                        evt = iv(clicktimes - trapwin/2,clicktimes + trapwin/2);
                        if isfield(cfg,'lfp')
                            evt.label = cfg.lfp.label;
                        end
                        
                    case 'unfixed'
                        if ~isempty(cfg.resume)
                           clicktimes.tstart = [clicktimes.tstart cfg.resume.tstart']; 
                           clicktimes.tend = [clicktimes.tend cfg.resume.tend']; 
                        end
                        [clicktimes.tstart,idx] = sort(clicktimes.tstart);
                        clicktimes.tend = clicktimes.tend(idx);
                        evt = iv(clicktimes.tstart,clicktimes.tend);  
                        if isfield(cfg,'lfp')
                            evt.label = cfg.lfp.label;
                        end
                end
                
                % check that there aren't any doubles
                discard = find(diff(evt.tstart)== 0);
                discard2 = find(diff(evt.tend) == 0);
                if any(discard2 ~= discard)
                    warning('Some intervals have the same start or end times.')
                elseif ~isempty(discard)
                    disp([mfilename,': doubles found, removing.'])
                    evt.tstart(discard) = [];
                    evt.tend(discard) = [];
                end
                        
                [~,name,~] = fileparts(pwd);
                uisave('evt',[name,'-manualIV']) % opens window for saving stuff
                
                RoboDuck
                
            case 'No'
                RoboDuck
                return
        end
    end % of saveme

% ~~~~~~ SPECTRAXIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function spectraxis(~,~)
        % create a new set of axes and plot a spectrogram behind the raster plot
        uipressed = 1;
        
        xlims = get(gca,'XLim'); 
        xticks = get(gca,'XTick');
        
        ugh = restrict2(cfg.lfp,xlims(1),xlims(2));
        
        if length(ugh.tvec) > 50000
            % outright refuse
            msgbox('I refuse to spectrogram this much data...it''s for your own good.','Bad things can happen','error')
            RoboDuck
            return
        end
        
        if length(ugh.tvec) > 20000
            choice = questdlg('Oh wow! Your computer figuratively almost caught fire because you wanted to spectrogram a lot of data.','This might take a long time...','I like fire','Oops!','Oops!');
          
            switch choice
                case 'I like fire'
                case 'Oops!'
                    RoboDuck
                    return  
            end
        end
        
        set(gca,'color','none')
        set(hfig,'CurrentAxes',ax_spec) % now change current axes
        
        fs = 2000; % sampling frequency
        
        % get frequencies of interest
        foi = str2double(get(frange(1),'String')):str2double(get(frange(2),'String')); % don't be evil..numbers only
        
        [~,F,T,P] = spectrogram(ugh.data,hanning(102),100,foi,fs);
        
        % get z scale option (see uicontrol zPop)
        string = get(zPop,'String'); choice = get(zPop,'Value');
        zscale = string{choice};
        switch zscale
            case 'root'
                P = sqrt(P);
            case 'decibel-watt'
                P = 10*log10(P);
            case 'raw'
                % do nothing
        end
        
        imagesc(T,F,P)
        set(ax_spec,'YAxisLocation','right','XTick',xticks,'XTickLabel',[],'YDir','normal','Visible','on')
        
        ylabel(ax_spec,'Frequency (Hz)')
        grid on
        uistack(ax_spec,'bottom') % send it to the very back; other axes in front
        % return to main axes (for navigation and stuff)
        set(hfig,'CurrentAxes',ax_main) 
        
        RoboDuck
        
    end % of spectraxis

end

