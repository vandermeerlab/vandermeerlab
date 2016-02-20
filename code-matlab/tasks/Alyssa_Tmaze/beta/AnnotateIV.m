function AnnotateIV(cfg_in,iv_in,varargin)
%ANNOTATEIV Add annotations to iv data
% ANNOTATEIV(cfg,iv,varargin) plots S (spiketrains) and/or csc (LFP)and
% provides uicontrols allowing the user to move through intervals while
% annotating or deleting them. The order that S and csc are passed in does not
% matter. *Do not use more than one S or more than one csc.
%
% ANNOTATEIV(cfg,iv,S) 
% ANNOTATEIV(cfg,iv,CSC) 
% ANNOTATEIV(cfg,iv,S,CSC) or ANNOTATEIV(cfg,iv,CSC,S)
%
% Use keyboard input to navigate to different areas of the plot (type "doc
% navigate" in the command window, or the letter h on your keyboard).
%
% Special keyboard commands for AnnotateIV (unassigned to navigate at the
% time this was written):
%
%  n - go to the next event
%  p - go to the previous event
%  s - spectrogram
%
%  Additional keyboard commands are available if cfg.EnableRobot = 1
%  t - (type) go to the annotate box (when in the annotation box, hit 'Enter' to
%      return focus to the figure window and continue using keyboard commands)
%  r - delete (remove) this event upon saving
%
%  Tip: if the keyboard commands are not working, try clicking on the
%  figure window to return focus.
%
%   OUTPUT
%    - is variable called "evt" that is saved to directory of your choice.
%      evt is iv data type with fields:
%         .cfg
%         .tstart
%         .tend
%         .label (CSC used, if applicable)
%         .usr.annotation: [nIntervals x 1] cell array containing 
%                 the annotations you have given each interval. 
%          Note: the usr field was created to hold any [nIntervals x 1] 
%                data pertaining to each interval. The dimensions of these
%                fields must be the same as .tstart and .tend.
%
%   CONFIG OPTIONS
%       cfg.mode = 'random'; In which order should the intervals be
%                                shown?
%                'random'  -  Show me the intervals in random order.
%                'chrono'  -  Show me the intervals in chronological order.
%
%       cfg.EnableRobot - Java robot automatically returns focus to figure
%                         window after user interaction with UI controls.
%                         It also makes complete keyboard control of
%                         AnnotateIV possible by performing mouse movements
%                         when the user clicks 't' and 'r'.
%                         If there is something strange try disabling the
%                         robot by setting cfg.EnableRobot = 0 to see if
%                         this fixes the problem.
%
%
% This function was written for testing detector performance against
% human-identified SWRs. It allows a user to rate the quality of the SWR so
% they can compare the intervals found by detectors to the intervals found
% by humans. Low ratings reflect low confidence and high ratings reflect high
% confidence. 
% Consider the following qualities when rating intervals:
% - frequency content (use the spectrogram button)
% - associated spiking activity 
% - number of ripples in SWRs / duration
% - amplitude of ripples in SWRs
% - the "bundle-iness" of the SWR 
%
% All text entered into the box must begin with a number (or
% otherwise, some sequence of characters establised by the user).
%
% Intervals flagged for deletion are removed upon saving; i.e. you can
% change your mind about deleting an interval at any time as long as the
% figure window remains open. 
%
% aacarey Oct 2015, edit Dec 2015

cfg_def.mode = 'random'; % 'chrono','random'
cfg_def.EnableRobot = 1;

% MR-specific cfg options
cfg_in.evt = iv_in; % don't let them put in the wrong evt via config
cfg_def.evt = iv_in;
cfg_def.ivColor = 'r';
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'spandex';
cfg_def.spkColor = 'k';
cfg_def.lfpColor = 'k';
cfg_def.lfpHeight = 15;
cfg_def.lfpMax = 15;
cfg_def.axislabel = 'on';
cfg_def.windowSize = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);

is_iv = CheckIV(iv_in);
if ~is_iv
    error('iv_in must be an iv data type.')
end

% check for iv doubles
checkstart = find(diff(iv_in.tstart)== 0);
checkend = find(diff(iv_in.tend) == 0);
if any(checkstart ~= checkend)
    warning('Some intervals have the same start or end times.')
elseif ~isempty(checkstart)
    error('iv data contains doubles.')
end

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
    S.t{1}(1,1) = cfg.lfp.tvec(1); % make a fake S because MultiRaster requires this as an input in order to work
    S.t{1}(2,1) = cfg.lfp.tvec(end);
    cfg.spkColor = 'w'; % make the fake S invisible, assuming figure background is white ^_^
elseif ~isempty(S) && isempty(cfg.lfp)
    cfg = rmfield(cfg,'lfp');
elseif isempty(S) && isempty(cfg.lfp)
    error('Require csc and/or S as inputs')
end

% keep track of intervals
switch cfg.mode
    case 'random' % user wants to be shown intervals in a random order
        s = rng; % get current rng state
        rng(1,'twister') % seed the random number generator
        idx_list_all = randperm(length(iv_in.tstart));
        rng(s); % return to original state
        
    case 'chrono' % user wants to be shown intervals in chronological order
        idx_list_all = 1:length(iv_in.tstart);
end

% indices into iv data in the order that we will see them. when "Show unseen" is
% selected, idx_list is updated to include only unseen intervals. when 'show
% all' is selected, idx_list is reset to idx_list_all
idx_list = idx_list_all;

ii = 1;
current_iv = idx_list(ii);
midpoints = IVcenters(iv_in);

if ~isfield(iv_in,'usr') || ~isstruct(iv_in.usr) || ~isfield(iv_in.usr,'annotation')
    % add empty annotations
    iv_in.usr.annotation = cell(length(iv_in.tstart),1);
end

del = zeros(1,length(iv_in.tstart)); % del will become 1 where the user wants to delete an interval.

% plot the stuff
hMR = MultiRaster(cfg,S); hold on;
hfig = gcf; 
set(hfig,'Name',mfilename,'KeyPressFcn',@keystuff,'WindowButtonDownFcn',@clickstuff,'CloseRequestFcn',@leaveme);

% handle focus IV plotting
h_curr = [];
windowSize = 1;

ax_main = gca;

uipressed = 0; % this helps track whether a UI button has been pressed when handling RoboDuck and human clicks in clickstuff function
inRedraw = 0; % whether we are currently redrawing an interval

% this keeps track of clicks for redrawing intervals (see clickstuff)
state = 'start';
% 'start' is before clicks
% 'once' is a single click
% 'twice' is two clicks
x1 = []; x2 = []; % click spots for redrawing intervals

% if the interval is from ducktrap, we need to know which mode was used
if isfield(iv_in,cfg) && isfield(iv_in.cfg,'history') && any(strcmp(iv_in.cfg.history.mfun,'ducktrap'))
    cfg_temp = [];
    cfg_temp.target = 'ducktrap';
    
    history = GetHistory(cfg_temp,iv_in);
    mode = history{end}.mode;
    
    % if mode was fixed, we need trapwin so all intervals have same
    % duration
    if strcmp(mode,'fixed')
        trapwin = history{end}.trapwin;
    end
    
else % if it's not from ducktrap, just make user click twice because w/e
    mode = 'unfixed';
end

updateProgress

%% ui control buttons

% save button
uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Units','normalized','Position', [0.925 0.175 0.05 0.05],...
    'FontUnits','normalized','Callback', @saveme);

% quit button
uicontrol('Style', 'pushbutton', 'String', 'Quit',...
    'TooltipString',['Exit ',mfilename],...
    'Units','normalized','Position', [0.925 0.125 0.05 0.05],...
    'FontUnits','normalized','Callback', @leaveme);

% next button
hnext = uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'TooltipString','Go to next interval',...
    'Units','normalized','Position', [0.92 0.75 0.075 0.05],...
    'FontUnits','normalized','Callback', @nextIV);

% prev button
hprev = uicontrol('Style', 'pushbutton', 'String', 'Previous',...
    'TooltipString','Go to previous interval',...
    'Units','normalized','Position', [0.92 0.88 0.075 0.05],...
    'FontUnits','normalized','Enable','off','Callback', @prevIV);

% collect annotation
hnote = uicontrol('Style', 'edit', 'String', '',...
    'TooltipString','Type annotation here',...
    'Units','normalized','Position', [0.92 0.83 0.075 0.05],...
    'FontUnits','normalized','Callback', @writestuff,'KeyPressFcn',@NotLetters);

% delete checkbox
hdel = uicontrol('Style', 'checkbox', 'String', 'Delete this interval',...
    'Units','normalized','Position', [0.92 0.8 0.075 0.03],...
    'Callback', @dontwant);
     set(hdel,'BackgroundColor',get(gcf,'Color')); % checkbox otherwise sits in a white rectangle, which looks bad
     
% show certain intervals drop down option
hChange = uicontrol('Style','popupmenu','String',{'Show all','Show unseen'},...
        'TooltipString','Which intervals are shown',...
        'Units','normalized','Position',[0.92 0.68 0.075 0.05],'Callback',@ChangeList);     
     
% redraw button
hRedraw = uicontrol('Style', 'pushbutton', 'String', 'Redraw',...
    'TooltipString','Redefine current interval boundaries',...
    'Units','normalized','Position', [0.92 0.64 0.075 0.05],...
    'FontUnits','normalized','Callback', @EnableRedraw);

% cancel redraw button
hCancelRedraw = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
    'TooltipString','Return to annotating without completing redraw',...
    'Units','normalized','Position', [0.92 0.59 0.075 0.05],...
    'FontUnits','normalized','Enable','off','Callback', @DisableRedraw);

if isfield(cfg,'lfp')
    % spectrogram button
    uicontrol('Style', 'pushbutton', 'String', 'Spectrogram',...
        'TooltipString','Plot spectrogram for current window',...
        'Units','normalized','Position', [0.01 0.88 0.1 0.05],...
        'FontUnits','normalized','Callback', @spectraxis);
    
    % spectrogram z scale drop down option
    zPop = uicontrol('Style','popupmenu','String',{'root','decibel-watt','raw'},...
        'TooltipString','Choose colour axis scaling',...
        'Units','normalized','Position',[0.01 0.76 0.1 0.05]);
    
    % spectrogram frequency range 
    frange(1) = uicontrol('Style','edit','String','50',...
        'TooltipString','Pass band lower frequency',...
        'Units','normalized','Position',[0.01 0.82 0.04 0.04]);
    frange(2) = uicontrol('Style','edit','String','300',...
        'TooltipString','Pass band higher frequency',...
        'Units','normalized','Position',[0.07 0.82 0.04 0.04]);
    
    % create exes for a spectrogram. These do not move with navigate, unlike the main
    % axes, and instead are set to invisible unless a spectrogram is plotted
    ax_spec = axes('Position', get(gca, 'Position'),'Visible','off');
    set(hfig,'CurrentAxes',ax_main)
end

%% this needs to be called AFTER the uicontrols because it makes reference to some of them

PlotCurrentIV

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    HELPER FUNCTIONS                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~~~~~~ UPDATE PROGRESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function updateProgress
        xLoc = 1.1; yLoc = 1.08; % the horizontal and vertical location of the progress text
        temp = ~cellfun(@isempty,iv_in.usr.annotation);
        temp = temp' | del;
        temp = sum(temp)/length(iv_in.usr.annotation)*100;
        
        txt = ['Progress: ',sprintf('%.2f',temp),'%'];
        text(xLoc,yLoc,txt,'Units','normalized',...
            'VerticalAlignment','top',...
            'HorizontalAlignment','right',...
            'FontSize',20,'BackgroundColor','w');
    end % of updateProgress

% ~~~~~~ UPDATE BUTTONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function updateButtons
        % update edit box text to show current annotation string
        if isempty(iv_in.usr.annotation{current_iv})
            set(hnote,'String','')
        else
            set(hnote,'String',iv_in.usr.annotation{current_iv})
        end
        
        % update delete interval choice
        set(hdel,'Value',del(current_iv))
    end

% ~~~~~~ PLOT CURRENT IV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function PlotCurrentIV 
        set(gca,'XLim',[midpoints(current_iv)-windowSize/2 midpoints(current_iv)+windowSize/2])
        delete(h_curr)
        ylims = get(gca,'Ylim');
        h_curr(1) = plot([iv_in.tstart(current_iv) iv_in.tstart(current_iv)],ylims);
        h_curr(2) = plot([iv_in.tend(current_iv) iv_in.tend(current_iv)],ylims);
        set(h_curr, 'Color',cfg.ivColor,'LineWidth',1.5)
        title(['Annotate event ',num2str(current_iv),' out of ',num2str(length(iv_in.tstart))],'FontSize',14)
        drawnow % this makes it wait for the plotting to happen before it updates the buttons (plotting very slow)
        updateButtons
    end

% ~~~~~~ ROBODUCK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function RoboDuck(obj_handle)
        % use java robot to help return focus to the axes after interacting with a uicontrol
        % RoboDuck and uipressed work together to produce the desired behaviour
        % RoboDuck also performs some clicking actions such as moving the
        % cursor to the edit and check boxes.
        % obj_handle is the handle of the object you want RoboDuck to click
        % for you
        
        if cfg.EnableRobot && ~inRedraw
            % get original location of mouse cursor
            ml_orig = get(0,'PointerLocation');
            
            % get location of figure window
            objLoc = get(obj_handle,'Position');
            
            type = get(obj_handle,'Type');
            
            switch type
                case 'figure'
                    % in this case, the handle position is with reference
                    % to the screen
                    ml_new = [objLoc(1)+objLoc(3)/2 objLoc(2)+objLoc(4)/2];
                    set(hfig,'Pointer','custom','PointerShapeCData',NaN(16,16)) % make cursor invisible
                    
                case 'uicontrol'
                    % in this case, the handle position is with reference
                    % to the parent (which is the figure)
                    
                    parentLoc = get(hfig,'Position'); % note can also do this by getting the parent from the uicontrol handle
                    ml_new = [parentLoc(1)+objLoc(1)*parentLoc(3)+(objLoc(3)*parentLoc(3))/2,parentLoc(2)+objLoc(2)*parentLoc(4)+(objLoc(4)*parentLoc(4))/2];         
            end
            
            % set new mouse location for the robot
            set(0,'PointerLocation',ml_new)
            
            % robot!
            drawnow
            robot = java.awt.Robot ;
            robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK);
            robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);
            
            if strcmp(type,'figure')
                set(hfig,'Pointer','arrow') % make cursor visible again
            end
            
            % now return pointer location to where the user last had it
            set(0,'PointerLocation',ml_orig)
            
        end
    end % of RoboDuck

% ~~~~~~ HIDE SPECTRAXIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function HideSpectraxis
        specObj = findobj(ax_spec); % use findobj, otherwise can't make imagesc invisible (just axes)
        set(specObj,'Visible','off'); set(ax_main,'Color','w')
        
        % return color of LFP/spikes/intervals to have nice contrast with white background
        %if isfield(hMR,'S'); set(hMR.S(:),'Color','k'); end
        if isfield(hMR,'LFP'); set(hMR.LFP,'Color','k'); end
        if isfield(hMR,'LFP_iv'); set(hMR.LFP_iv,'Color','r'); end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    CALLBACK FUNCTIONS                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ~~~~~~ NEXT INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function nextIV(~,~)
        
        uipressed = 1;
        
        HideSpectraxis
        if ii == length(idx_list) % we've gone through all of them
            set(hnext,'Enable','off')
            title('End of interval list','FontSize',14)
        else
            ii = ii+1;
            current_iv = idx_list(ii);
            PlotCurrentIV
        end  
        set(hprev,'Enable','on')
        
        RoboDuck(hfig)
    end

% ~~~~~~ PREVIOUS INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function prevIV(~,~)
        
        uipressed = 1;
        
        HideSpectraxis
        if ii-1 <= 0 % we're at the beginning
            set(hprev,'Enable','off')
        else
            ii = ii-1;
            current_iv = idx_list(ii);
            PlotCurrentIV
        end
        set(hnext,'Enable','on')
        
        RoboDuck(hfig)
    end

% ~~~~~~ CHANGE LIST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ChangeList(~,~)
        % change idx_list between the whole set of intervals and the
        % non-annotated subset
        uipressed = 1;
        string = get(hChange,'String'); selected = get(hChange,'Value');
        choice = string{selected};
        switch choice
            case 'Show unseen'
                % remove seen intervals from idx_list
                temp = iv_in.usr.annotation;
                temp(logical(del)) = [];
                unseen = cellfun(@isempty,temp);
                unseen = unseen(idx_list);
                idx_list = idx_list(unseen);
                
                ii = 1; current_iv = idx_list(ii);
                PlotCurrentIV
                
            case 'Show all'
                idx_list = idx_list_all;
                ii = find(idx_list == current_iv);
        end
        RoboDuck(hfig)
    end % of changelist

% ~~~~~~ KEY STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function keystuff(source,event)
        if strcmp(event.Key,'n')
            nextIV(source,event)
        elseif strcmp(event.Key,'p')
            prevIV(source,event)
        elseif strcmp(event.Key,'t')
            RoboDuck(hnote)
        elseif strcmp(event.Key,'r')
            RoboDuck(hdel)
        elseif strcmp(event.Key,'s')
            spectraxis(source,event)
        else
            HideSpectraxis
            navigate(source,event)
        end
    end

% ~~~~~~ NOT LETTERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function NotLetters(~, event)
        % keypress function for keys other than letter/character keys,
        % intended just to detect when Enter is pressed for the annotation
        % and delete box so that focus can be returned to the figure window
        if strcmp(event.Key,'return')
            RoboDuck(hfig)
        end
    end

% ~~~~~~ WRITE STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function writestuff(source,~)
        % add annotation to iv
        
        uipressed = 1;
        
        str = strtrim(get(source, 'String')); 
        % remove leading and trailing whitespace
        % before str = '    '; now str = ''; or before str = ' 5 '; now str = '5';
        iv_in.usr.annotation{current_iv} = str;
        
        updateProgress
    end

% ~~~~~~ CLICK STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function clickstuff(~,~)
        % handle mouse clicks
        if ~uipressed
            HideSpectraxis
        end
        uipressed = 0;
        
        % right now this is just for hiding spectrogram axis, but can add
        % more mouse things here
        
    end

% ~~~~~~ ENABLE REDRAW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function EnableRedraw(~,~)
        % changes the current window button down function to redraw
        set(ax_main,'XLim',[midpoints(current_iv)-windowSize/2 midpoints(current_iv)+windowSize/2])
        set(hfig,'WindowButtonDownFcn',@redraw)
        set(hRedraw,'Enable','off')
        set(hCancelRedraw,'Enable','on')
        switch mode
            case 'fixed'
                title('Redraw enabled: click center of interval','FontSize',14)
            case 'unfixed'
                title('Redraw enabled: click beginning of interval','FontSize',14)
        end
    end % of EnableRedraw

% ~~~~~~ DISABLE REDRAW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function DisableRedraw(~,~)
        % changes the current window button down function to clickstuff
        uipressed = 1;
        inRedraw = 0; % whether we are currently redrawing an interval
        state = 'start';
        set(hfig,'WindowButtonDownFcn',@clickstuff)
        set(hRedraw,'Enable','on') % allow user to click button for next time
        set(hCancelRedraw,'Enable','off')
        PlotCurrentIV
        RoboDuck(hfig)
    end % of DisableRedraw

% ~~~~~~ REDRAW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function redraw(source,event)
        % handles mouse input for redrawing interval (this is modified from
        % ducktrap's clickstuff
        
        inRedraw = 1; % whether we are currently redrawing an interval
        ylims = get(ax_main,'Ylim');
        
        clickpoint = get(ax_main,'CurrentPoint');
        x = clickpoint(1,1);
        
        switch state
            case 'start'
                switch mode
                    case 'fixed' % interval is defined by trapwin around click location
                        state = 'twice';
                        x1 = x;
                        % assign new iv boundaries
                        iv_in.tstart(current_iv) = x1-trapwin/2;
                        iv_in.tend(current_iv) = x1+trapwin/2;
                        
                        % update plot of current iv
                        PlotCurrentIV
                        
                        DisableRedraw(source,event)
                        
                    case 'unfixed' % interval start and stop are defined by user
                        state = 'once';
                        x1 = x;
                        delete(h_curr(1))
                        h_curr(1) = plot(ax_main,[x1 x1],ylims,'Color',cfg.ivColor,'LineWidth',1.5);
                        title('Redraw enabled: click end of interval','FontSize',14)
                        
                    otherwise
                        error('Unrecognized mode')
                end
                
            case 'once'
                assert(strcmp(mode,'unfixed'))
                state = 'twice';
                x2 = x;
                delete(h_curr(2))
                h_curr(2) = plot(ax_main,[x2 x2],ylims,'Color',cfg.ivColor,'LineWidth',1.5);
                
                % reassign iv boundaries
                clicktimes = sort([x1 x2]);
                iv_in.tstart(current_iv) = clicktimes(1);
                iv_in.tend(current_iv) = clicktimes(2);
                drawnow
                DisableRedraw(source,event)
                
            case 'twice'
                % if user has clicked twice, cycle is complete - return
                % to 'start' state
                state = 'start';
        end
    end % of redraw

% ~~~~~~ DON'T WANT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function dontwant(source,~)
        % keeps a record of which intervals to delete
        
        uipressed = 1;
        
        del(current_iv) = get(source, 'Value');
        
        updateProgress 
        RoboDuck(hfig)        
    end

% ~~~~~~ LEAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function leaveme(~,~)
        % quit request function
        
        uipressed = 1;
        
        temp = iv_in.usr.annotation;
        temp(logical(del)) = [];
        nEmpties = length(temp(cellfun(@isempty,temp)));
        
        if nEmpties > 0
            choice = questdlg(['You are missing ',num2str(nEmpties),' annotations. Are you sure you want to close the figure? Any unsaved data will be lost.'],'Quit Requested','Quit anyway','Cancel','Cancel');
            switch choice
                case 'Cancel'
                    RoboDuck(hfig)
                    return
                case 'Quit anyway'
                    delete(gcf)
            end
        else  
            choice = questdlg('Are you sure you want to close the figure? Any unsaved data will be lost.','Quit Requested','Yes','No','No');
            switch choice
                case 'Yes'
                    delete(gcf)
                case 'No'
                    RoboDuck(hfig)
                    return
            end
        end
    end % of leaveme

% ~~~~~~ SAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function saveme(~,~)
        % save request function
        
        uipressed = 1;
        
        choice = questdlg('Would you like to save the data?','Save Requested','Yes','No','No');
        switch choice
            case 'Yes'
                evt = iv_in;
                evt.tstart(logical(del)) = [];
                evt.tend(logical(del)) = [];
                evt.usr.annotation(logical(del)) = [];
                
                % check that there aren't any doubles
                discard = find(diff(evt.tstart)== 0);
                discard2 = find(diff(evt.tend) == 0);
                if any(discard2 ~= discard)
                    warning('Some intervals have the same start or end times.')
                elseif ~isempty(discard)
                    disp([mfilename,': doubles found, removing.'])
                    evt.tstart(discard) = [];
                    evt.tend(discard) = [];
                    evt.usr.annotation(discard) = [];
                end
                
                [~,name,~] = fileparts(pwd);
                uisave('evt',[name,'-manualIV']) % opens window for saving stuff
                
                RoboDuck(hfig)
                
            case 'No'
                RoboDuck(hfig)
                return
        end
    end % of saveme

% ~~~~~~ SPECTRAXIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function spectraxis(~,~)
        % create a new set of axes and plot a spectrogram behind the raster plot
        
        uipressed = 1;
        
        set(gca,'color','none')
        xlims = get(gca,'XLim'); 
        xticks = get(gca,'XTick');
        set(hfig,'CurrentAxes',ax_spec) % now change current axes
        
        nSamples = 102; % the number of samples to use in the spectrogram
       
        CSCr = restrict(cfg.lfp,xlims(1),xlims(2)); % restrict to count the numbers of samples in viewing window
        
        fs = length(CSCr.data)/(xlims(2)-xlims(1)); % get the local approx sampling frequency
        
        buffer = (nSamples/2)/fs; % how much time buffer is needed so that the spectrogram lines up with the data in viewing window
        CSCr = restrict(cfg.lfp,xlims(1)-buffer,xlims(2)+buffer); % re-restrict with buffer
        
        % get frequencies of interest
        foi = str2double(get(frange(1),'String')):str2double(get(frange(2),'String')); % don't be evil..numbers only
        
        [~,F,T,P] = spectrogram(CSCr.data,hanning(nSamples),100,foi,fs);
        
        % get z scale option (see uicontrol zPop)
        string = get(zPop,'String'); choice = get(zPop,'Value');
        zscale = string{choice};
        switch zscale
            case 'root'
                P = sqrt(P);
                col = [0.09*10^-6 18*10^-6]; % some arbitrary range for the color scaling, keep everything relative or w/e
                
                % change colors of LFP/spikes/intervals for contrast with spectrogram
                %if isfield(hMR,'S'); set(hMR.S,'Color','w'); end
                if isfield(hMR,'LFP'); set(hMR.LFP,'Color','w'); end
                if isfield(hMR,'LFP_iv'); set(hMR.LFP_iv,'Color','r'); end
                
            case 'decibel-watt'
                P = 10*log10(P);
                col = [-170 -80];
                
                % change colors of LFP/spikes/intervals for contrast with spectrogram
                %if isfield(hMR,'S'); set(hMR.S(:),'Color','k'); end
                if isfield(hMR,'LFP'); set(hMR.LFP,'Color','k'); end
                if isfield(hMR,'LFP_iv'); set(hMR.LFP_iv,'Color','b'); end
                
            case 'raw'
                col = [0.09*10^-10 3*10^-10]; 
                
                % change colors of LFP/spikes/intervals for contrast with spectrogram
                %if isfield(hMR,'S'); set(hMR.S(:),'Color','w'); end
                if isfield(hMR,'LFP'); set(hMR.LFP,'Color','w'); end
                if isfield(hMR,'LFP_iv'); set(hMR.LFP_iv,'Color','r'); end
        end
        
        imagesc(T,F,P,col); 
        set(ax_spec,'YAxisLocation','right','XTick',xticks,'XTickLabel',[],'YDir','normal','Visible','on')
        
        ylabel(ax_spec,'Frequency (Hz)')
        grid on
        uistack(ax_spec,'down') % send backwards
        % return to main axes (for navigation and stuff)
        set(hfig,'CurrentAxes',ax_main) 
        
        RoboDuck(hfig)
        
    end % of spectraxis

end

