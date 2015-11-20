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
% navigate" in the cammand window, or the letter h on your keyboard).
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
%       cfg_def.mode = 'random'; In which order should the intervals be
%                                shown?
%                'random'  -  Show me the intervals in random order.
%                'chrono'  -  Show me the intervals in chronological order.
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
% aacarey Oct 2015

cfg_def.mode = 'random'; % 'chrono','random'

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
        idx_list = randperm(length(iv_in.tstart));       
    case 'chrono' % user wants to be shown intervals in chronological order
        idx_list = 1:length(iv_in.tstart);
end
ii = 1;
current_iv = idx_list(ii);
midpoints = IVcenters(iv_in);

if ~isfield(iv_in,'usr') || ~isstruct(iv_in.usr) || ~isfield(iv_in.usr,'annotation')
    % add empty annotations
    iv_in.usr.annotation = cell(length(iv_in.tstart),1);
end

del = zeros(1,length(iv_in.tstart)); % del will become 1 where the user wants to delete an interval.

% plot the stuff
MultiRaster(cfg,S); hold on;
hfig = gcf; 
set(hfig,'Name',mfilename,'KeyPressFcn',@keystuff,'WindowButtonDownFcn',@clickstuff,'CloseRequestFcn',@leaveme);

% handle focus IV plotting
h_curr = [];
windowSize = 1;

ax_main = gca;


updateProgress

%% ui control buttons

% save button
uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Units','normalized','Position', [0.93 0.25 0.05 0.05],...
    'FontUnits','normalized','Callback', @saveme);

% quit button
uicontrol('Style', 'pushbutton', 'String', 'Quit',...
    'Units','normalized','Position', [0.93 0.125 0.05 0.05],...
    'FontUnits','normalized','Callback', @leaveme);

% next button
hnext = uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Units','normalized','Position', [0.9075 0.75 0.09 0.05],...
    'FontUnits','normalized','Callback', @nextIV);

% prev button
hprev = uicontrol('Style', 'pushbutton', 'String', 'Previous',...
    'Units','normalized','Position', [0.9075 0.88 0.09 0.05],...
    'FontUnits','normalized','Enable','off','Callback', @prevIV);

% collect annotation
hnote = uicontrol('Style', 'edit', 'String', '[enter text]',...
    'Units','normalized','Position', [0.9075 0.83 0.09 0.05],...
    'FontUnits','normalized','Callback', @writestuff);

% delete checkbox
hdel = uicontrol('Style', 'checkbox', 'String', 'Delete this interval',...
    'Units','normalized','Position', [0.9075 0.8 0.09 0.03],...
    'Callback', @dontwant);
     set(hdel,'BackgroundColor',get(gcf,'Color')); % checkbox otherwise sits in a white rectangle, which looks bad

if isfield(cfg,'lfp')
    % spectrogram button
    uicontrol('Style', 'pushbutton', 'String', 'Spectrogram',...
        'Units','normalized','Position', [0.01 0.88 0.1 0.05],...
        'FontUnits','normalized','Callback', @spectraxis);
    
    % spectrogram z scale drop down option
    zPop = uicontrol('Style','popupmenu','String',{'root','decibel-watt','raw'},...
        'Units','normalized','Position',[0.01 0.76 0.1 0.05]);
    
    % spectrogram frequency range 
    frange(1) = uicontrol('Style','edit','String','50',...
        'Units','normalized','Position',[0.01 0.82 0.04 0.04]);
    frange(2) = uicontrol('Style','edit','String','300',...
        'Units','normalized','Position',[0.07 0.82 0.04 0.04]);
    
    % create exes for a spectrogram. These do not move with navigate, unlike the main
    % axes, and instead are set to invisible unless a spectrogram is plotted
    ax_spec = axes('Position', get(gca, 'Position'),'Visible','off');
    set(hfig,'CurrentAxes',ax_main)
end

%% this needs to be called AFTER the uicontrols because it makes reference to some of them

PlotCurrentIV

%% helper functions

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
            set(hnote,'String','[enter text]')
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
        h_curr(1) = plot([iv_in.tstart(current_iv) iv_in.tstart(current_iv)],ylims,'Color',cfg.ivColor);
        h_curr(2) = plot([iv_in.tend(current_iv) iv_in.tend(current_iv)],ylims,'Color',cfg.ivColor);
        title(['Annotate event ',num2str(current_iv),' out of ',num2str(length(iv_in.tstart))],'FontSize',14)
        drawnow % this makes it wait for the plotting to happen before it updates the buttons (plotting very slow)
        updateButtons
    end

% ~~~~~~ HIDE SPECTRAXIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function HideSpectraxis
        mafak = findobj(ax_spec); % use findobj, otherwise can't make imagesc invisible (just axes)
        set(mafak,'Visible','off'); set(ax_main,'Color','w')
    end

%% callback functions

% ~~~~~~ NEXT INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function nextIV(~,~)
        HideSpectraxis
        if ii == length(iv_in.tstart) % we've gone through all of them
            set(hnext,'Enable','off')
            title('All intervals have been annotated','FontSize',14)
        else
            ii = ii+1;
            current_iv = idx_list(ii);
            PlotCurrentIV
        end  
        set(hprev,'Enable','on')  
    end

% ~~~~~~ PREVIOUS INTERVAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function prevIV(~,~)
        HideSpectraxis
        if ii-1 <= 0 % we're at the beginning
            set(hprev,'Enable','off')
        else
            ii = ii-1;
            current_iv = idx_list(ii);
            PlotCurrentIV
        end
        set(hnext,'Enable','on')       
    end

% ~~~~~~ KEY STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function keystuff(source,event)
        if strcmp(event.Key,'n')
            nextIV(source,event)
        else
            HideSpectraxis
            navigate(source,event)
        end
    end

% ~~~~~~ WRITE STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function writestuff(source,~)
        % add annotation to iv
        str = strtrim(get(source, 'String')); 
        % remove leading and trailing whitespace
        % before str = '    '; now str = ''; or before str = ' 5 '; now str = '5';
        iv_in.usr.annotation{current_iv} = str;
        
        updateProgress
    end

% ~~~~~~ CLICK STUFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function clickstuff(~,~)
        % handle mouse clicks
        HideSpectraxis
        
        % right now this is just for hiding spectrogram axis, but can add
        % more mouse things here
        
    end

% ~~~~~~ DON'T WANT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function dontwant(source,~)
        % keeps a record of which intervals to delete
        del(current_iv) = get(source, 'Value');
        
        updateProgress
    end

% ~~~~~~ LEAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function leaveme(~,~)
        % quit request function
        temp = iv_in.usr.annotation;
        temp(logical(del)) = [];
        nEmpties = length(temp(cellfun(@isempty,temp)));
        
        if nEmpties > 0
            choice = questdlg(['You are missing ',num2str(nEmpties),' annotations. Are you sure you want to close the figure? Any unsaved data will be lost.'],'Quit Requested','Quit anyway','Cancel','Cancel');
            switch choice
                case 'Cancel'
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
                    return
            end
        end
    end % of leaveme

% ~~~~~~ SAVE ME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function saveme(~,~)
        % save request function
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
                
            case 'No'
                return
        end
    end % of saveme

% ~~~~~~ SPECTRAXIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function spectraxis(~,~)
        % create a new set of axes and plot a spectrogram behind the raster plot
        set(gca,'color','none')
        xlims = get(gca,'XLim'); 
        xticks = get(gca,'XTick');
        set(hfig,'CurrentAxes',ax_spec) % now change current axes
        ugh = restrict2(cfg.lfp,xlims(1),xlims(2));
        
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
                col = [0.09*10^-6 18*10^-6]; % some arbitrary range for the color scaling, keep everything relative or w/e
            case 'decibel-watt'
                P = 10*log10(P);
                col = [-200 -80];
            case 'raw'
                col = [0.09*10^-10 3*10^-10];
        end
        
        imagesc(T,F,P,col); 
        set(ax_spec,'YAxisLocation','right','XTick',xticks,'XTickLabel',[],'YDir','normal','Visible','on')
        
        ylabel(ax_spec,'Frequency (Hz)')
        grid on
        uistack(ax_spec,'down') % send backwards
        % return to main axes (for navigation and stuff)
        set(hfig,'CurrentAxes',ax_main) 
    end % of spectraxis

end

