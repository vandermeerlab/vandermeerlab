function SWRtimes = ducktrap(cfg_in,S,csc)
% DUCKTRAP Manually identify windows containing sharp wave-ripple
% complexes. The default window is 0.08 seconds (80 milliseconds).
%
% function SWRtimes = ducktrap(cfg,S,csc) plots S (spiketrain) and csc (LFP) 
% and allows the user to manually identify SWRs by clicking on the figure to 
% produce a containment window (trapwin). Upon clicking, trapwin will be
% plotted on the figure as a blue horizontal line. To keep the SWR time,
% press enter. Trapwin will change to green and remain on the plot. The
% function will return output when the user collects the number of SWR
% times specified by cfg.num.
%
% User can move viewing window via keyboard input in the same general manner as 
% MultiRaster (type "help navigate" in command window). *** some navigate
% features may not compatible with ducktrap.
%
% For more information, open ducktrap and the read the section "About ducktrap".
%
%
%   INPUTS:
%       cfg_in: input cfg
%       S: spiketrain struct. Output from LoadSpikes.
%       csc: tsd object (LFP). Output from LoadCSC.
%
%   OUTPUT
%       SWRtimes: iv struct with fields:
%           .tcent   - actual click times; this is the important one
%           .tstart  - [num x 1 double] start times
%           .tend    - [num x 1 double] end times
%           .trapwin - window width used to "catch" SWRs.
%           .label   - the csc used for identification
%
%       ducktrap_backup - [1 x numKept double]; transient variable
%              containing x-axis values clicked by the user. It exists in 
%              the base workspace while the function is running, and 
%              remains only if the function does not reach completion.  
%
%   CFG OPTIONS:
%
%       cfg.trapwin - default 0.08 seconds (80 milliseconds). The window for
%           containing individual SWRs manually identified by the user. It
%           is the same window that will be used in SWRfreak and detectSWR.
%
%       cfg.num - default 50. The number of SWRs you would like to
%           manually identify.
%
%       cfg.evt - default NONE
%           Event data in iv format. The output from getSWR or a similar function.
%           Useful for viewing "loud" SWRs with few false positives if the
%           detection threshold was properly set; may help identify
%           positives. 
%
%       cfg.lfpHeight - default 15 
%           The height of the lfp (csc) from maximum to minium in y-axis units. Vertically 
%           stretches the lfp to improve visibility of oscillations.  
%
%       cfg.lfpMax - default 15
%           Values greater than lfpMax times above mean(abs(lfp.data)) will be plotted 
%           as NaNs. Improves visibilty by cutting off high-amplitude [noise] blips 
%           in raw lfp.
%
% ACarey, Jan 2015. 
% plotting code modified from youkitan's MultiRaster 

%% About ducktrap 

% Ducktrap was designed to be used as the first step in identifying SWRs
% using the discrete Fourier transform. The identified SWRs are used by
% the next function to find Fourier cofficients for SWRs. 
%
% SETTING UP DUCKTRAP:
%
% please = cfg.useClustersFile = 0;
% S = LoadSpikes(please);
% cfg.fc = {'R050-2014-04-02-CSC07a.ncs};
% csc = LoadCSC(cfg);
% evt_getSWR = getSWR(cfg,S,csc);
% cfg.evt = evt_getSWR;
% SWRtimes = ducktrap(cfg,S,csc);

% HOW TO USE DUCKTRAP:
%
% SWRtimes = ducktrap(cfg,S,csc);
%
% A figure will open that displays the entire csc and S.
% Use keyboard input to navigate to different areas of the plot ("help
% navigate"). If you've plotted evt produced by power thresholding, this
% can help you identify positives. Without evt, nagivate will move the
% window in units of seconds using the keys A and D. 
% When you see a SWR you want to keep, click on its centerpoint and trapwin
% will appear on the screen as a horizontal blue bar. If the trapwin
% seems properly placed, type Enter and the trapwin will change
% to green. If you do not want to keep the selection, then press another 
% key or click elsewhere. Once you've hit Enter, there's no way to undo the
% choice, so don't mess up :P The figure legend keeps track of how many
% SWR times you've collected. If you accidentally close the figure, you can
% access your previous clicks in ducktrap_backup, but you need to -/+
% trapwin/2 to get the interval.
%
% Ducktrap's output can be sent into SWRfreak to get the Fourier
% coefficients of the sum of your manually identified SWRs. 
%
% WHY DO I HAVE TO SPECIFY CFG.NUM?
%
% Because AC sucks at programming, and couldn't get the function to work
% unless there was some way to terminate the while loop cleanly. If you
% don't know how many SWRs you want to get, then specify a large num
% and close the figure when you're done. This will generate errors. But you
% can still get the SWR center times from ducktrap_backup in the base workspace.
%
% WHY IS IT CALLED DUCKTRAP?
% It's a joke based on an analogy that EGrant gave me about power 
% thresholding finding only the loudest SWRs; it's like recording from ducks
% but only being able to detect the loudest quacks... 
%
%%
% DUCKTRAP: THE STORY
%
% There once was a lowly biologist who liked to stalk ducks. She recorded 
% the sounds of the ducks, which were neighbours to the geese and the noisy
% waterfall. Later, she wanted to count how many times the ducks quacked, 
% but more often than not her detection would only find the loudest ducks, 
% and sometimes it would also detect the pesky geese (that incidentally also 
% liked to shit everywhere with no regard for public decency). So, the lowly 
% biologist complained to the arrogant mathematician in passing.
% 
% And he told her, “You biologists are so dumb, you’re probably doing 
% something like filtering the signal and looking for high power in a single 
% ripple band. If you’re just thresholding over a single frequency band you’ll 
% miss all the quiet ducks and detect all the loud geese as false positives”.
% 
% And her response was something along the lines of, “I know: I’m a lowly 
% biologist. I’m worse than lowly. In fact, I’m a plebian—no, a platyhelminth! 
% An intestinal parasite!”
% 
% Somehow this activated the arrogant mathematician’s math superiority complex, 
% “Look, if you want to find the quiet ducks, too, here’s how you should do it. 
% First, find me at least 100 real ducks.”
% 
% Supremely grateful (but also a smidgen annoyed by the arrogant mathematician’s 
% arrogance), the lowly [intestinal] biologist set out to rob ducks of
% their freedom.
%
%% PLOTTING SECTION
% duckplotter is based on MultiRaster, but is probably more prone to errors
    function trapwin = duckplotter(cfg_in,S,csc)    
        %% Set cfg parameters and check inputs

        cfg_def.trapwin = 0.08;
        cfg_def.num = 50;
        cfg_def.SpikeHeight = 0.4;
        cfg_def.axisflag = 'tight';
        cfg_def.spkColor = 'k';
        cfg_def.lfpColor = 'r';
        cfg_def.lfpHeight = 15;
        cfg_def.lfpMax = 15;
        cfg_def.axislabel = 'on';
        cfg_def.windowSize = 1;
        cfg = ProcessConfig2(cfg_def,cfg_in); 
        
        % send trapwin back out
        trapwin = cfg.trapwin;

        %% Setup navigate

        global evtTimes windowSize time 
        windowSize = cfg.windowSize;
        time = csc.tvec;

        % Initialize events
        if isfield(cfg,'evtTimes')
            evtTimes = cfg.evtTimes;
        %elseif ~isfield(cfg.evt,'tstart')
        else
            % navigate evt by each second of recording?
            evtTimes = time(1:2000:end); % since length(time) is nSamples, step up by fs to get seconds? 
        end
    
        figure('KeyPressFcn',@navigate);
        hold on;
%assignin('base','time',time)
        %% Decide plotMode

        if isfield(cfg,'evt') % add evt to plot
            plotMode = 2;

        else %default S and csc
            plotMode = 1;

        end

        %% Plot
        switch plotMode
            case 1 % S + csc
                PlotSpikeRaster2(cfg,S); hold on;
                abscsc = abs(csc.data);
                nans_here = abscsc > cfg.lfpMax*mean(abscsc);
                csc.data(nans_here) = NaN;
                csc.data = rescale(csc.data,-cfg.lfpHeight,0);
                plot(csc.tvec,csc.data,'k')
                %ylims(1) = -cfg.lfpHeight - 1;
       
            case 2 % S + csc + evt
                evtTimes = (cfg.evt.tstart + cfg.evt.tend)./2;
                PlotSpikeRaster2(cfg,S); hold on;
                abscsc = abs(csc.data);
                nans_here = abscsc > cfg.lfpMax*mean(abscsc);
                csc.data(nans_here) = NaN;
                csc.data = rescale(csc.data,-cfg.lfpHeight,0);
                cfg_temp.display = 'tsd';
        
                PlotTSDfromIV(cfg_temp,cfg.evt,csc);
                %ylims(1) = -cfg.lfpHeight - 1;
        end

        %% Set Axis limits
        xlims = get(gca,'XLim');
        ylims = get(gca,'Ylim');

        xlim([xlims(1) xlims(2)]);
        ylim([ylims(1) ylims(2)]);

        hold off;
      
    end

%% DUCKTRAP SECTION
trapwin = duckplotter(cfg_in,S,csc);
hold on;

numKept = 0;
while numKept < cfg.num
    
    [x,y] = ginput_ax(gca,1);
    x_window = [x-trapwin/2 x+trapwin/2];
    H(1) = plot(x_window,[y y],'o','Color',[30/255 144/255 1],'MarkerSize',8); H(2) = plot(x_window,[y y],'LineWidth',2,'Color',[30/255 144/255 1]);

    keydown = waitforbuttonpress;
    
    if keydown
        key = get(gcf, 'CurrentCharacter');
 
        if key == 13; % the code for enter;
            set(H(1),'Color','g') 
            set(H(2),'Color','g') 
            numKept = numKept + 1;
            clicktimes(numKept) = x;
            legmess = [num2str(numKept),'/',num2str(cfg.num)];
            legend(legmess);
            assignin('base','ducktrap_backup',clicktimes);
        
        else
            if exist('H','var')
                delete(H);
            end
        end
    
    else
        delete(H); 
    end  
end

clicktimes = sort(clicktimes);
SWRtimes = iv(clicktimes-trapwin/2,clicktimes+trapwin/2);
SWRtimes.tcent = clicktimes';
SWRtimes.trapwin = trapwin;
SWRtimes.label = csc.label;

evalin('base','clear ducktrap_backup')
end