function ncfs = SWRfreak(cfg_in,SWRtimes,csc)
%SWRFREAK get frequency spectrum for manually identified sharp wave-ripples
%
%   ncfs = SWRfreak(SWRtimes,csc,)
%
%   ncfs = SWRfreak(SWRtimes,csc,NOISEtimes) % not implemented yet
%
%   cfg.weightby = 'power'; or 'amplitude'
%                  'power' re-weight the spectrum to "unbias" the voltage
%                  over each frequency 
%                  'amplitude' raw spectrum
%   cfg.hiPassCutoff = 100; in Hz, disregard all frequencies below
%   cfg.fs = 2000; in Hz, the sampling frequency
%   cfg.win1 = 0.06; in s, the window size
%   cfg.win2 = []; specify if you want two windows to be used
%   Suggested settings for two windows:
%         cfg.win1 = 0.08; in s, the "noise reduction" window
%         cfg.win2 = 0.04; in s, the "precision" window
% 
%   cfg.showfig = 0; 0 do not display the figures; 1 display the figures
%          There is a hidden fig config if you open the function and read
%          the section called "Parse cfg parameters".
%   cfg.openNewFig = 1; if using SWRfreak with subplot, set this to 0
%   cfg.verbose = 1; If 1 talk to me, if 0 don't
%
%   OUTPUT
%   
%   ncfs - noise-corrected frequency spectrum: Fourier coefficients for 
%              SWRs minus the Fourier coefficients for HC background noise
%         .freqs1      - frequency spectrum based on win1
%         .freqs2      - (if win2 specified)frequency spectrum based on win2
%         .label       - csc used
%         .cfg         - record of cfg history
%         .parameters  - all parameters used to generate the output
% 
%   For additional information, open SWRfreak and read the ABOUT section
%
% Elyot Grant, Jan 2015 (math and original code)
% ACarey, Jan 2015 (responsible for poor function design >_<)
% - edit Mar 2015 -- additional config options, plotting

%% ABOUT SWRFREAK

% Why are there two windows, and why are 40 ms and 80 ms the defaults?

% The noise reduction and precision windows were chosen based on how well
% the plotted scores (from amSWR) characterized the data in the csc and the 
% spiketrains.
% 120, 100, 80, 60, and 40 ms windows were used for FFTing the data. 120
% and 100 merged nearby events and were rejected, despite lower rates
% of false positives. 80 ms still merged nearby events, but to
% a lesser degree than 100 and 120. It also had decently low scores for false
% positives. 40 ms had the best separation for nearby events, but was the
% worst for peaking at false positives. The geometric mean of the two of
% them did a good job, so it was decided that two sets of ncfs
% were better than one. 
% But then it was re-decided that a single 60 ms window was fine. It's also
% faster for amSWR. Just use 60.

%% Parse cfg parameters

cfg_def.verbose = 1;
cfg_def.weightby = 'power'; 
cfg_def.win1 = 0.06;
cfg_def.win2 = [];
cfg_def.hiPassCutoff = 100; %We want to delete all frequencies below 100 Hz
cfg_def.fs = 2000; 
cfg_def.outputSWRfreqs = 0; 
cfg_def.outputNOISEfreqs = 0;
cfg_def.showfig = 0;
cfg_def.openNewFig = 1; 
    % this makes it easier to delete the fig config from history:
    cfg_def.fig.SWRcolor = 'k'; % line colour for raw SWR freqs
    cfg_def.fig.NOISEcolor = [0.67 0.67 0.67]; % line colour for noise freqs 
    cfg_def.fig.FREQScolor = 'r'; % line colour for noise-corrected SWR freqs
    cfg_def.fig.LineWidth = 2; 
    cfg_def.fig.xlim = [0 600]; % xvals go all the way up to 1000 Hz (nyquist), but the freqs are basically flat after 600 Hz
    cfg_def.fig.TitleFontSize = 14;
    cfg_def.fig.LegendFontSize = 12;
    cfg_def.fig.xyFontSize = 11;
    cfg_def.fig.LegendBoxAspectRatio = []; % [xsize ysize zsize]
    if isfield(cfg_in,'fig')
        fig = ProcessConfig2(cfg_def.fig,cfg_in.fig);
    end
    
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if isfield(cfg_in,'fig')
    cfg.fig = fig;
end
%% check if csc is the same as the one used in ducktrap

if ~strcmp(SWRtimes.label{1}(1:15),csc.label{1}(1:15))
    error('CSC must be from the same day as the one used for manual identification')
end

if ~strcmp(SWRtimes.label,csc.label)
    warning('CSC is different from the one used for manual SWR identification')
end

%% internal function to do the thing

    function freqs = freakHelper(SWRtimes,csc,timewin)
        sampwin = timewin*cfg.fs; % the window size in nSamples = timewindow * sampling frequency 
        %tstart = tcent - timewin/2;
        %tend = tcent + timewin/2;
        tcent = IVcenters(SWRtimes);
        midIndices = nearest_idx3(tcent,csc.tvec);
        SWRsum = 0;
        for iSWR = 1:length(tcent) 
            nextFFT = windowedFFT(cfg,csc.data,sampwin,midIndices(iSWR));
            SWRsum = SWRsum + nextFFT;
        end

        midIndices = nearest_idx3(tcent+2,csc.tvec); % Add 2 seconds to each clicked time -> random time (this could error if a clicked time was closer than 2 s away from the end of recording?)
        SWRsumNoise = 0;
        for iSWR = 1:length(tcent)
            nextFFT = windowedFFT(cfg,csc.data,sampwin,midIndices(iSWR));
            SWRsumNoise = SWRsumNoise + nextFFT;
        end

        SWRsum = conv(SWRsum,[0.1,0.2,0.4,0.2,0.1]); % narrow smoothing kernel
        SWRsum = SWRsum(3:length(SWRsum)-2);
        SWRsum = SWRsum./sum(SWRsum); %Normalize
        %plot(SWRsum);

        SWRnoise = conv(SWRsumNoise,[0.1,0.2,0.4,0.2,0.1]);
        SWRnoise = SWRnoise(3:length(SWRnoise)-2);
        SWRnoise = SWRnoise./sum(SWRnoise);
        %figure;plot(SWRnoise);

        freqs = (SWRsum - SWRnoise);
        freqs = conv(freqs,[0.1,0.2,0.4,0.2,0.1]);
        freqs = freqs(3:length(freqs)-2);
        %freqs = max(freqs,zeros(size(freqs))); %%%%%%%%%%%%%%%%%%%%%%
        %figure;plot(freqs);
        
        if cfg.showfig == 1 
            % plot some stuff
            
            frequency = (1:timewin*1000)./timewin;
            
            maxSWR = max(SWRsum);
            minSWR = min(SWRsum);
            maxNOISE = max(SWRnoise);
            minNOISE = min(SWRnoise);
            maxDIFF = max(freqs);
            minDIFF = min(freqs);
            
            ymax = max([maxSWR maxNOISE maxDIFF]);
            ymin = min([minSWR minNOISE minDIFF]);
            
            if cfg.openNewFig
                figure; 
            end
            hold on;
            plot(frequency,SWRsum,'Color',cfg.fig.SWRcolor,'LineWidth',cfg.fig.LineWidth)
            plot(frequency,SWRnoise,'Color',cfg.fig.NOISEcolor,'LineWidth',cfg.fig.LineWidth)
            plot(frequency,freqs,'Color',cfg.fig.FREQScolor,'LineWidth',cfg.fig.LineWidth)
            plot([cfg.fig.xlim(1) cfg.fig.xlim(2)],[0 0],'Color','k','LineStyle',':','LineWidth',1)
            
            xlabel('Frequency (Hz)','Fontsize',cfg.fig.xyFontSize); 
            if strcmp(cfg.weightby,'power')
                type = 'Power weighted, ';
                y_label = '';
            elseif strcmp(cfg.weightby,'amplitude')
                type = 'Amplitude weighted, ';
                y_label = 'Fourier coefficients';
            end
            xlim(cfg.fig.xlim); ylim([ymin-0.01 ymax+0.01])
            title([type,sprintf('%d SWRs, %d ms window',length(tcent),timewin*1000)],'FontSize',cfg.fig.TitleFontSize);
            ylabel(y_label,'Fontsize',cfg.fig.xyFontSize);
            box on
            hL = legend('SWR','background noise','SWR - noise');
            if ~isempty(cfg.fig.LegendBoxAspectRatio)
                set(hL,'PlotBoxAspectRatio',cfg.fig.LegendBoxAspectRatio)
            end
            set(hL,'FontSize',cfg.fig.LegendFontSize)
        end
    end

%% generate the output

if ~isempty(cfg.win2)
    freqs1 = freakHelper(SWRtimes,csc,cfg.win1);
    freqs2 = freakHelper(SWRtimes,csc,cfg.win2);
else 
    freqs1 = freakHelper(SWRtimes,csc,cfg.win1);
    freqs2 = [];
end

%% return the output
parameters = struct('weightby',cfg.weightby,'win1',cfg.win1,'win2',cfg.win2','fs',cfg.fs,'hiPassCutoff',cfg.hiPassCutoff,'csc',csc.label);

ncfs = struct('freqs1',freqs1,'freqs2',freqs2);

ncfs.parameters = parameters;
ncfs.label = csc.label;

% keep a record
cfg = rmfield(cfg,'fig'); % because who cares what the figure settings were?
ncfs.cfg.history.mfun = mfilename;
ncfs.cfg.history.cfg = {cfg};

end

