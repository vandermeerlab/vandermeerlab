function out = CategorizeStriatumWave(cfg_in,S)
% function out = CategorizeStriatumWave(cfg,S)
%
% striatal data cell type classification based on waveform and firing rate
% following Berke et al. Neuron 2004
%
% NOTE: *wv and *ClusterQual files MUST be present (use Create_CQ_File.m)
%
% INPUTS:
%
% S: ts struct with spike data
%
% CONFIGS:
%
% cfg.MinID = 20; % minimum isolation distance for cell to be considered
% cfg.MaxLr = 0.1; % maximum L-ratio for cell to be considered
% cfg.MinAmpl = 200; % minimum peak (max) amplitude for cell to be considered
% cfg.cMethod = 'Berke'; % selects set of classification criteria to use
%   (alternate: 'MvdM', i.e. as in van der Meer & Redish, Front Int Nsci 2009)
%
% OUTPUTS:
%
% out.msn: idx into S.t for cells classified as putative medium spiny
%   neuron
% out.fsi: idx into S.t for cells classified as putative fast-spiking
%   interneuron
% out.other: idx into S.t for cells classified as other (neither MSN not
%   FSI)
% out.wv: waveforms for each neuron (taken directly from *wv files)
% out.cq: cluster quality metrics for each neuron (taken directly from *cq
%   files)
% out.fr: firing rate (mean(1/ISI)) for each neuron
% out.pw: waveform peak width (ms) for each neuron
% out.vw: waveform valley width (ms) for each neuro
%
% MvdM 2015-02-05 

cfg_def = [];
cfg_def.MinID = 20;
cfg_def.MaxLr = 0.1;
cfg_def.MinAmpl = 600;
cfg_def.cMethod = 'Berke'; % alt: 'MvdM'

cfg = ProcessConfig2(cfg_def,cfg_in);


nCells = length(S.t);
for iC = nCells:-1:1

    cellfname = S.label{iC};

    clufname = regexprep(cellfname,'_','-');
    [fp clufname] = fileparts(clufname);
    clufname = cat(2,clufname,'-wv');
    
    try
        eval(['load(''' clufname ''');']);
    catch
        try
            clufname = regexprep(clufname,'TT0','TT');
            eval(['load(''' clufname ''');']);
        catch
            error('Unable to load wv file %s',clufname);
        end
    end

    wvfname = regexprep(cellfname,'_','-');
    [fp wvfname] = fileparts(wvfname);
    wvfname = cat(2,wvfname,'-ClusterQual');
    
    try
        eval(['load(''' wvfname ''');']);
    catch
        try
            wvfname = regexprep(wvfname,'TT0','TT');
            eval(['load(''' wvfname ''');']);
        catch
            error('Unable to load CQ file %s',wvfname);
        end
    end

    [pw(iC),vw(iC)] = waveWidths(mWV);
    wv(iC,:) = mWV(:);
    cq.ampl(iC) = max(mWV(:));
    
    % could consider adding alternative ways to compute average firing rate
    
    %fr(iC) = length(Data(sd.S{iC}))/(sd.ExpKeys.TimeOffTrack-sd.ExpKeys.TimeOnTrack);
    fr(iC) = 1./mean(diff(S.t{iC}));
    
    cq.id(iC) = CluSep.IsolationDist; if isnan(cq.id(iC)), cq.id(iC) = -Inf; end 
    cq.lr(iC) = CluSep.Lratio;

end


switch cfg.cMethod
    case 'Berke'

        % Berke et al. classification
        fsi = find(vw < 0.265 & fr > 2 & pw < 0.12 & cq.id >= cfg.MinID & cq.lr <= cfg.MaxLr & cq.ampl >= cfg.MinAmpl);
        msn = find(vw > 0.300 & fr < 5 & cq.id >= cfg.MinID & cq.lr <= cfg.MaxLr & cq.ampl >= cfg.MinAmpl);
        other = setdiff(1:nCells,fsi);
        other = setdiff(other,msn);

    case 'MvdM'

        % alternative classification that fits our data a little better...
        fsi = find(vw < 0.35 & fr > 2 & pw < 0.15 & cq.id >= cfg.MinID & cq.lr <= cfg.MaxLr & cq.ampl >= cfg.MinAmpl);
        msn = find(vw >= 0.35 & fr < 5 & cq.id >= cfg.MinID & cq.lr <= cfg.MaxLr & cq.ampl >= cfg.MinAmpl);
        other = setdiff(1:nCells,fsi);
        other = setdiff(other,msn);

end

out.fsi = fsi;
out.msn = msn;
out.other = other;
out.wv = wv;
out.cq = cq;
out.fr = fr;
out.pw = pw;
out.vw = vw;