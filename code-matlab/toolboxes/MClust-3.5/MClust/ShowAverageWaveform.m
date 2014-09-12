function ShowAverageWaveform(iClust,varargin)

% ShowWaveformDensity(iClust)
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% ADR 2003
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
% Extensively modified by ADR to accomodate new ClusterOptions methodology
%
% MvdM 2013-08-16 added PlotHere varargin to use with
% PlotAllAverageWaveformsS

        global MClust_TTData MClust_Clusters MClust_FeatureData
        global MClust_Colors
        
        run_avg_waveform = 1;
        PlotHere = 0;
        
        extract_varargin;

        [f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
        

        if length(f) == 0
            run_avg_waveform = 0;
            msgbox('No points in cluster.')
        end
 
        if run_avg_waveform
            clustTT = ExtractCluster(MClust_TTData, f);
            [mWV sWV] = AverageWaveform(clustTT);
            nWVSamples = size(mWV,2);
            
            if PlotHere
                AveWVFig = gca;
            else
                AveWVFig = figure;
            end
            
            for it = 1:4
                xrange = ((nWVSamples + 2) * (it-1)) + (1:nWVSamples);
                
                if ~PlotHere, figure(AveWVFig); end
                
                hold on;

                plot(xrange, mWV(it,:) + sWV(it,:),':w','LineWidth',2);
                plot(xrange, mWV(it,:) - sWV(it,:),':w','LineWidth',2);
                h = plot(xrange, mWV(it,:),'LineWidth',2);

                set(h,'Color',MClust_Colors(iClust + 1,:));
            end
            
            axis off
            axis([0 4*(nWVSamples+2) min([-21000 min(mWV(:) - sWV(:))]) max([21000 max(mWV(:) + sWV(:))])])
       
            title(['Average Waveform: Cluster ' num2str(iClust)]);
        end
