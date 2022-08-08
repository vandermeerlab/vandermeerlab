function AMPX_phase_distrib_plot(phase, trial_data, data, varargin)
%AMPX_pow_distrib_plot
%  Plot the phase distributions from AMPX_phase_distrib
%
%  This will create a plot of the phase differences for each reference
%  channel speficied in the input (default: all channels).
%
%  Saving is optional
%
%
% EC v1.0 12/09/2014
%% extract and define the variables
run(FindFile('*Keys.m'));
save_fig = 'on';
min_freq = phase.params.freq(1);
max_freq = phase.params.freq(2);
trials_def = 1:length(trial_data.trial);
ref_channels = ExpKeys.Ref_chan;
extract_varargin

if min_freq<65; type = 'low_gamma'; clim = [0 900]; type_title = 'low gamma';
elseif min_freq>66 && max_freq<100; type = 'high_gamma'; type_title = 'high gamma'; clim = [10 300];
elseif 100<=min_freq; type = 'high_freq'; type_title = 'high freq with randomized iei'; clim = [0.5 2.5];
end

% make an 8x8 matrix of the channel labels.  used for putting text on the
% figure
chan_map = reshape(ExpKeys.Probe_layout, 8, 8);
%% Make sure that the folder exists
if exist(['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' type '_phase_figures'],'dir')==0
    mkdir(['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' type '_phase_figures'])
end
    


%%  Create the plots

for itrial = trials_def
    for ref = ref_channels
        [r,c] =find(chan_map==ref);
        phase_dist = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,4,1:2)
        imagesc(phase.phase_distrib{itrial}{ref}, [-.15 .15])
        title(['Phase difference for trial:' num2str(itrial) ' (ref: ' num2str(ref) ')'], 'fontsize', 14)
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
        text(c-.2,r,num2str(ref),'FontSize',16)
        colorbar('location','southoutside');
        
        subplot(1,4,3:4)
        imagesc(phase.phase_distrib_avg{ref}, [-.15 .15])
        title(['Average phase difference (ref: ' num2str(ref) ')'], 'fontsize', 14)
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
        text(c-.2,r,num2str(ref),'FontSize',16)
        colorbar('location','southoutside');
        
        %save if requested.
        if strcmp(save_fig,'on'); print(gcf,'-dpng','-r300',['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' date '\' type '_phase_figures\' data.hdr.Filename(5:end-4) '_' num2str(itrial) '_' num2str(ref) '_spec_test.png']);end
        close(phase_dist)
    end
end

end

