%% config
params.trodesDir = 'C:\Trodes_2-3-2_Windows64';
params.dataDir = 'C:\Data\r204_screening_rec16\r204_screening_rec16.spikes'; % NOTE assumes that spikes have been extracted (use WORKFLOW_process_trodes_for_MClust.m)
params.githubDir = 'C:\Users\mvdm\Documents\GitHub';

params.plotstring = '.k';
params.MarkerSize = 0.5;
params.FontSize = 6;

%% set up path
restoredefaultpath;
addpath(genpath(cat(2, params.githubDir, '\vandermeerlab\code-matlab\shared')));
addpath(cat(2, params.trodesDir, '\Resources\TrodesToMatlab\'));
%addpath(params.trodesDir);

%% some features - these will take [nSpikes x 40] waveform data matrix as input
feature_peak_fun = @(x) max(x, [], 2);
feature_energy_fun = @(x) sum(x.^2, 2);

%%
cd(params.dataDir);
fd = FindFiles('*.dat');


for iF = 1:length(fd) % one figure per tetrode

    data = readTrodesExtractedDataFile(fd{iF});

    [fp fn fe] = fileparts(fd{iF});
    
    clear feature_peak feature_energy;
    for iCh = 4:-1:1
       
        feature_peak(iCh, :) = feature_peak_fun(data.fields(iCh + 1).data);
        feature_energy(iCh, :) = feature_energy_fun(data.fields(iCh + 1).data);
        
    end
    
    figure(iF);
    
    subplot(231)
    plot3(feature_peak(1,:), feature_peak(2,:), feature_peak(3,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('peak 1'); ylabel('peak 2'); zlabel('peak 3');
    
    subplot(232)
    plot3(feature_peak(2,:), feature_peak(3,:), feature_peak(4,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('peak 2'); ylabel('peak 3'); zlabel('peak 4');
    
    subplot(233)
    plot3(feature_peak(3,:), feature_peak(4,:), feature_peak(1,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('peak 3'); ylabel('peak 4'); zlabel('peak 1');
    
    subplot(234)
    plot3(feature_energy(1,:), feature_energy(2,:), feature_energy(3,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('energy 1'); ylabel('energy 2'); zlabel('energy 3');
    
    subplot(235)
    plot3(feature_energy(2,:), feature_energy(3,:), feature_energy(4,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('energy 2'); ylabel('energy 3'); zlabel('energy 4');
    
    subplot(236)
    plot3(feature_energy(3,:), feature_energy(4,:), feature_energy(1,:), params.plotstring, 'MarkerSize', params.MarkerSize);
    grid on; set(gca, 'FontSize', params.FontSize);
    set(gca, 'XTick', '', 'YTick', '', 'ZTick', '')
    xlabel('energy 3'); ylabel('energy 4'); zlabel('energy 1');
    
    drawnow;
    
    f_out = cat(2, fn, '_scatter.png');
    WriteFig(gcf, f_out, 1);
    
    close(gcf);

end