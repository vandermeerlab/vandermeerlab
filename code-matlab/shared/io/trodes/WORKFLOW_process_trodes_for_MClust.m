%% define where stuff lives on your machine
params.trodesDir = 'C:\Trodes_2-3-2_Windows64';
params.dataDir = 'C:\data\r204_screening_rec16';
params.githubDir = 'C:\Users\mvdm\Documents\GitHub';

%% set paths, discover filenames
restoredefaultpath;
addpath(genpath(cat(2, params.githubDir, '\vandermeerlab\code-matlab\shared')));
addpath(cat(2, params.trodesDir, '\Resources\TrodesToMatlab\'));
addpath(params.trodesDir);

full_fn = FindFile('*.rec', 'StartingDirectory', params.dataDir);
[fp, fn, fe] = fileparts(full_fn);
cd(params.dataDir);

%% call trodes utilities to export data from .rec file
% NOTE: default uses reference channel defined in recording environment; to
% re-reference, do ???
% NOTE: spike detection parameters defined in recording environment; to
% re-reference, do ???

eval(cat(2, '! ', params.trodesDir, filesep, 'exportspikes.exe -rec ', fn, fe, ' -output ', fn));
eval(cat(2, '! ', params.trodesDir, filesep, 'exportLFP.exe -rec ', fn, fe, ' -output ', fn));

%% now batch-run KlustaKwik
restoredefaultpath;
addpath(genpath(cat(2, params.githubDir, '\vandermeerlab\code-matlab\toolboxes\Mclust-3.5')));
addpath(cat(2, params.trodesDir, '\Resources\TrodesToMatlab')); % NOTE make sure that you rename Trodes's MATLAB loading engine so that MClust's one takes precedence

eval(cat(2, '! copy ', which('Batch_Trodes.txt'), ' ', fn, '.spikes\Batch.txt')); % copy basic Trodes template

pushdir(cat(2, fn, '.spikes'));
RunClustBatch;

%% ready to sort in MClust; this will generate .t files

%% example to load and plot data
cd('C:\Data\r204_screening_rec16\r204_screening_rec16.LFP')

out = readTrodesExtractedDataFile('r204_screening_rec16.LFP_nt11ch1.dat');
out_ts = readTrodesExtractedDataFile('r204_screening_rec16.timestamps.dat');
lfp(1) = tsd(double(out_ts.fields.data)./out_ts.clockrate, double(out.fields.data)'); % ugly and hard to test, should make wrapped loading function

out = readTrodesExtractedDataFile('r204_screening_rec16.LFP_nt24ch1.dat');
lfp(2) = tsd(double(out_ts.fields.data)./out_ts.clockrate, double(out.fields.data)');
clear out out_ts;

% spikes
cd('C:\Data\r204_screening_rec16\r204_screening_rec16.spikes')

please = [];
please.uint = '64';
S = LoadSpikes(please);

% plot
cfg_mr = [];
cfg_mr.lfp = lfp;
MultiRaster(cfg_mr, S);