function spike_images(fname, varargin)

% Spike_images will acquire screen shots (in groups of 8) of each channel
% during a given epoch.  
% Inputs:
%
% fname: this is the name of the file to be used.  
%
% type: this will dictate what epoch will be used.  
%    'pre' -  will use the entire pre-record session
%    'post' -  will use the entire post-record session
%    'maim' [default] - will result in a 

%% Identify the variables.  
type = 'pre';
channels_to_load = [10:13];
extract_varargin;

%% set the values for a given input type
if strcmp(type,'pre')
    fname = [fname '-pre.dat'];
    t_end = 200;
elseif strcmp(type,'post')
   fname = [fname '-post.dat'];
     t_end = 200;
elseif strcmp(type,'main')
  t_end = 300;
  fname = [fname '.dat'];
end
       
%% load wideband data to spikesort; can be a single channel, or up to 4 channels
% using at least two is a good idea to separate noise from spikes!
% data = AMPX_loadData(fname,channels_to_load);

%% optional: get correlation matrix of decimated data to determine appropriate referencing
data_ref = AMPX_loadData(fname,1:54,20); % note this will be slow -- it is loading a LOT of data!

corrMat = nan(54);

for ii = 1:54
    for jj = ii:54
   
        cc = corrcoef(data_ref.channels{ii},data_ref.channels{jj});
        corrMat(ii,jj) = cc(1,2);
        
    end
end  
  
ref = figure;
imagesc(corrMat);
saveas(ref,['D:\DATA\R036\Ref_Matricies\' fname(1:end-4)],'png')
close(ref)
%% do the re-referencing // note this should be based on sensible channels
% % tic
% % av = AMPX_getAverage(fname,[1:8 10:15 17 19:23 26:30 32]); % also slow!
% % data = AMPX_reref(data,av);
% % toc

%% filter for spikes, and place in matrix appropriate for UMS2K
% clear D;
% invert_waveforms = 0; % if you have spikes bigger in the positive direction
% 
% tic
% disp('filter')
% for ii = length(data.channels):-1:1
%     D(ii,:) = filter_for_spikes(data.channels{ii}, data.hdr.Fs);
% end
% if invert_waveforms, D = -D; end
% toc
% 
% %% take a look
% % t_end = 200; % number of seconds to plot
% 
% figure(1)
% clf
% 
% plot(D(:,1:t_end*data.hdr.Fs)')