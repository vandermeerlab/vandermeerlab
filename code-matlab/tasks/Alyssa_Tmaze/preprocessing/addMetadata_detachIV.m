% get iv for detachment events R044
% amplitude thresholding. expand intervals.
clear
threshold = 0.0015; % the amplitude in volts, above which the data will be excluded
cfg.d = [-0.025 0.025]; % ms, amount to expand detachIV

cfg.rats = {'R044'};
fd = sort(getTmazeDataPath(cfg)); % get all session directories 

for iFD = 1:length(fd)
    
    cd(fd{iFD});
    LoadExpKeys
    cfg_temp.fc = ExpKeys.goodSWR(1);
    csc = LoadCSC(cfg_temp);
    
    tstart = csc.tvec(diff(abs(csc.data) > threshold) == 1);
    tend = csc.tvec(diff(abs(csc.data) > threshold) == -1);
    detachIV = iv(tstart,tend);
    detachIV = ResizeIV(cfg,detachIV);
    detachIV = MergeSingleIV([],detachIV);
    
    % make it work with restrict(), so the iv has the stuff we want to keep
    tstart = [csc.tvec(1); detachIV.tend];
    tend = [detachIV.tstart; csc.tvec(end)];
    detachIV = iv(tstart,tend);

    csc2 = restrict(csc,detachIV);
   % figure; hold on
   % plot(csc.tvec,csc.data);
   % plot(csc2.tvec,csc2.data,'r');
    
    % Save fields in metadata
       % first check if metadata exists yet
    
    loaded = LoadMetadata2;
    
    if ~loaded % then it doesn't exist yet, so make a metadata struct
        metadata.detachIV = detachIV;
        metadata.detachIV = detachIV;
    else % it does, so add a new field
        metadata.detachIV = detachIV;
        metadata.detachIV = detachIV;
    end
    
    [~,name,~] = fileparts(pwd); % pwd is your current folder, we just want its namepart
    
       % now save
    savename = strcat(name,'-metadata.mat'); % use the folder's name, but concatenate it with '-metadata'
    save(savename,'metadata'); % this saves the specified variables under the given [save]name
    
    clearvars -except cfg fd iFD threshold
    
end