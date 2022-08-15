function S = LoadSpikes(cfg_in)
% LOADSPIKES Load .t files or .tt files containing spike timestamps 
%
%   S = LOADSPIKES(cfg) Loads spiketrains into a ts struct. Spikes can
%   originate as .t files (spike sorted in MClust) or .tt files (all spikes 
%   recorded by a tetrode, see MakeTTFiles). To create an S struct from
%   spiketrains with another format (such as from another lab), load the
%   files and use the ts datatype constructor.
%
%   CFG OPTIONS:
%
%       cfg.fc = {}; Cell array containing filenames to load. If 
%                no filenames are specified, loads all *.t files in the
%                current directory.
%       cfg.load_questionable_cells = 0; Load *._t files if set to 1
%       cfg.getTTnumbers = 1; If 1, includes the tetrode number each 
%                spiketrain originated from in S.usr.tt_num. If 0 there is 
%                no such field.
%       cfg.getRatings = 0; If 1, includes cluster ratings in S.usr.rating. 
%                If 0, (or if .tt files are loaded) there is no such field.
%       cfg.min_cluster_quality = 5; Minimum cluster quality to load for .t
%                files (applies if cfg.getRatings = 1).
%       cfg.tsflag = 'sec'; Units to use for the timestamps
%                 'sec'  -  seconds
%                  'ms'  -  milliseconds
%                  'ts'  -  exact (as exists in the file)
%       cfg.verbose = 1; If 1, allow command window text. If 0, suppress
%                command window text.
%
%   OUTPUT
%
%       S: spike data ts struct
%
% see also ts, MakeTTFiles
%
% MvdM 2014-06-17 based on ADR LoadSpikes(), 25 edit to use cfg_in
% aacarey edit Nov 2015
% youkitan edit Aug 2016 fixed getRatings/min_cluster_quality issue

cfg_def.fc = {};
cfg_def.tsflag = 'sec';
cfg_def.load_questionable_cells = 0;
cfg_def.min_cluster_quality = 5;
cfg_def.getRatings = 0;
cfg_def.getTTnumbers = 1;
cfg_def.verbose = 1;
cfg_def.uint = '32';

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % this takes fields from cfg_in and puts them into cfg

if ~isempty(cfg.fc)
    % can't get ratings for .tt files, also can't load questionable cells
    ext_tt = false(size(cfg.fc));
    for iFN = 1:length(ext_tt)
        [~,~,ext] = fileparts(cfg.fc{iFN});
        if strcmp(ext,'.tt')
            ext_tt(iFN) = true;
        end
    end
    
    if any(ext_tt)
        % reset cfg.getRatings to 0 both so that it won't attempt to get the
        % ratings and so that the config history shows that ratings were not
        % added to usr field
        if cfg.verbose && cfg.getRatings; fprintf('WARNING in %s: at least one file is a .tt file, ratings cannot be included in usr\n',mfun); end
        if cfg.verbose && cfg.load_questionable_cells; fprintf('WARNING in %s: at least one file is a .tt file, cfg.load_questionable_cells does not apply\n',mfun); end
        cfg.getRatings = 0;
        cfg.load_questionable_cells = 0;
        cfg.min_cluster_quality = [];
    end
end

if isempty(cfg.fc)
        
        cfg.fc = FindFiles('*.t');
        
        if cfg.load_questionable_cells
            if cfg.verbose; fprintf('%s: WARNING: loading questionable cells\n',mfun); end
           cfg.fc = cat(1,cfg.fc,FindFiles('*._t'));
        end
        
else
        
        if ~isa(cfg.fc,'cell')
            error('LoadSpikes: cfg.fc should be a cell array.');
        end
        
end

cfg.fc = sort(cfg.fc);
nFiles = length(cfg.fc);

if cfg.verbose; fprintf('%s: Loading %d files...\n',mfun,nFiles); end

% for each tfile
% first read the header, then read a tfile

S = ts; % initialize new ts struct

for iF = 1:nFiles
    tfn = cfg.fc{iF};
    
    if ~isempty(tfn)
        
        tfp = fopen(tfn, 'rb','b');
        if (tfp == -1)
            error(['LoadSpikes: Could not open tfile ' tfn]);
        end
        
        ReadHeader(tfp);
        switch cfg.uint
            case '32'
                S.t{iF} = fread(tfp,inf,'uint32'); % read as 32 bit ints
            case '64'
                S.t{iF} = fread(tfp,inf,'uint64');	% read as 64 bit ints
        end
        
        % set appropriate time units
        switch cfg.tsflag
            case 'sec'
                S.t{iF} = S.t{iF}/10000;
            case 'ts'
                S.t{iF} = S.t{iF};
            case 'ms'
                S.t{iF} = S.t{iF}*10000*1000;
            otherwise
                error('LoadSpikes: invalid tsflag.');
        end
        
        % add filenames
        [~,fname,fe] = fileparts(tfn);
        S.label{iF} = cat(2,fname,fe);
        
        fclose(tfp);
        
    end 		% if tfn valid
end		% for all files

% this makes sure that if you specify a min cluster quality it will select appropriate clusters
% otherwise you can have specified a min cluster quality without using the getRatings flag
% and not have LoadSpikes restrict the cluster quality
if cfg.min_cluster_quality < 5 && cfg.getRatings == 0
    cfg.getRatings = 1;
end
    
if cfg.getRatings
        
        warning('off','MATLAB:unknownElementsNowStruc');
        
        for iC = 1:length(S.label)
            
            curr_fn = S.label{iC};
            
            % get tt name and cell number of this neuron
            tok = regexp(curr_fn,'(.*)\_([0-9]+)','tokens');
            
            % load clusters
            clu_fn = cat(2,tok{1}{1},'.clusters');
            
            if exist(clu_fn,'file')
                load(clu_fn,'-mat');
            else
                S.usr.rating(iC) = nan;
                warning('File %s does not exist. Cannot get rating.',clu_fn);
                continue
            end
            
            % get rating
            cellno = str2double(tok{1}{2});
            
            if exist('MClust_Clusters','var')
                clu_name = MClust_Clusters{cellno}.name;
            elseif exist('Clusters','var')
                clu_name = Clusters{cellno}.name;
            else
                error('what clusters??'); % snarky LoadSpikes is so appalled
            end
            
            if ischar(clu_name(1))
                S.usr.rating(iC) = str2double(clu_name(1)); % first character of name should be rating
            else
                S.usr.rating(iC) = nan;
                warning('Cluster %d has no rating (name is %s). Cannot get rating.',cellno,clu_name);               
            end
            
        end % of filenames to process
        
        warning('on','MATLAB:unknownElementsNowStruc');
        
        % select only cells that meet criterion
        
        keep_idx = find(S.usr.rating <= cfg.min_cluster_quality);
        
        S.t = S.t(keep_idx);
        S.label = S.label(keep_idx);
        S.usr.rating = S.usr.rating(keep_idx);
        
        if cfg.verbose; fprintf('%s: Outputting %d spiketrains that meet criteria...\n',mfun,length(S.t)); end 
end

% check if ExpKeys available
keys_f = FindFiles('*keys.m');
if ~isempty(keys_f)
    try
    run(keys_f{1});
    catch
       disp('Failed to load keys file.'); 
    end
    S.cfg.ExpKeys = ExpKeys;
end

% add tt numbers
if cfg.getTTnumbers
    for iC = length(S.t):-1:1
        
        out = regexp(S.label{iC},'TT(\d\d)','tokens');
        
        if isempty(out)
            warning('Filename %s does not contain TTxx, cannot extract tetrode number.',S.label{iC});
        else
            S.usr.tt_num(iC) = str2double(cell2mat(out{1}));
        end
    end
end

% add sessionID
[~,S.cfg.SessionID,~] = fileparts(pwd);

% housekeeping
S = History(S,mfun,cfg);