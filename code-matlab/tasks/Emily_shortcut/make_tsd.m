function [ pos_tsd,tvec ] = make_tsd( pos_x,pos_y,Timestamps )

% Make tsd
pos_tsd = tsd;

tvec = Timestamps * 10^-6;
pos_tsd.tvec = tvec;

pos_tsd.data(1,:) = pos_x;
pos_tsd.label{1} = 'x';

pos_tsd.data(2,:) = pos_y;
pos_tsd.label{2} = 'y';

% check if ExpKeys available
keys_f = FindFiles('*keys.m');
if ~isempty(keys_f)
    run(keys_f{1});
    pos_tsd.cfg.ExpKeys = ExpKeys;
end

mfun = mfilename;
cfg = [];

% add sessionID
[~,pos_tsd.cfg.SessionID,~] = fileparts(pwd);

% housekeeping
pos_tsd.cfg.history.mfun = cat(1,pos_tsd.cfg.history.mfun,mfun);
pos_tsd.cfg.history.cfg = cat(1,pos_tsd.cfg.history.cfg,{cfg});

end