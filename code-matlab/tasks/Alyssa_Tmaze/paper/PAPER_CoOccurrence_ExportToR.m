%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                Export CoOccurrence Data to R                         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tmaze project
%
% This script loads co-occurrence data produced by PAPER_CoOccurrence and
% organizes it in a format that is useable by R (for statistical
% analysis).
%
% Ouput is saved to the same directory as the input
%
%   data_out.p
%       [n x 1] double containing all the p values specified by cfg.whichP.
%       The p values are the output from CoOccurQ.
%   data_out.arm
%       [n x 1] double containing the number 1 for p values belonging to cell
%       pairs on the left arm, and the number 2 for cell pairs on the right
%       arm.
%   data_out.restriction
%       [n x 1] double containing the number 1 for p values from food days,
%       and the number 2 for p values from water days
%       arm.
%   data_out.rat
%       [n x 1] double containing the numbers 1-4 identifying which rat the
%       p value belongs to: 1 - R042, 2 - R044, 3 - R050, 4 - R064.
%
% aacarey Jan 2015 from code by MvdM

%clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Where is the data kept?
fd = [pwd,'\data']; % E:\Documents\TmazePaper\data

% Which p value do you want to collect?
cfg.whichP = 'p4'; % Options: 'p0' (proportion SWR single cell),'p4' (zscore cooccurrence),'p5' (shuffled data), also 'p1','p2','p3'

% Which event category do you want to collect data for?
cfg.whichEvents = 'prerecord'; % Options:  'all', 'prerecord', 'task', 'postrecord', 'allITI', 'equalBehaviorITI'


%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.cc);
end

%% initialize some things

rats = TmazeRats;
restr = {'food','water'};
arms = {'L','R'};

if ~any(strcmp(cfg.whichEvents,{'all', 'prerecord', 'task', 'postrecord', 'allITI', 'equalBehaviorITI'}))
    error('cfg.whichEvents contains an unrecognized event category')
end

if ~any(strcmp(cfg.whichP,{'p0','p1','p2','p3','p4','p5'}))
    error('cfg.whichP contains an unrecognized p option')
end

% assemble filenames
fn = ['coOccurrence_',cfg.whichEvents];
fn_out = [fn,'_',cfg.whichP,'_R'];

% load the data structure
originalFolder = pwd;
cd(fd);
data_in = loadpop(fn,'coocData');

%%
data_out.p = [];
data_out.arm = [];
data_out.restriction = [];
data_out.rat = [];

for iRat = 1:length(rats)
    
    % go through each session
    for iSession = 1:length(data_in.(rats{iRat}))
        
        % get restriction type and assign an id for it
        switch data_in.(rats{iRat})(iSession).restrictionType
            case 'food'
                restr = 1;
            case 'water'
                restr = 2;
        end
        
        for iArm = 1:length(arms)
            
            % get the target p field
            p = data_in.(rats{iRat})(iSession).(arms{iArm}).ALLp.(cfg.whichP);
            
            % remove NaNs
            p = p(~isnan(p));
            
            % grow the output structure
            data_out.p = [data_out.p; p];
            data_out.arm = [data_out.arm; repmat(iArm,size(p))];
            data_out.restriction = [data_out.restriction; repmat(restr,size(p))];
            data_out.rat = [data_out.rat; repmat(iRat,size(p))];
            
        end % of arms
    end % of sessions   
end % of rats

%% remove NaNs
% n = isnan(data_out.p);
% fprintf('%d NaNs removed.\n',sum(n));
% data_out.p = data_out.p(~n);
% data_out.arm = data_out.arm(~n);
% data_out.restriction = data_out.restriction(~n);
% data_out.rat = data_out.rat(~n);

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

save(fn_out,'data_out');
cd(originalFolder)

disp(' ')
disp('~~~~~~~ END OF COOCCURRENCE EXPORT ~~~~~~~')