function ExpKeys = ParLoadExpKeys
%PARLOADEXPKEYS Load ExpKeys into base workspace 'transparently'
%for use in parfor loops 
%   No need for any inputs; it checks the current directory automatically 
%
% ACarey May 2015 

%%

%[~,name,~] = fileparts(pwd); 

fn = FindFiles('*keys.m');
if isempty(fn)
    ExpKeys = [];
    disp(['ParLoadExpKeys: No files matching',' ''*keys.m'' ','were found in', [' ',pwd]])
    disp('ExpKeys not loaded')
elseif length(fn) > 1
    ExpKeys = [];
    disp(['ParLoadExpKeys: More than one file matching',' ''*keys.m'' ','was found in', [' ',pwd]])
    disp('ExpKeys not loaded')
else
    ExpKeys = [];
    run(fn{1});
end

end

