function passflag = TEST_restrict(cfg_in)
% function passflag = TEST_restrict(cfg_in)
%
% unit tests for restrict()
%
% MvdM 2015-11-30 basic initial version
% youkitan 2016-10-15 include edge tests
% youkitan 2019-01-27 update test readability

%% parse cfg

cfg_def = [];
cfg_def.verbose = 1;
cfg_def.visuals = 0;
cfg_def.version = 'restrict3';

if nargin > 1
    error('TEST_restrict only takes a config file')
elseif nargin < 1
    cfg_in = [];
    disp('WARNING: there is no input config file!');
end

cfg = ProcessConfig(cfg_def,cfg_in);

passflag = 1; % set to 1 by default, make 0 if any test fails

% preliminaries: report which version of restrict() we are using
fprintf('\n*** UNIT TESTS FOR RESTRICT ***\n');

fp = which(cfg.version);
if isempty(fp)
	error(sprintf(['No version of restrict exists with the name ''',cfg.version,'''']))
else
    fprintf('\nTesting %s\n\n',fp);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 1: TESTS FOR PRESERVING INPUT FORMAT CORRECTLY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% test TS input

% construct some testing data - TS
empty_ts = ts;

one_ts_data = (0:1:10)';
one_ts = ts;
one_ts.t{1} = one_ts_data;
one_ts.label{1} = 'label1';

rvec = [0 10]; curr_iv = iv(rvec);

%%%%%%%%%%%%% test empty TS, bounding value arguments %%%%%%%%%%%%%
switch cfg.version
    case 'restrict'
        empty_ts_out = restrict(empty_ts,rvec(1),rvec(2));
    case 'restrict2'
        empty_ts_out = restrict2(empty_ts,rvec(1),rvec(2));
    case 'restrict3'        
        empty_ts_out = restrict3(empty_ts,rvec(1),rvec(2));

end

empty_ts_out_flag = CheckTS(empty_ts_out);
switch empty_ts_out_flag
    case 1
        fprintf('PASSED: empty TS, raw restrict argument\n');
    case 0
        fprintf('FAILED: empty TS, raw restrict argument\n');
        passflag = 0;
end

%%%%%%%%%%%%% test empty TS, iv restrict argument %%%%%%%%%%%%%
switch cfg.version
    case 'restrict'
        empty_ts_out = restrict(empty_ts,curr_iv);
    case 'restrict2'
        empty_ts_out = restrict2(empty_ts,curr_iv);
    case 'restrict3'
        empty_ts_out = restrict3(empty_ts,curr_iv);
end

empty_ts_out_flag = CheckTS(empty_ts_out);
switch empty_ts_out_flag
    case 1
        fprintf('PASSED: empty TS, iv restrict argument\n');
    case 0
        fprintf('FAILED: empty TS, iv restrict argument\n');
        passflag = 0;
end

%%%%%%%%%%%%% test one t TS, iv restrict argument %%%%%%%%%%%%%
switch cfg.version
    case 'restrict'
        one_ts_out = restrict(one_ts,curr_iv);
    case 'restrict2'
        one_ts_out = restrict2(one_ts,curr_iv);
    case 'restrict3'
        one_ts_out = restrict3(one_ts,curr_iv);
end
switch CheckTS(one_ts_out)
    case 1
        fprintf('PASSED: one t TS, iv restrict argument\n');
    case 0
        fprintf('FAILED: one t TS, iv restrict argument\n');
        passflag = 0;
end

%% test TSD input

% construct some testing data - TSD

%% test IV input

% construct some testing data - IV


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 2: TESTS FOR CORRECT PROCESSING OF EDGE CASES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test edge cases for TS input
curr_iv = iv(rvec);

one_ts_out_old = restrict_old(one_ts,curr_iv);
switch tscmp(one_ts_out,one_ts_out_old)
    case 1
        fprintf('PASSED: one t TS, iv argument, same as old restrict \n');
    case 0
        fprintf('FAILED: one t TS, iv argument, not same as old restrict \n');
        passflag = 0;
end

%%%%%%%%%%%%% one t TS, ceiling effect %%%%%%%%%%%%%
under_edge = [0.999999999999 1.999999999999]; curr_iv = iv(under_edge);
switch cfg.version
    case 'restrict'
        under_edge_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        under_edge_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        under_edge_ts = restrict3(one_ts,curr_iv);
end
switch vertcat(under_edge_ts.t{:}) == 1
    case 1
        fprintf('PASSED: one t TS, iv no ceil effect \n');
    case 0
        fprintf('FAILED: one t TS, iv ceil effect \n');
        passflag = 0;
end

%%%%%%%%%%%%% one t TS, floor effect %%%%%%%%%%%%%
over_edge = [1.000000000001 2.000000000001]; curr_iv = iv(over_edge);
switch cfg.version
    case 'restrict'
        over_edge_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        over_edge_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        over_edge_ts = restrict3(one_ts,curr_iv);
end
switch vertcat(over_edge_ts.t{:}) == 2
    case 1
        fprintf('PASSED: one t TS, iv no floor effect \n');
    case 0
        fprintf('FAILED: one t TS, iv floor effect \n');
        passflag = 0;
end

%%%%%%%%%%%%% one t TS, iv wraps outside edges %%%%%%%%%%%%%
sandwich_buns = [.999999999999 2.000000000001]; curr_iv = iv(sandwich_buns);
switch cfg.version
    case 'restrict'
        sandwich_buns_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        sandwich_buns_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        sandwich_buns_ts = restrict3(one_ts,curr_iv);
end
if isequal(vertcat(sandwich_buns_ts.t{:}),[1;2])
    fprintf('PASSED: one t TS, iv wrap outside \n');
else
    fprintf('FAILED: one t TS, iv wrap outside \n');
    passflag = 0;
end

%%%%%%%%%%%%% one t TS, iv wraps inside edges %%%%%%%%%%%%%
sandwich_meat = [1.000000000001 1.999999999999]; curr_iv = iv(sandwich_meat);
switch cfg.version
    case 'restrict'
        sandwich_meat_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        sandwich_meat_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        sandwich_meat_ts = restrict3(one_ts,curr_iv);
end
if isempty(vertcat(sandwich_meat_ts.t{:})) 
    fprintf('PASSED: one t TS, iv wrap inside \n');
else
    fprintf('FAILED: one t TS, iv wrap inside \n');
    passflag = 0;
end

%%%%%%%%%%%%% one t TS, iv goes across edge %%%%%%%%%%%%%
across_edge = .87:.05:1.02;
across_edge2 = .97:.05:1.12;
curr_iv = iv(across_edge,across_edge2);
if cfg.visuals; PlotIV([],curr_iv); end
switch cfg.version
    case 'restrict'
        across_edge_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        across_edge_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        across_edge_ts = restrict3(one_ts,curr_iv);
end
if isequal(vertcat(across_edge_ts.t{:}),1)
    fprintf('PASSED: one t TS, iv across edge \n');
else
    fprintf('FAILED: one t TS, iv across edge \n');
    passflag = 0;
end

%%%%%%%%%%%%% one t TS, iv sits on edge %%%%%%%%%%%%%
on_edge = .9:0.1:1;
on_edge2 = 1:0.1:1.1;
curr_iv = iv(on_edge,on_edge2);
switch cfg.version
    case 'restrict'
        on_edge_ts = restrict(one_ts,curr_iv);
    case 'restrict2'
        on_edge_ts = restrict2(one_ts,curr_iv);
    case 'restrict3'
        on_edge_ts = restrict3(one_ts,curr_iv);
end
if cfg.visuals; PlotIV([],curr_iv); end
if isequal(vertcat(on_edge_ts.t{:}),1)
    fprintf('PASSED: one t TS, iv on edge \n');
else
    fprintf('FAILED: one t TS, iv on edge \n');
    passflag = 0;
end

%%%%%%%%%%%%% two t TS, compare to original %%%%%%%%%%%%%
multiple_ts_data = (0:1:2)';
two_ts = ts;
two_ts.t{1} = multiple_ts_data; two_ts.label{1} = 'label1';
two_ts.t{2} = multiple_ts_data; two_ts.label{2} = 'label2';

all_edges = [under_edge; over_edge; sandwich_buns; sandwich_meat; [across_edge' across_edge2']; [on_edge' on_edge2']];
curr_iv = iv(all_edges);
s1 = restrict_old(two_ts,curr_iv);
switch cfg.version
    case 'restrict'
        s2 = restrict(two_ts,curr_iv);
    case 'restrict2'
        s2 = restrict2(two_ts,curr_iv);
    case 'restrict3'
        s2 = restrict3(two_ts,curr_iv);
end
if tscmp(s1,s2)
    fprintf('PASSED: multiple t TS, same result as previous generation \n')
else
    fprintf('FAILED: multiple t TS, same result as previous generation \n')
    passflag = 0;
end
    
%% Test edge cases for TSD input


%% Test edge cases for IV input


%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 2: SPEED TESTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end