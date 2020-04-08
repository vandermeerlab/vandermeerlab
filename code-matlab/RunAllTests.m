import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TAPPlugin;
import matlab.unittest.plugins.ToFile;
import('matlab.unittest.plugins.CodeCoveragePlugin');
import('matlab.unittest.plugins.codecoverage.CoberturaFormat');

% pathsep
if ispc
    pathsep = ';';
elseif isunix
    pathsep = ':';
else
   error ('Undefined path separator.');
end

% find out hostname
[~, hostname] = system('hostname'); hostname = hostname(1:end - 1);

try
    % make sure we start from where RunAllTests.m is located
    whoami = mfilename('fullpath');
    whereami = fileparts(whoami);
    cd(whereami);
    
    % set up path containing all folders to be tracked for code coverage
    src = fullfile(ws, 'shared');
    p = genpath(src);
    addpath(p); % this returns a single string with ; as separator between folders
    
    % path contining all tests to be run
    tests = fullfile(ws, 'tests');
    suite = testsuite(tests, 'IncludeSubfolders', true);
    
    runner = TestRunner.withTextOutput();
    
    % add TAP
    tapFile = fullfile(getenv('WORKSPACE'), 'testResults.tap');
    runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(tapFile)));
    
    % Add Cobertura: need to add each tracked folder separately, apparently
    fprintf('\nRunAllTests.m: Adding folders to cover:\n')
    sep = strfind(p,pathsep); idx = 1; % idx "cursor" tracks position in path string
    for iF = 1:length(sep)
        this_folder = p(idx:sep(iF) - 1);
        disp(this_folder);
        
        coverageFile = fullfile(getenv('WORKSPACE'), sprintf('coverage%d.xml',iF));
        runner.addPlugin(CodeCoveragePlugin.forFolder(this_folder,'Producing', CoberturaFormat(coverageFile)));
        
        idx = sep(iF)+1; % update cursor to start of next path
    end

    % Run the tests
    results = runner.run(suite);
    display(results);
catch e
    fprintf('\n*********************\nRunAllTests.m failed!\n*********************\n');
    disp(getReport(e,'extended'));
    

    if strcmp(hostname,'mvdmlab-athena'), exit(1); end % hack to only quit on CI machine
end
% if running on CI machine, exit
if strcmp(hostname,'mvdmlab-athena'), exit; end 