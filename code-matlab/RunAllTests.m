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
    
    % clean
    if ispc
        system('del *.xml');
    elseif isunix
        system('rm *.xml');
    end
    
    % set up path containing all folders to be tracked for code coverage
    src = fullfile(pwd, 'shared');
    p = genpath(src); p_tok = regexp(p, pathsep, 'split'); % tokenize path so we can exclude folders later
    
    addpath(p);
    
    % remove non-tracked folders from path
    rmpath(genpath(fullfile(pwd, 'shared\io\open_ephys')));
    rmpath(genpath(fullfile(pwd, 'shared\io\neuralynx')));
    rmpath(genpath(fullfile(pwd, 'shared\viz\matlab-plot-big')));
    rmpath(genpath(fullfile(pwd, 'shared\viz\export_fig')));
    rmpath(genpath(fullfile(pwd, 'shared\viz\linspecer')));
    rmpath(genpath(fullfile(pwd, 'shared\viz\numSubplots')));
    
    fp = path; fp_tok = regexp(fp, pathsep, 'split');
    
    tracked_folders = intersect(p_tok, fp_tok);
    
    % path contining all tests to be run
    tests = fullfile(pwd, 'tests');
    suite = testsuite(tests, 'IncludeSubfolders', true);
    
    runner = TestRunner.withTextOutput();
    
    % add TAP
    tapFile = fullfile(getenv('WORKSPACE'), 'testResults.tap');
    runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(tapFile)));
    
    % Add Cobertura: need to add each tracked folder separately, apparently
    fprintf('\nRunAllTests.m: Adding folders to cover:\n')
    for iF = 1:length(tracked_folders)
        this_folder = tracked_folders{iF};
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