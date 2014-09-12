function OK = WriteTfiles(self)

% Write T files (MClustData)

MCS = MClust.GetSettings();
MCD = self;

% ADR 2013-12-12
if isempty(MCD.Clusters)
    msgbox('There are no clusters to write.', 'WriteTFiles','warn');
    OK = false;
    return 
end

% some of the clusters may be _t
if MCS.UseUnderscoreT
    names = cellfun(@(x)x.name, MCD.Clusters, 'UniformOutput', false);
    underscoreTclusters = listdlg(...
        'ListString', names, ...
        'Name', 'Files to save as _t', ...
        'PromptString', 'Which clusters should be saved as with an "_t" extension?', ...
        'OKString', 'DONE', 'CancelString', 'No _t files.', ...
        'InitialValue', []);    
else
    underscoreTclusters = [];
end

nClust = length(MCD.Clusters);

% ADR 2008 - delete .t set if replacing
fcT = FindFiles([MCD.TTfn '_*.t'], 'StartingDirectory', MCD.TTdn, 'CheckSubdirs', 0);
fc_T = FindFiles([MCD.TTfn '_*._t'], 'StartingDirectory', MCD.TTdn, 'CheckSubdirs', 0);
fc = cat(1, fcT, fc_T);

if ~isempty(fc)
	reply = questdlg({'There are already .t or ._t files for this tetrode.','Do you want to replace them?'},...
		'Overwrite t files?', 'Yes', 'Cancel', 'Yes');
	if streq(reply, 'Yes')
		for iC = 1:length(fc)
			delete(fc{iC});
		end
	else
		OK = false;
		return
	end
end
for iC = 1:nClust
   spikes = MCD.Clusters{iC}.GetSpikes();
   if ~isempty(spikes)
       
      tSpikes = MCD.FeatureTimestamps(spikes);
      
      if ismember(iC, underscoreTclusters)
          fn = [MCD.TfileBaseName(iC) '._' MCS.tEXT];
      else
          fn = [MCD.TfileBaseName(iC) '.' MCS.tEXT];
      end            
      
      fp = fopen(fn, 'wb', 'b');
      if (fp == -1)
         errordlg(['Could not open file"' fn '".']);
      end
      MClust.WriteHeader(fp, ...
          'T-file', ...
          'Output from MClust', ...
          'Time of spiking stored in timestamps (tenths of msecs)',...
          'as unsigned integer: uint64');
      switch MCS.tEXT
          case 't64'
              tSpikes = uint64(tSpikes*10000); % NEED TO CONVERT TO NEURALYNX's .t format save in integers of 0.1ms
              fwrite(fp, tSpikes, 'uint64');
          case 't32'
              tSpikes = uint64(tSpikes*10000); % NEED TO CONVERT TO NEURALYNX's .t format save in integers of 0.1ms
              fwrite(fp, tSpikes, 'uint64');
          case 't'
              tSpikes = uint64(tSpikes*10000); % NEED TO CONVERT TO NEURALYNX's .t format save in integers of 0.1ms
              fwrite(fp, tSpikes, 'uint64');
          otherwise
              error('MClust::tEXT', 'Unknown extension for t files');
      end
       
      fclose(fp);
   end
end
OK = true;