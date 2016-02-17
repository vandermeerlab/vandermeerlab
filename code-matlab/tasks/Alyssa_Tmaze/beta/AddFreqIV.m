function IV = AddFreqIV(cfg,IV,CSC)
%ADDFREQIV Add usr data containing information about frequency to iv data
%
%   IV = ADDFREQIV(cfg,IV,CSC) computes the Fast Fourier Transform on
%   segments of LFP and adds information about peak frequency or peak
%   amplitude to iv usr data.
%   Filtering the signal in the frequency range you are interested in is
%   recommended. For example, if you are interested in ripples you should
%   exclude low frequencies since they often have the highest amplitude in
%   unfiltered HC signals.
%
%    INPUTS
%       cfg: config struct with fields controlling function behaviour
%        IV: interval data to add usr field to
%       CSC: local field potential data as recorded by a continuously
%            sampled channel (csc)
%
%    OUTPUT
%        IV: interval data with usr field containing frequency information
%  
%    CONFIG OPTIONS
%       cfg.output = 'fpeak'; 
%             'fpeak' - Return the frequency that has the highest amplitude
%             'fampl' - Return the highest amplitude (the max Fourier
%                       coefficient)
%       cfg.label = cfg.output; The usr field name. If you do not provide a
%                       label, the default is the string in cfg.output.
%
%       cfg.verbose = 1; If 1, prints informative text to the command
%             window; if 0, doesn't
%
% aacarey Feb 2016

cfg_def.output = 'fpeak'; % 'fpeak','fampl'
cfg_def.label = [];
cfg_def.verbose = 1;

mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg,mfun);

if ~CheckIV(IV,mfun)
   error('IV is poorly formed. See the iv data type constructor iv()') 
end

if ~CheckTSD(CSC)
   error('CSC is poorly formed. See the tsd data type constructor tsd()') 
end

if cfg.verbose; fprintf('%s: Adding frequency data to IV.usr...\n',mfun); end

% If user does not specify an alternative label, call it the same as cfg.output
if isempty(cfg.label); cfg.label = cfg.output; end

for iInterval = length(IV.tstart):-1:1
    % restrict LFP to the current interval
   CSCr = restrict(CSC,IV.tstart(iInterval),IV.tend(iInterval));
   % compute fft and keep magnitude
   fftedData = abs(fft(CSCr.data));
   % keep the first half
   fftedData = fftedData(1:round(length(fftedData)/2));
   % get peak
   [maxPeak,idx] = max(fftedData(:));
   
   switch cfg.output
       case 'fpeak'
           % return peak frequency
           timeWin = IV.tend(iInterval) - IV.tstart(iInterval);
           IV.usr.(cfg.label)(iInterval) = idx/timeWin;
           
       case 'fampl'
           % return peak fourier coeff
           IV.usr.(cfg.label)(iInterval) = maxPeak;
           
       otherwise
           error('Unrecognized option specified in cfg.output')
   end
   
end

% write config history
IV = History(IV,mfun,cfg);

end

