function nS = BarSpikes(S,varargin)

% BarSpikes(S)
% 
% input: S -- cell array of spike times (class ts)
% plots bar graph of total spikes fired for each cell in S
% also works with non-class ts

% ADR 
% version L4.0
% status: PROMOTED

display = 0;
Extract_varargin;

if (~isa(S, 'cell'))
  error('Expects cell array as input.');
end

v = zeros(size(S));
for iS=1:length(S)
   if (isa(S{iS}, 'ts'))
      v(iS) = length(Data(S{iS}));
   elseif (isa(S{iS}, 'tsd'))
      v(iS) = length(Data(S{iS}));      
   else
      v(iS) = length(S{iS});
   end   
end
if display
	bar(v);
end

nS = v;