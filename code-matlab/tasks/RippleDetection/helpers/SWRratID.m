function rat = SWRratID(sessionID)
%SWRRATID Identify rat ID from session ID
%   rat = SWRratID(sessionID)
%
%   INPUTS
%       [~,sessionID,~] = fileparts(pwd) when current directory is set to
%           the folder for a given recording day.
%   OUTPUTS
%       rat - a string specifying the rat's ID (ex: 'R064' or 'bon')
%
%   NOTE: this function uses SWRrats(). To add a new rat, open SWRrats and
%   edit the rat list.
%
% aacarey Dec 2017

rats = SWRrats;

rat = {};

for iRat = 1:length(rats)
    if strfind(sessionID,rats{iRat})
        rat = rats{iRat};
    end
end

if isempty(rat)
   error('Unable to identify rat') 
end
end

